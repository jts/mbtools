//---------------------------------------------------------
// Copyright 2021 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
extern crate clap;
use std::time::{Instant};
//use fnv::{FnvHashMap, FnvHashSet};
use hashbrown::{HashMap, HashSet};
use clap::{Arg, App, SubCommand, value_t};
use rust_htslib::{bam, bam::Read, bam::record::Aux, bam::record::Cigar::*};

//
pub struct AlignedPair
{
    reference_index: usize,
    read_index: usize
}

fn calculate_aligned_pairs(record: &bam::Record) -> Vec::<AlignedPair> {
    let mut aligned_pairs = Vec::<AlignedPair>::new();

    let mut current_read_index :usize = 0;
    let mut current_reference_index :usize = record.pos() as usize;

    for c in record.cigar().iter() {

        let (aligned, reference_stride, read_stride) = match *c {
            Match(_) | Equal(_) | Diff(_) => (true, 1, 1),
            Del(_) => (false, 1, 0),
            Ins(_) => (false, 0, 1),
            SoftClip(_) => (false, 0, 1),
            HardClip(_) | Pad(_) => (false, 0, 0),
            RefSkip(_) => (false, 1, 0)
        };

        for _i in 0 .. c.len() {
            if aligned {
                aligned_pairs.push( AlignedPair{reference_index:current_reference_index, read_index:current_read_index} );
            }

            current_reference_index += reference_stride;
            current_read_index += read_stride;
        }
    }
    aligned_pairs
}

// struct storing the modifications 
pub struct ReadModifications
{
    canonical_base: char,
    modified_base: char,
    modification_indices: Vec<usize>,
    modification_probabilities: Vec<f64>,
    read_to_reference_map: HashMap<usize, usize>
}

impl ReadModifications
{
    pub fn from_bam_record(record: &bam::Record) -> Self {
        //println!("Parsing bam record");
        let mut rm = ReadModifications {
            canonical_base: 'x',
            modified_base: 'x',
            modification_indices: vec![],
            modification_probabilities: vec![],
            read_to_reference_map: HashMap::new()
        };

        // TODO: remove when we switch bam_seq to match original strand of read
        assert_eq!(record.is_reverse(), false);

        //
        let bam_seq = record.seq().as_bytes();
        let mut modified_seq = bam_seq.clone();
        //let qname = std::str::from_utf8(record.qname()).unwrap();

        // parse probabilities
        if let Ok(Aux::ArrayU8(array)) = record.aux(b"Ml") {
            for encoded_probability in array.iter() {
                rm.modification_probabilities.push(encoded_probability as f64 / 255.0);
            }
        }

        if let Ok(Aux::String(mm_str)) = record.aux(b"Mm") {

            rm.canonical_base = mm_str.as_bytes()[0] as char;
            rm.modified_base = mm_str.as_bytes()[2] as char;

            // calculate the index in the read of each canonical base
            let mut canonical_indices = Vec::<usize>::new();
            for (index, base) in bam_seq.iter().enumerate() {
                if *base as char == rm.canonical_base {
                    canonical_indices.push(index)
                }
            }

            // parse the modification string and transform it into indices in the read
            let mut canonical_count : usize = 0;
            if mm_str.len() > 4 {
                for token in mm_str[4..].split(",") {
                    canonical_count += token.parse::<usize>().unwrap();
                    rm.modification_indices.push(canonical_indices[canonical_count]);
                    canonical_count += 1;
                }
            }

            // temporary set of reference positions with a modification
            let mut read_modification_set = HashSet::<usize>::new();
            for i in &rm.modification_indices {
                read_modification_set.insert(*i);
            }

            // extract the alignment from the bam record
            let aligned_pairs = calculate_aligned_pairs(record);

            rm.read_to_reference_map.reserve(rm.modification_indices.len());
            for t in &aligned_pairs {
                if read_modification_set.contains(&t.read_index) {
                    rm.read_to_reference_map.insert(t.read_index, t.reference_index);
                }
            }

            // construct modified sequence (debug)
            for index in &rm.modification_indices {
                modified_seq[*index] = rm.modified_base as u8;
            }
        }

        //println!("Parsed {} -> {} modifications for {}\n{}", rm.canonical_base, rm.modified_base, qname, std::str::from_utf8(&modified_seq).unwrap());
        rm
    }
}

fn main() {
    let matches = App::new("mbtools")
        .version("0.1")
        .author("Jared Simpson <jared.simpson@oicr.on.ca>")
        .about("Toolkit for working with modification bam files")
        .subcommand(SubCommand::with_name("modification-frequency")
                .about("calculate the frequency of modified bases per position of the genome")
                .arg(Arg::with_name("probability-threshold")
                    .short("t")
                    .long("probability-threshold")
                    .takes_value(true)
                    .help("only use calls where the probability of being modified/not modified is at least t"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process")))
        .get_matches();


    if let Some(matches) = matches.subcommand_matches("modification-frequency") {

        // TODO: set to nanopolish default LLR
        let threshold = value_t!(matches, "probability-threshold", f64).unwrap_or(0.8);
        calculate_modification_frequency(threshold, matches.value_of("input-bam").unwrap())
    }
}

fn calculate_modification_frequency(threshold: f64, input_bam: &str) {
    println!("calculating modification frequency with t:{} on file {}", threshold, input_bam);

    let mut bam = bam::Reader::from_path(input_bam).unwrap();
    //let header = bam::Header::from_template(bam.header());

    // map from (tid, position) -> (methylated_reads, total_reads)
    let mut reference_modifications = HashMap::<(i32, usize), (usize, usize)>::new();

    //
    let start = Instant::now();
    let mut reads_processed = 0;
    for r in bam.records() {
        let record = r.unwrap();
        let tid = record.tid();
        let rm = ReadModifications::from_bam_record(&record);
        //for (mod_index, mod_probability) in zip(rm.modification_indices.iter(), rm.modification_probabilities.iter()) {
        for (mod_index, mod_probability) in rm.modification_indices.iter().zip(rm.modification_probabilities.iter()) {

            if (*mod_probability > threshold) || ((1.0 - mod_probability) > threshold) {
                let reference_position = rm.read_to_reference_map.get(mod_index).expect("Unable to find reference position for read base");
                let mut e = reference_modifications.entry( (tid, *reference_position) ).or_insert( (0, 0) );
                if *mod_probability > 0.5 {
                    (*e).0 += 1;
                }
                (*e).1 += 1;
                //println!("processing {} i:{} p:{} ({} {})", reference_position, mod_index, mod_probability, e.0, e.1);
            }
        }

        reads_processed += 1;
    }

    //
    for ( (tid, position), (methylated_reads, total_reads) ) in reference_modifications {
        println!("{}\t{}\t{}\t{}\t{}", tid, position, methylated_reads, total_reads, methylated_reads as f64 / total_reads as f64);
    }
    let duration = start.elapsed();
    eprintln!("Processed {} reads in {:?}", reads_processed, duration);
}
