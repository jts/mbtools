//---------------------------------------------------------
// Copyright 2021 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
extern crate clap;
use std::time::{Instant};
use hashbrown::{HashMap, HashSet};
use itertools::Itertools;
use clap::{Arg, App, SubCommand, value_t};
use rust_htslib::{bam, bam::Read, bam::record::Aux, bam::record::Cigar::*};
use bio::alphabets;

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

// fill in modification indices/modification probabilities
// so it has an entry for every position in canonical indices
fn fill_untagged_bases(canonical_indices: &Vec<usize>,
                       modification_indices: &mut Vec<usize>,
                       modification_probabilities: &mut Vec<f64>)
{
    let mut curr_mod_index = 0;
    let tmp_indices = modification_indices.clone();
    let tmp_probabilities = modification_probabilities.clone();

    modification_indices.clear();
    modification_probabilities.clear();

    for i in 0 .. canonical_indices.len() {

        modification_indices.push(canonical_indices[i]);
        if curr_mod_index < tmp_indices.len() && canonical_indices[i] == tmp_indices[curr_mod_index] {
            // copy probability from the input vector
            modification_probabilities.push(tmp_probabilities[curr_mod_index]);
            curr_mod_index += 1;
        } else {
            modification_probabilities.push(0.0);
        }
    }

    assert_eq!(curr_mod_index, tmp_indices.len());
}

// struct storing the modifications 
pub struct ReadModifications
{
    canonical_base: char,
    modified_base: char,
    strand: char,
    modification_indices: Vec<usize>,
    modification_probabilities: Vec<f64>,
    read_to_reference_map: HashMap<usize, usize>
}

impl ReadModifications
{
    pub fn from_bam_record(record: &bam::Record, assume_canonical: bool) -> Option<Self> {
        
        // records that are missing the SEQ field cannot be processed
        if record.seq().len() == 0 {
            return None;
        }

        //println!("Parsing bam record");
        let mut rm = ReadModifications {
            canonical_base: 'x',
            modified_base: 'x',
            strand: '+',
            modification_indices: vec![],
            modification_probabilities: vec![],
            read_to_reference_map: HashMap::new()
        };

        // this SEQ is in the same orientation as the reference,
        // revcomp it here to make it the original orientation that the instrument read
        let mut instrument_read_seq = record.seq().as_bytes();
        if record.is_reverse() {
            instrument_read_seq = alphabets::dna::revcomp(instrument_read_seq);
            rm.strand = '-';
        }

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
            for (index, base) in instrument_read_seq.iter().enumerate() {
                if *base as char == rm.canonical_base {
                    canonical_indices.push(index)
                }
            }

            // parse the modification string and transform it into indices in the read
            let mut canonical_count : usize = 0;
            assert_eq!(mm_str.matches(';').count(), 1);
            for token in mm_str.split(';').next().unwrap().split(',').skip(1) {
                canonical_count += token.parse::<usize>().unwrap();
                rm.modification_indices.push(canonical_indices[canonical_count]);
                canonical_count += 1;
            }

            if assume_canonical {
                fill_untagged_bases(&mut canonical_indices, &mut rm.modification_indices, &mut rm.modification_probabilities);
            }

            // extract the alignment from the bam record
            let mut aligned_pairs = calculate_aligned_pairs(record);

            // aligned_pairs stores the index of the read along SEQ, which may be
            // reverse complemented. switch the coordinates here to be the original
            // sequencing direction to be consistent with the Mm/Ml tag.
            let seq_len = instrument_read_seq.len();
            if record.is_reverse() {
                for t in &mut aligned_pairs {
                    t.read_index = seq_len - t.read_index - 1;
                }
            }

            // temporary set of read positions with a modification
            let mut read_modification_set = HashSet::<usize>::new();
            for i in &rm.modification_indices {
                read_modification_set.insert(*i);
            }

            rm.read_to_reference_map.reserve(rm.modification_indices.len());
            for t in &aligned_pairs {
                if read_modification_set.contains(&t.read_index) {
                    rm.read_to_reference_map.insert(t.read_index, t.reference_index);
                }
            }
        }

        //println!("Parsed {} -> {} modifications for {}\n{}", rm.canonical_base, rm.modified_base, qname, std::str::from_utf8(&modified_seq).unwrap());
        Some(rm)
    }
}

fn main() {
    let matches = App::new("mbtools")
        .version("0.1")
        .author("Jared Simpson <jared.simpson@oicr.on.ca>")
        .about("Toolkit for working with modification bam files")
        .subcommand(SubCommand::with_name("modification-frequency")
                .about("calculate the frequency of modified bases per position of the genome")
                .arg(Arg::with_name("collapse-strands")
                    .short("c")
                    .long("collapse-strands")
                    .takes_value(false)
                    .help("merge the calls for the forward and negative strand into a single value"))
                .arg(Arg::with_name("assume-canonical")
                    .short("a")
                    .long("assume-canonical")
                    .takes_value(false)
                    .help("assume bases not present in the Mm tag are canonical (unmodified)"))
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
        calculate_modification_frequency(threshold,
                                         matches.is_present("collapse-strands"),
                                         matches.is_present("assume-canonical"),
                                         matches.value_of("input-bam").unwrap())
    }
}

fn calculate_modification_frequency(threshold: f64, collapse_strands: bool, assume_canonical: bool, input_bam: &str) {
    eprintln!("calculating modification frequency with t:{} on file {}", threshold, input_bam);

    let mut bam = bam::Reader::from_path(input_bam).expect("Could not read input bam file:");
    let header = bam::Header::from_template(bam.header());

    // map from (tid, position, strand) -> (methylated_reads, total_reads)
    let mut reference_modifications = HashMap::<(i32, usize, char), (usize, usize)>::new();

    //
    let start = Instant::now();
    let mut reads_processed = 0;
    for r in bam.records() {
        let record = r.unwrap();

        if let Some(rm) = ReadModifications::from_bam_record(&record, assume_canonical) {
            
            for (mod_index, mod_probability) in rm.modification_indices.iter().zip(rm.modification_probabilities.iter()) {
                let is_modified_call = *mod_probability > 0.5;
                let probability_correct = if is_modified_call { *mod_probability } else { 1.0 - *mod_probability };
                let map_lookup = rm.read_to_reference_map.get(mod_index);
                if probability_correct > threshold && map_lookup.is_some() {
                    let mut reference_position = map_lookup.unwrap().clone();
                    let mut strand = rm.strand;
                    if collapse_strands && strand == '-' {
                        reference_position -= 1; // TODO: this only works for CpG
                        strand = '+';
                    }
                    let mut e = reference_modifications.entry( (record.tid(), reference_position, strand) ).or_insert( (0, 0) );
                    (*e).0 += is_modified_call as usize;
                    (*e).1 += 1;
                }
            }
        }

        reads_processed += 1;
    }

    //
    let mut sum_reads = 0;
    let mut sum_modified = 0;

    let header_view = bam::HeaderView::from_header(&header);
    println!("chromosome\tposition\tstrand\tmodified_reads\ttotal_reads\tmodified_frequency");
    for key in reference_modifications.keys().sorted() {
        let (tid, position, strand) = key;
        let contig = String::from_utf8_lossy(header_view.tid2name(*tid as u32));
        let (modified_reads, total_reads) = reference_modifications.get( key ).unwrap();
        println!("{}\t{}\t{}\t{}\t{}\t{:.3}", contig, position, strand, modified_reads, total_reads, *modified_reads as f64 / *total_reads as f64);
    
        sum_reads += total_reads;
        sum_modified += modified_reads;
    }
    
    let mean_depth = sum_reads as f64 / reference_modifications.keys().len() as f64;
    let mean_frequency = sum_modified as f64 / sum_reads as f64;
    eprintln!("Processed {} reads in {:?}. Mean depth: {:.2} mean modification frequency: {:.2}", reads_processed, start.elapsed(), mean_depth, mean_frequency);
}
