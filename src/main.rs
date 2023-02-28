//---------------------------------------------------------
// Copyright 2021 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
extern crate clap;
use std::time::{Instant};
use hashbrown::{HashMap};
use itertools::Itertools;
use clap::{Arg, App, SubCommand, value_t};
use rust_htslib::{bam, faidx, bam::Read, bam::record::Aux, bam::record::Cigar::*, bam::ext::BamRecordExtensions};
use bio::alphabets;
use intervaltree::IntervalTree;
use core::ops::Range;
use std::cmp::{min, max};

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

    let cigar = record.cigar();
    if cigar.leading_hardclips() > 0 || cigar.trailing_hardclips() > 0 {
        eprintln!("Error: cannot process alignments with hardclips");
        std::process::exit(1);
    }

    for c in cigar.iter() {

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

// fill in call indices/modification probabilities
// so it has an entry for every position in canonical indices
fn fill_untagged_bases(canonical_indices: &Vec<usize>,
                       read_to_reference_map: HashMap::<usize, usize>,
                       calls: &mut Vec<ModificationCall>)
{
    let mut curr_mod_index = 0;
    let tmp_calls = calls.clone();

    calls.clear();

    for i in 0 .. canonical_indices.len() {

        let mut call = ModificationCall {
            read_index: canonical_indices[i],
            modification_probability: 0.0,
            reference_index: read_to_reference_map.get(&canonical_indices[i]).cloned()
        };

        if curr_mod_index < tmp_calls.len() && canonical_indices[i] == tmp_calls[curr_mod_index].read_index {
            // copy probability
            call.modification_probability = tmp_calls[curr_mod_index].modification_probability;
            curr_mod_index += 1;
        } 
        calls.push(call);
    }

    assert_eq!(curr_mod_index, tmp_calls.len());
}

// 
#[derive(Copy, Clone)]
pub struct ModificationCall
{
    read_index: usize,
    modification_probability: f64,
    reference_index: Option<usize>
}

impl ModificationCall
{
    pub fn is_modified(self) -> bool {
        return self.modification_probability > 0.5;
    }

    pub fn get_probability_correct(self) -> f64 {
        let p = if self.is_modified() { self.modification_probability } else { 1.0 - self.modification_probability };
        return p
    }

    pub fn is_confident(self, threshold: f64) -> bool {
        self.get_probability_correct() > threshold
    }
}

// struct storing the modifications 
pub struct ReadModifications
{
    canonical_base: char,
    modified_base: char,
    strand: char,
    modification_calls: Vec<ModificationCall>
}

// this wraps uncertainty around whether the draft tag (Mm) is used
pub fn get_mm_tag(record: &bam::Record) -> Result<Aux, rust_htslib::errors::Error> {
    let r = match record.aux(b"MM") {
        Ok(v) => Ok(v),
        Err(_) => record.aux(b"Mm")
    };
    return r;
}

pub fn get_ml_tag(record: &bam::Record) -> Result<Aux, rust_htslib::errors::Error> {
    let r = match record.aux(b"ML") {
        Ok(v) => Ok(v),
        Err(_) => record.aux(b"Ml")
    };
    return r;
}

impl ReadModifications
{

    pub fn from_bam_record(record: &bam::Record) -> Option<Self> {
        
        // records that are missing the SEQ field cannot be processed
        if record.seq().len() == 0 {
            return None;
        }

        let mut rm = ReadModifications {
            canonical_base: 'x',
            modified_base: 'x',
            strand: '+',
            modification_calls: vec![]
        };

        // this SEQ is in the same orientation as the reference,
        // revcomp it here to make it the original orientation that the instrument read
        let mut instrument_read_seq = record.seq().as_bytes();
        if record.is_reverse() {
            instrument_read_seq = alphabets::dna::revcomp(instrument_read_seq);
            rm.strand = '-';
        }

        // do I need to nest these?
        if let Ok(Aux::String(mm_str)) = get_mm_tag(record) {
            if let Ok(Aux::ArrayU8(probability_array)) = get_ml_tag(record) {

                // TODO: handle multiple mods
                if mm_str.matches(';').count() != 1 {
                    return None;
                }

                let first_mod_str = mm_str.split(';').next().unwrap();
                let mod_meta = first_mod_str.split(',').next().unwrap().as_bytes();

                rm.canonical_base = mod_meta[0] as char;
                rm.modified_base = mod_meta[2] as char;
                let mut assume_canonical = true;
                let optional_flag_idx = mod_meta.len() - 1;
                if mod_meta[optional_flag_idx] as char == '?' {
                    assume_canonical = false;
                }

                // calculate the index of each canonical base in the read
                // ACTACATA -> (1,4) if C is canonical
                let mut canonical_indices = Vec::<usize>::new();
                for (index, base) in instrument_read_seq.iter().enumerate() {
                    if *base as char == rm.canonical_base {
                        canonical_indices.push(index);
                    }
                }
                
                // make a map from read positions to where they map on the reference genome
                let mut read_to_reference_map = HashMap::<usize, usize>::new();

                // aligned_pairs stores the index of the read along SEQ, which may be
                // reverse complemented. switch the coordinates here to be the original
                // sequencing direction to be consistent with the Mm/Ml tag.
                let mut aligned_pairs = calculate_aligned_pairs(record);
                let seq_len = instrument_read_seq.len();
                if record.is_reverse() {
                    for t in &mut aligned_pairs {
                        t.read_index = seq_len - t.read_index - 1;
                    }
                }

                for t in &aligned_pairs {
                    read_to_reference_map.insert(t.read_index, t.reference_index);
                }

                // parse the modification string and transform it into indices in the read
                let mut canonical_count : usize = 0;

                for (token, encoded_probability) in first_mod_str.split(',').skip(1).zip(probability_array.iter()) {
                    canonical_count += token.parse::<usize>().unwrap();
                    let i = canonical_indices[canonical_count];
                    
                    let call = ModificationCall {
                        read_index: i,
                        modification_probability: encoded_probability as f64 / 255.0,
                        reference_index: read_to_reference_map.get(&i).cloned()
                    };

                    rm.modification_calls.push(call);
                    canonical_count += 1;
                }

                if assume_canonical {
                    fill_untagged_bases(&mut canonical_indices, read_to_reference_map, &mut rm.modification_calls);
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
        .subcommand(SubCommand::with_name("reference-frequency")
                .about("calculate the frequency of modified bases per position of the genome")
                .arg(Arg::with_name("collapse-strands")
                    .short("c")
                    .long("collapse-strands")
                    .takes_value(false)
                    .help("merge the calls for the forward and negative strand into a single value"))
                .arg(Arg::with_name("probability-threshold")
                    .short("t")
                    .long("probability-threshold")
                    .takes_value(true)
                    .help("only use calls where the probability of being modified/not modified is at least t"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process")))
        .subcommand(SubCommand::with_name("read-frequency")
                .about("calculate the frequency of modified bases per read")
                .arg(Arg::with_name("probability-threshold")
                    .short("t")
                    .long("probability-threshold")
                    .takes_value(true)
                    .help("only use calls where the probability of being modified/not modified is at least t"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process")))
        .subcommand(SubCommand::with_name("read-region-frequency")
                .about("calculate the frequency of modified bases per read filtered across regions within the BED file")
                .arg(Arg::with_name("probability-threshold")
                    .short("t")
                    .long("probability-threshold")
                    .takes_value(true)
                    .help("only use calls where the probability of being modified/not modified is at least t"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process"))
                .arg(Arg::with_name("region-bed")
                    .short("r")
                    .long("region-bed")
                    .required(true)
                    .takes_value(true)
                    .help("bed file containing the regions to calculate modification frequencies for")))
        .subcommand(SubCommand::with_name("region-frequency")
                .about("calculate the frequency of modified bases for regions provided within the BED file")
                .arg(Arg::with_name("probability-threshold")
                    .short("t")
                    .long("probability-threshold")
                    .takes_value(true)
                    .help("only use calls where the probability of being modified/not modified is at least t"))
                .arg(Arg::with_name("reference-genome")
                    .short("g")
                    .long("reference-genome")
                    .takes_value(true)
                    .required(false)
                    .help("path to reference genome"))
                .arg(Arg::with_name("cpg")
                    .long("cpg")
                    .takes_value(false)
                    .required(false)
                    .help("only include calls at CpG dinucleotides"))
                .arg(Arg::with_name("region-bed")
                    .short("r")
                    .long("region-bed")
                    .takes_value(true)
                    .required(true)
                    .help("bed file containing the regions to calculate modification frequencies for"))
                .arg(Arg::with_name("classification-threshold")
                    .short("c")
                    .long("classification-threshold")
                    .takes_value(true)
                    .help("threshold to classify reads as U/X/M using Loyfer criteria"))
                .arg(Arg::with_name("classification-min-sites")
                    .short("m")
                    .long("classification-min-sites")
                    .takes_value(true)
                    .help("minimum number of sites needed to classify a read"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process")))
        .get_matches();

    // TODO: set to nanopolish default LLR
    let calling_threshold = value_t!(matches, "probability-threshold", f64).unwrap_or(0.8);

    if let Some(matches) = matches.subcommand_matches("reference-frequency") {

        calculate_reference_frequency(calling_threshold,
                                      matches.is_present("collapse-strands"),
                                      matches.value_of("input-bam").unwrap())
    }
    
    if let Some(matches) = matches.subcommand_matches("read-region-frequency") {

        // TODO: set to nanopolish default LLR
        calculate_read_region_frequency(calling_threshold,
                                    matches.value_of("region-bed").unwrap(),
                                    matches.value_of("input-bam").unwrap())
    }
    if let Some(matches) = matches.subcommand_matches("read-frequency") {

        // TODO: set to nanopolish default LLR
        calculate_read_frequency(calling_threshold,
                                    matches.value_of("input-bam").unwrap())
    }

    if let Some(matches) = matches.subcommand_matches("region-frequency") {

        // TODO: set to nanopolish default LLR
        let classification_threshold = value_t!(matches, "classification-threshold", f64).unwrap_or(0.25);
        let classification_min_sites = value_t!(matches, "classification-min-sites", usize).unwrap_or(4);
        calculate_region_frequency(calling_threshold,
                                   classification_threshold,
                                   classification_min_sites,
                                   matches.value_of("region-bed").unwrap(),
                                   matches.value_of("input-bam").unwrap(),
                                   matches.is_present("cpg"),
                                   matches.value_of("reference-genome").unwrap_or(""))
    }
}

fn calculate_reference_frequency(threshold: f64, collapse_strands: bool, input_bam: &str) {
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

        if let Some(rm) = ReadModifications::from_bam_record(&record) {
            
            for call in rm.modification_calls {
                if call.is_confident(threshold) && call.reference_index.is_some() {
                    let mut reference_position = call.reference_index.unwrap().clone();
                    let mut strand = rm.strand;
                    if collapse_strands && strand == '-' {
                        reference_position -= 1; // TODO: this only works for CpG
                        strand = '+';
                    }
                    let mut e = reference_modifications.entry( (record.tid(), reference_position, strand) ).or_insert( (0, 0) );
                    (*e).0 += call.is_modified() as usize;
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

fn calculate_read_frequency(threshold: f64, input_bam: &str) {

    let mut bam = bam::Reader::from_path(input_bam).expect("Could not read input bam file:");
    let header_view = bam.header().clone();
    let start = Instant::now();
    let mut reads_processed = 0;
    println!("read_name\tchromosome\tstart_position\tend_position\talignment_length\tstrand\tmapping_quality\ttotal_calls\tmodified_calls\tmodification_frequency");
    let mut summary_total = 0;
    let mut summary_modified = 0;

    eprintln!("calculating read modifications with t:{} on file {}", threshold, input_bam);
    for r in bam.records() {
        let record = r.unwrap();

        if record.is_unmapped() {
            continue;
        }

        let mut total_calls = 0;
        let mut total_modified = 0;

        if let Some(rm) = ReadModifications::from_bam_record(&record) {

            for call in rm.modification_calls {
                if call.is_confident(threshold) {
                    total_calls += 1;
                    total_modified += call.is_modified() as usize;
                }
            }
        }

        let qname = std::str::from_utf8(record.qname()).unwrap();
        let strand = if record.is_reverse() { '-' } else { '+' };
        let mod_frequency = if total_calls > 0 { total_modified as f64 / total_calls as f64 } else { f64::NAN };
        let contig = String::from_utf8_lossy(header_view.tid2name(record.tid() as u32));
        let start_position = record.pos();
        let end_position = record.reference_end();
        let alignment_length = end_position - start_position + 1;
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}", 
            qname, contig, record.pos(), end_position, alignment_length, 
            strand, record.mapq(), total_calls, total_modified, mod_frequency);
    
        // for the ending summary line
        reads_processed += 1;
        summary_total += total_calls;
        summary_modified += total_modified;
    }
    let summary_frequency = if summary_total > 0 { summary_modified as f64 / summary_total as f64 } else { f64::NAN };
    eprintln!("Processed {} reads in {:?}. Mean modification frequency: {:.2}", reads_processed, start.elapsed(), summary_frequency);
}
struct RegionStats {
    chromosome: String,
    start: usize,
    end: usize,
    mod_calls: usize,
    total_calls: usize,
    called_reads: usize,
    total_prob: f64,
    u_reads: usize,
    m_reads: usize,
    x_reads: usize
}
fn interval_tree_from_bed(region_bed: &str, header_view: bam::HeaderView) -> (HashMap<u32, IntervalTree<usize, usize>>, Vec<RegionStats>) {
    // Read bed file into a data structure we can use to make intervaltrees from
    // this maps from tid to a vector of intervals, with an interval index for each
    let mut bed_reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(region_bed).expect("could not open bed file");
    let mut region_desc_by_chr = HashMap::<u32, Vec<(Range<usize>, usize)>>::new();

    // this stores the modified/total counts for each interval
    let mut region_data = Vec::new();

    for r in bed_reader.records() {
        let record = r.expect("could not parse bed record");
        if let Some(tid) = header_view.tid(record[0].as_bytes()) {
            let start: usize = record[1].parse().unwrap();
            let end: usize = record[2].parse().unwrap();
            let region_desc = region_desc_by_chr.entry(tid).or_insert( Vec::new() );
            region_desc.push( (start..end, region_data.len()) );
            let init = RegionStats {
                chromosome: record[0].to_string(),
                start: start,
                end: end,
                mod_calls: 0,
                total_calls: 0,
                called_reads: 0,
                total_prob: 0.0,
                u_reads: 0,
                m_reads: 0,
                x_reads: 0 
            };
            region_data.push( init );
        }
    }
    // build tid -> intervaltree map
    // the intervaltree allows us to look up an interval_idx for a given chromosome and position
    let mut interval_trees = HashMap::<u32, IntervalTree<usize, usize>>::new();
    for (tid, region_desc) in region_desc_by_chr {
        interval_trees.insert(tid, region_desc.iter().cloned().collect());
    }
    (interval_trees, region_data)
}

fn calculate_read_region_frequency(threshold: f64, region_bed: &str, input_bam: &str) {

    eprintln!("calculating read modifications with t:{} on file {}\n Modifications filtered to those within the regions of {}", threshold, input_bam, region_bed);
    println!("read_name\tchromosome\tstart_position\tend_position\talignment_length\tstrand\tmapping_quality\ttotal_calls\tmodified_calls\tmodification_frequency");
    let mut bam = bam::Reader::from_path(input_bam).expect("Could not read input bam file:");
    let header_view = bam.header().clone();
    let start = Instant::now();
    let mut reads_processed = 0;
    let mut summary_total = 0;
    let mut summary_modified = 0;
    let (interval_trees, _) = interval_tree_from_bed(region_bed, header_view.clone());
        
    // iterate over bam records
    for r in bam.records() {
        let record = r.unwrap();
        if record.is_unmapped() {
            continue;
        }
        // Check if the read overlaps any of the regions
        // for every region the read overlaps, we write a record for that region
        // counting only the modifications calls in that region
        let tid = record.tid() as u32;
        if let Some(tree) = interval_trees.get(&tid) {
            let start = record.pos() as usize;
            let end = record.reference_end() as usize;
            if let Some(rm) = ReadModifications::from_bam_record(&record) {
                for region in tree.query(start..end) {
                     let mut total_calls = 0;
                     let mut total_modified = 0;
                     let start_position = max(start, region.clone().range.start);
                     let end_position = min(end, region.clone().range.end);
                     // get the read modifications for this read
                     // for each modification call, check if it is in the region
                     for call in &rm.modification_calls {
                        if call.is_confident(threshold) && call.reference_index.is_some() {
                            let reference_position = call.reference_index.unwrap().clone();
                            if reference_position >= start_position && reference_position <= end_position {
                                total_calls += 1;
                                total_modified += call.is_modified() as usize;
                            }
                        }
                    }
                    let qname = std::str::from_utf8(record.qname()).unwrap();
                    let strand = if record.is_reverse() { '-' } else { '+' };
                    let mod_frequency = if total_calls > 0 { total_modified as f64 / total_calls as f64 } else { f64::NAN };
                    let contig = String::from_utf8_lossy(header_view.tid2name(record.tid() as u32));
                    let alignment_length = record.reference_end() - record.pos() + 1;
                    println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}", 
                        qname, contig, start_position, end_position, alignment_length, 
                        strand, record.mapq(), total_calls, total_modified, mod_frequency);
                    summary_total += total_calls;
                    summary_modified += total_modified;
                }
            }
        }
        // for the ending summary line
        reads_processed += 1;
    }
    let summary_frequency = if summary_total > 0 { summary_modified as f64 / summary_total as f64 } else { f64::NAN };
    eprintln!("Processed {} reads in {:?}. Mean modification frequency: {:.2}", reads_processed, start.elapsed(), summary_frequency);
}

fn calculate_region_frequency(call_threshold: f64, 
                              classification_threshold: f64, 
                              classification_min_sites: usize, 
                              region_bed: &str, 
                              input_bam: &str, 
                              filter_to_cpg: bool, 
                              reference_genome: &str) {
    eprintln!("calculating modification frequency for regions from {} on file {}", region_bed, input_bam);
    let faidx = match filter_to_cpg {
        true => Some(faidx::Reader::from_path(reference_genome).expect("Could not read reference genome:")),
        false => None
    };
    
    let mut bam = bam::Reader::from_path(input_bam).expect("Could not read input bam file:");
    let header = bam::Header::from_template(bam.header());
    let header_view = bam::HeaderView::from_header(&header);

    let (interval_trees, mut region_data) = interval_tree_from_bed(region_bed, header_view.clone());

    let mut curr_chromosome_id = -1;
    let mut curr_chromosome_seq = String::new();
    let mut curr_chromosome_length = 0;

    //
    let start = Instant::now();
    let mut reads_processed = 0;
    for r in bam.records() {
        let record = r.unwrap();
        if record.is_unmapped() {
            continue
        }

        // if this record is on a chromosome we don't have in memory, load it
        if filter_to_cpg && record.tid() != curr_chromosome_id {
            let contig = String::from_utf8_lossy(header_view.tid2name(record.tid() as u32));
            curr_chromosome_length = header_view.target_len(record.tid() as u32).unwrap() as usize;
            
            curr_chromosome_id = record.tid();
            curr_chromosome_seq = faidx.as_ref().expect("faidx not found").fetch_seq_string(contig, 0, curr_chromosome_length).unwrap();
            curr_chromosome_seq.make_ascii_uppercase();
            
        }

        if let Some(rm) = ReadModifications::from_bam_record(&record) {

            // keep track of read-level stats per interval
            let mut read_stats_by_interval = HashMap::new();
                     
            for call in rm.modification_calls {
                if call.is_confident(call_threshold) && call.reference_index.is_some() {
                    let reference_position = call.reference_index.unwrap().clone();
                    
                    // For the cpg filter
                    let mut motif_position = reference_position;
                    if rm.strand == '-' && reference_position > 0 {
                        motif_position -= 1;
                    }

                    // 
                    if filter_to_cpg && 
                        ( rm.strand == '+' || reference_position > 0 ) &&
                        motif_position < curr_chromosome_length - 1 && 
                        &curr_chromosome_seq[motif_position..motif_position+2] != "CG" {
                        continue
                    }
                    //eprintln!("sub: {}", reference_subseq);

                    if let Some(tree) = interval_trees.get( &(record.tid() as u32)) {
                        for element in tree.query_point(reference_position) {
                            let interval_idx = element.value;
                            
                            let is = &mut read_stats_by_interval.entry(interval_idx).or_insert( (0, 0, 0.0) );
                            is.0 += call.is_modified() as usize;
                            is.1 += 1;
                            is.2 += call.get_probability_correct();
                        }
                    }
                }
            }
            
            for (idx, read_data) in read_stats_by_interval {
                let data = &mut region_data[idx];
                data.mod_calls += read_data.0;
                data.total_calls += read_data.1;
                data.total_prob += read_data.2;

                data.called_reads += (read_data.1 > 0) as usize;
            
                if read_data.1 >= classification_min_sites {
                    let f = read_data.0 as f64 / read_data.1 as f64;
                    if f < classification_threshold {
                        data.u_reads += 1;
                    } else if f > (1.0 - classification_threshold)  {
                        data.m_reads += 1;
                    } else {
                        data.x_reads += 1;   
                    }
                }
            }
        }

        reads_processed += 1;
    }

    println!("chromosome\tstart\tend\tnum_called_reads\tmodified_calls\ttotal_calls\tmodification_frequency\tmean_probability_correct\tu_reads\tm_reads\tx_reads");
    for d in region_data {
        let f = d.mod_calls as f64 / d.total_calls as f64;
        let p = d.total_prob as f64 / d.total_calls as f64;
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{}\t{}\t{}", 
            d.chromosome, d.start, d.end, d.called_reads, d.mod_calls, d.total_calls, f, p, d.u_reads, d.m_reads, d.x_reads);
    }
    eprintln!("Processed {} reads in {:?}", reads_processed, start.elapsed());
}

