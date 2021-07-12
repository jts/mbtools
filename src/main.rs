extern crate clap;
use clap::{Arg, App, SubCommand, value_t};

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
    println!("calculating modification frequency with t:{} on file {}", threshold, input_bam)
}
