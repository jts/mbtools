mbtools - toolkit for working with modification BAM files
---------------------------------------------------------

## Compiling

Building mbtools requires rust, which can be installed following [these instructions](https://doc.rust-lang.org/cargo/getting-started/installation.html)

Then you can clone and build this program:

```
git clone https://github.com/jts/mbtools
cd mbtools
cargo build --release
```

The compiled program can be found in `target/release/mbtools`

## Programs

Current tools:

```
# calculate the modification frequency at every position of the reference genome
mbtools reference-frequency input.bam > output.tsv

# calculate the modification frequency for each read
mbtools read-frequency input.bam > output.tsv

# calculate the modification frequency for each region in a bed file
mbtools region-frequency -r regions.bed input.bam > output.tsv

# calculate the modification frequency for each read segmented by regions in a bed file (good for long reads spanning multiple regions)
mbtools region-read-frequency -r regions.bed input.bam > output.tsv
```
