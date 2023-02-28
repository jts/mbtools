mbtools - toolkit for working with modification BAM files
---------------------------------------------------------

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
