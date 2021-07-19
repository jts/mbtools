mbtools - toolkit for working with modification BAM files
---------------------------------------------------------

Current tools:

```
# calculate the modification frequency at every position of the reference genome
# any calls with P(modified | data) < THRESHOLD or P(not modified | data) < THRESHOLD
# will be ignored
mbtools modification-frequency -t THRESHOLD input.bam > output.tsv
```
