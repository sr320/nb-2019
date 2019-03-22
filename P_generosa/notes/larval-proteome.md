

A TruSeq Stranded mRNA library (Illumina) was created using geoduck larvae (5 dpf) and sequenced on NovaSeq Platform (2 lanes). Sequence reads were assembled using Trinity.

```
Trinity \
--seqType fq \
--SS_lib_type RF \
--left NR021_S8_L001_R1_001.fastq,NR021_S8_L002_R1_001.fastq \
--right NR021_S8_L001_R2_001.fastq,NR021_S8_L002_R2_001.fastq \
--CPU 50 --trimmomatic --max_memory 500G
```

[output](http://owl.fish.washington.edu/halfshell/bu-mox/analyses/0804_1818/trinity_out_dir/0804_Pgen_larvae.fasta)
