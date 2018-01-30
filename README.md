
# Simply Looking At Pairwise Comparisons for Optimized Parsing

ie SLAPCHOP

This is a program to take a big fastq file, align it to expected subsequences
in the reads, then make a decision of filtration based on this as to keep
the reads or now, then also to use this matching to extract out an
interspersed UMI from an undetermined location in the read, found by alignment.

Currently just runs for SoBaSeq amplicons (low-input barseq of the yeast deletion collection).

Goals for the future:

- Generalize the template design so you can configure via arguments of
    a template sequence, a regular expression, and the ability to chain these internally.
- Report generation, so making summary plots from these reports
- Different levels of report summary generation, ie per-read statistics or
    just distributions, or nothing
- Unit tests on choice examples from the sequencing.
- A better understanding of controlling queuesize and total memory footprint
    to make that super reliable for tuning to the max memory.

An example:

    python3 slapChop.py --processes 26 --biteSize 10000 \
      --maxQueueSize 50 --inputFastq input.fastq \
      --outputFastq output.chopped.fastq \
      --fixed1pattern GTCCACGAGGTCTCT 
      --fixed2pattern CGTACGCTGCAGGTCGAC \
      --umipattern CGTACGCTGCAGGTCGACXGXAXGXGXGXGAT \
      --outputReport output.report --logFile output.log \
      --filters "AlignmentStart1==5 and Score1 > 26 and Score2 > 28" 


