
# Simply Looking At Pairwise Comparisons for Optimized Parsing

ie SLAPCHOP

This is a program to take a big fastq file, align it to expected subsequences
in the reads, then make a decision of filtration based on this as to keep
the reads or now, then also to use this matching to extract out an
interspersed UMI from an undetermined location in the read, found by alignment.

Currently just runs for SoBaSeq amplicons (low-input barseq).

Goals for the future:

- Tests, to make sure it's doing as we expect it to do
- A better understanding of controlling queuesize and total memory footprint
    to make that super reliable for tuning to the max memory.
- Different levels of report summary generation, ie per-read statistics or
    just distributions, or nothing
- Report generation, so making summary plots from these reports
    (seperate R script??? something in python???)
- Generalize the template design so you can configure via arguments of
    a template sequence, a regular expression, and the ability to chain these.
- Come up with a funny bit, describing possible applications of the tool
    using the lingo and style of the commercial

An example:

    python3 slapChop.py --processes 26 --biteSize 10000 \
      --maxQueueSize 50 --inputFastq ${i} \
      --outputFastq ${outbase}.chopped.fastq \
      --fixed1pattern GTCCACGAGGTCTCT 
      --fixed2pattern CGTACGCTGCAGGTCGAC \
      --umipattern CGTACGCTGCAGGTCGACXGXAXGXGXGXGAT \
      --outputReport ${outbase}.report --logFile ${outbase}.log \
      --filters "AlignmentStart1==5 and Score1 > 26 and Score2 > 28" 


