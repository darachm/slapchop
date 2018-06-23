
# **S**imply **L**ooking **A**t **P**airwise **C**omparisons for **O**ptimized **P**arsing

ie **SLAPCHOP**

This is a program to take a big fastq file, align it to an expected 
pattern of information, then make a filtering decision based on this 
alignment, and use this matching to extract out various elements from
the read (like an interspersed UMI from an undetermined location).

Goals for the development:

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


