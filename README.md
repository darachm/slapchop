
# URE2PRION - Using Regular Expressions 2 Parse Reads Into Other New Sequences

__(formerly SLAPCHOP)__

This is a script that takes each read from a big fastq file, and
uses fuzzy regular expressions (regex module in pypi)
to look for anchor sequences and use those to cut the matching
sequence groups into groups. It then repeats this for as many actions
as you specify, and then writes a new FASTQ file as specified from
the regex capture groups.

This is designed as a tool to use fixed sequences from an amplicon
sequencing experiment to find sample, strain, and UMI barcodes
from an amplicon where these are of indeterminate locations (due to
known length-variation in the strain barcodes), put the UMI into
the ID behind an underscore, and spit out the sample and strain
barcodes concatenated together for a demuxer to handle.
Here is that as an example with comments, 
kinda complicated but there's simpler
examples to follow this one:

    # First we call the script with the input fastq and an output
    # basename of `sobaseqout`
    python3 slapchop.py example_barseq.fastq sobaseqout \
        # Then we specify the first "operation". This one is called
        # "getSample" which is just used for logs and reporting. 
        # It takes "input" which always exists and matches the regex
        # following the ">".  That has several named groups called
        # "sample", "fixed1", and "rest".
        # Note the fuzzyness introduced by `{e<=1}`, this is detailed
        # later.
        -o "getSample: input > (?P<sample>[ATCG]{5})(?P<fixed1>GTCCACGAGGTC){e<=1}(?P<rest>TCT.*){e<=1}" \
        # Since those groups were captured (or not), we can access
        # them by their name. Here, "getStrainUMI" is an operation 
        # that matches this regex to the results from "rest" from
        # the previous operation, and finds the groups "tag" and
        # "strain".
        -o "getStrainUMI: rest > (?P<tag>TCT){e<=1}(?P<strain>[ATCG]{10,26})CGTACGCTGCAGGTCGAC" \
        # We can re-use the group "rest" to find "fixed2" and 
        # "aroundUMI". See how we can use wildcards even with the
        # fuzziness of `{s<=2}` specifying just two or less 
        # subsitutions.
        -o "getAroundUMI: rest > (?P<fixed2>CGTACGCTGCAGGTC)(?<aroundUMI>GAC[ATCG]G[ATCG]A[ATCG]G[ATCG]G[ATCG]G[ATCG]GAT){s<=2}" \
        #-o "getUMI: aroundUMI > (?P<wholeSet>.{16}){e<=2}" \
        # Then we can take this captured group, which may be in 
        # different locations, and we can now match each of the UMI
        # degenerate bases here as individual capture groups "umi1",
        # "umi2", "umi3", etc.
        -o "getUMI: aroundUMI > (.(?P<umi1>[ATCG]).(?<umi2>[ATCG]).(?<umi3>[ATCG]).(?<umi4>[ATCG]).(?<umi5>[ATCG]).(?<umi6>[ATCG]).){e<=2}" \
        # Now we can filter based on these expressions that are
        # going to be evaluated as expressions internally. For each
        # group, theres `*_length`, `*_start`, `*_end` calcuated,
        # so here we make sure the "sample" group is of exactly 
        # length 5. If not, it's tossed to the log report. Yes,
        # this is redundant with the regex. But what about the start
        # location? Here, I make sure it's at 0, you could say < 3
        # or whatever.
        --filter "sample_length == 5" --filter "sample_start == 0" \
        # Here we write a report. This takes resources.
        --write-report \
        # Here's how we specify resources and parralelism. Each 
        # process worker hands back and forth a filename and
        # position in the file. It takes a bite of records of size
        # specified, and then puts the new position with the 
        # filename back into the queue for the next worker.
        --bite-size 100 --processes 4 \
        # This is the output specification. What is the sequence
        # field? Here, capture groups "sample" and "strain" 
        # concatenated
        --output-seq "sample+strain" \
        # What about the ID? Here we take the original and add on
        # and underscore and then the six UMI bases. Yes, six is too
        # few, but hey, we know now.
        --output-id "input.id+'_'+umi1.seq+umi2.seq+umi3.seq+umi4.seq+umi5.seq+umi6.seq" 








Todo:

- speed ups, benchmarking
- extensive how-tos and tutorials
- Report generation, so making summary plots from these reports
- Different levels of report summary generation, ie per-read statistics or
    just distributions, or nothing
- Unit tests on choice examples from the sequencing.
