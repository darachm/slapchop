# SLAPCHOP

SLAPCHOP.py parses Illumina reads using patterns to extract barcodes.

By using fuzzy regular expressions we can chop robustly barcodes from 
indeterminate positions, filter the results based on sequence or match 
properties, and reassemble a fastq record from the results.

Available as 
[a singularity containter](https://www.singularity-hub.org/collections/1361)!
So you if you have 
[Singularity](https://github.com/sylabs/singularity/releases)
installed you can just use it (without worrying about dependencies) with: 
`singularity run shub://darachm/slapchop:latest -h` (to download and show the 
argument help message for example). Then you use it like
`singularity run shub://darachm/slapchop:latest whatever.fastq arguments added`

More completely, this tool is a python script. You give it a FASTQ(Z) file and 
some operations to do, and it'll do the following:

    - read chunks of Illumina-format sequencing reads
    - apply a series of operations:
        - match fuzzy regular expression to original sequence or previous
            capture groups
        - extract capture groups and start next operation
    - apply pythonic filters (pass/fail) on sequence or quality properties 
        (like average quality or group length)
    - apply pythonic constructors to construct new FASTQ read from the capture
        groups (so ID plus the last four bases of the UMI plus length of whatever)
    - write out these reads to new files of pass and fail

For tuning/debugging/designing it has some verbosity modes to spill the gory
details of each operation in stats files, and should still have some memory
profiling functionality to debug memory leaks (fixed that one).

For a very very verbose example of a debugging run:

    ./slapchop.py \
        input.fastqz -z \
        output_basename \
        --bite-size 10 --processes 3  \
        --write-report --limit 10000 \
        -o "Sample:  input   > (?P<sample>[ATCG]{5})(?P<fixed1>GTCCACGAGGTC){e<=1}(?P<rest>TCT.*){e<=1}" \
        -o "Strain:  rest    > (?P<tag>TCT){e<=1}(?P<strain>[ATCG]{10,26})CGTACGCTGCAGGTCGAC"  \
        -o "UMITail: rest    > (?P<fixed2>CGTACGCTGCAGGTC)(?<UMItail>GAC[ATCG]G[ATCG]A[ATCG]G[ATCG]G[ATCG]G[ATCG]GAT){s<=2}"  \
        -o "UMI:     UMItail > (GAC(?P<umi1>[ATCG])G(?<umi2>[ATCG])A(?<umi3>[ATCG])G(?<umi4>[ATCG])G(?<umi5>[ATCG])G(?<umi6>[ATCG])G){e<=2}"  \
        --output-seq "strain" \
        --output-id "input.id+'_umi='+umi1.seq+umi2.seq+umi3.seq+ \
            umi4.seq+umi5.seq+umi6.seq+'_sample='+sample.seq" \
        --filter "sample_length == 5 and rest_start >= 16 and ( min(strain.letter_annotations['phred_quality']) >= 30 )"\
        --verbose --verbose --verbose

That invocation:

    - Takes records from the `input.fastq`
    - Starts three processes that each take bites of 10 records
    - Applies the four operations to cut up the read
    - Writes the full detailed report including json reports for 
        each read, so we limit it to the first 10,000 bytes
        of the file (about 50 records). This is for debugging.
    - Filters the records on having a `sample` barcode of 5 bases 
        and having the `rest` sequence match starting at least past
        index 16 (so base 15 in english).
    - Re-assembles the records that pass this filter, making the ID
        of the fastq record having the original ID plus a UMI 
        sequence and the sample barcode, then the sequence is just
        the match to the strain barcode context. This is suitable for
        feeding into `bwa` for example.
    - We've got three levels of verbosity, so a per-record verbosity
        for debugging purposes.

Note that the output formatting is a bit funny. This is directly evaluated
(because security is what?) on BioPython SequenceRecords, so you need to specify
just the name of the capture group(s) for the outputs so it can access the
`.seq` and qualities. For the ID, etc, you can access `.seq` or `.id`.

Then if we like our thresholds we'd re-run, and drop the `--limit`
and `--write-report` flags. This will turn records like this:

    @NB501157:100:H5J5LBGX2:1:11101:10000:6068 1:N:0:
    CTACTGTCCACGAGGTCTCTGCAGATAATACACTGTCACCCGTACGCTGCAGGTCGACCGTAGGAGGGAGATGTG
    +
    AAAAAEEEE/AEE<EEEEEEEEAEEAEEAEEEEE/EEE/EEEEEEEEE/EEEEEEEEEEEEE/EEEEEEEEEEEE

into records like this:

    @NB501157:100:H5J5LBGX2:1:11101:10000:6068_umi=CTGAGA_sample=CTACT
    GCAGATAATACACTGTCACC
    +
    EEAEEAEEAEEEEE/EEE/E

The sample barcode is the first five, the strain barcode starts after
the `TCT`, and the UMI is interspersed downstream. This is modified
yeast BarSeq, btw.

---

This script depends strongly upon (uses) the work of 
[regex](https://pypi.org/project/regex/)
and
[Biopython](https://pypi.org/project/biopython/). Thanks! Check them out...

