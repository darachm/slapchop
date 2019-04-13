
.PHONY: container examples clean


example: tmp/barseq_shortrun_pass.fastq tmp/barseq_longrun_pass.fastq

container_test: tmp/barseq_shortrun_pass_container.fastq tmp/barseq_longrun_pass_container.fastq

tmp/barseq_shortrun_pass_container.fastq: example-data/barseq.fastq slapchop.singularity
	mkdir -p tmp
	singularity run slapchop.singularity \
		example-data/barseq.fastq tmp/barseq_shortrun \
		--bite-size 10 --processes 3 \
		--write-report --limit 1000 \
		-o "Sample:  input   > (?P<sample>[ATCG]{5})(?P<fixed1>GTCCACGAGGTC){e<=1}(?P<rest>TCT.*){e<=1}" \
		-o "Strain:  rest	> (?P<tag>TCT){e<=1}(?P<strain>[ATCG]{10,26})CGTACGCTGCAGGTCGAC" \
		-o "UMITail: rest    > (?P<fixed2>CGTACGCTGCAGGTC)(?<UMItail>GAC[ATCG]G[ATCG]A[ATCG]G[ATCG]G[ATCG]G[ATCG]GAT){s<=2}" \
		-o "UMI:     UMItail > (GAC(?P<umi1>[ATCG])G(?<umi2>[ATCG])A(?<umi3>[ATCG])G(?<umi4>[ATCG])G(?<umi5>[ATCG])G(?<umi6>[ATCG])G){e<=2}" \
		--output-seq "sample+spacer+strain" \
		--output-id "input.id+'_umi='+umi1.seq+umi2.seq+umi3.seq+ \
			umi4.seq+umi5.seq+umi6.seq+'_sample='+sample.seq" \
		--filter "sample_length == 5 and rest_start >= 16" \
		--verbose --verbose --verbose --verbose

tmp/barseq_longrun_pass_container.fastq: example-data/barseq.fastq slapchop.singularity
	mkdir -p tmp
	singularity run slapchop.singularity \
		example-data/barseq.fastq tmp/barseq_longrun \
		--processes 3 \
		-o "Sample:  input   > (?P<sample>[ATCG]{5})(?P<fixed1>GTCCACGAGGTC){e<=1}(?P<rest>TCT.*){e<=1}" \
		-o "Strain:  rest	 > (?P<tag>TCT){e<=1}(?P<strain>[ATCG]{10,26})CGTACGCTGCAGGTCGAC" \
		-o "UMITail: rest    > (?P<fixed2>CGTACGCTGCAGGTC)(?<UMItail>GAC[ATCG]G[ATCG]A[ATCG]G[ATCG]G[ATCG]G[ATCG]GAT){s<=2}" \
		-o "UMI:     UMItail > (GAC(?P<umi1>[ATCG])G(?<umi2>[ATCG])A(?<umi3>[ATCG])G(?<umi4>[ATCG])G(?<umi5>[ATCG])G(?<umi6>[ATCG])G){e<=2}" \
		--output-seq "sample+spacer+strain" \
		--output-id "input.id+'_umi='+umi1.seq+umi2.seq+umi3.seq+ \
			umi4.seq+umi5.seq+umi6.seq+'_sample='+sample.seq" \
		--filter "sample_length == 5 and rest_start >= 16" \
		--verbose --verbose

tmp/barseq_shortrun_pass.fastq: example-data/barseq.fastq slapchop.singularity
	mkdir -p tmp
	python3 slapchop.py \
		example-data/barseq.fastq tmp/barseq_shortrun \
		--bite-size 10 --processes 3 \
		--write-report --limit 1000 \
		-o "Sample:  input   > (?P<sample>[ATCG]{5})(?P<fixed1>GTCCACGAGGTC){e<=1}(?P<rest>TCT.*){e<=1}" \
		-o "Strain:  rest	> (?P<tag>TCT){e<=1}(?P<strain>[ATCG]{10,26})CGTACGCTGCAGGTCGAC" \
		-o "UMITail: rest    > (?P<fixed2>CGTACGCTGCAGGTC)(?<UMItail>GAC[ATCG]G[ATCG]A[ATCG]G[ATCG]G[ATCG]G[ATCG]GAT){s<=2}" \
		-o "UMI:     UMItail > (GAC(?P<umi1>[ATCG])G(?<umi2>[ATCG])A(?<umi3>[ATCG])G(?<umi4>[ATCG])G(?<umi5>[ATCG])G(?<umi6>[ATCG])G){e<=2}" \
		--output-seq "sample+spacer+strain" \
		--output-id "input.id+'_umi='+umi1.seq+umi2.seq+umi3.seq+ \
			umi4.seq+umi5.seq+umi6.seq+'_sample='+sample.seq" \
		--filter "sample_length == 5 and rest_start >= 16" \
		--verbose --verbose --verbose --verbose

tmp/barseq_longrun_pass.fastq: example-data/barseq.fastq slapchop.singularity
	mkdir -p tmp
	python3 slapchop.py \
		example-data/barseq.fastq tmp/barseq_longrun \
		--bite-size 10 --processes 3 \
		-o "Sample:  input   > (?P<sample>[ATCG]{5})(?P<fixed1>GTCCACGAGGTC){e<=1}(?P<rest>TCT.*){e<=1}" \
		-o "Strain:  rest	 > (?P<tag>TCT){e<=1}(?P<strain>[ATCG]{10,26})CGTACGCTGCAGGTCGAC" \
		-o "UMITail: rest    > (?P<fixed2>CGTACGCTGCAGGTC)(?<UMItail>GAC[ATCG]G[ATCG]A[ATCG]G[ATCG]G[ATCG]G[ATCG]GAT){s<=2}" \
		-o "UMI:     UMItail > (GAC(?P<umi1>[ATCG])G(?<umi2>[ATCG])A(?<umi3>[ATCG])G(?<umi4>[ATCG])G(?<umi5>[ATCG])G(?<umi6>[ATCG])G){e<=2}" \
		--output-seq "sample+spacer+strain" \
		--output-id "input.id+'_umi='+umi1.seq+umi2.seq+umi3.seq+ \
			umi4.seq+umi5.seq+umi6.seq+'_sample='+sample.seq" \
		--filter "sample_length == 5 and rest_start >= 16" \
		--verbose --verbose \
		-m 

clean:
	rm -r tmp || echo ""
	rm -r .nextflow* || echo ""
	rm -r nextflow* || echo ""
	rm -r work || echo ""
	
slapchop.singularity: Singularity
	sudo rm -r $@ || echo "already gone"
	sudo singularity build $@ $<

/tmp/slapchop-test-simg: Singularity
	sudo rm -r $@ || echo "already gone"
	sudo singularity build --sandbox $@ $<

