# align_annotate
alignment and annotation pipeline for short read genome sequencing, with
particular filtering and complementation strategies to enrich for homozygous
mutations from forward genetic screens in C. elegans.

Matt Rich, University of Utah 2021

## pipeline description
this is, as far as I know, a very standard alignment and annotation pipeline. it
performs basic steps: read alignment to a reference genome, filtering for high
quality alignments, and annotation of mutations (both for single nucleotide
variants as well structural variants). 

the basic pipeline proceeds as such:
1. alignment using either Bowtie2 (aln_bt2.sh) or bwa-mem (aln_bwa.sh)
2. post-processing with samtools (process_alignment.sh)
3. call variants using either GATK4 (call_variants_gatk.sh) or samtools (call_variants.sh)
4. call indels using smoove/lumpy (call_indels_smoove.sh)
5. annotate variants using snpeff, and output a human-readable tab-delimited
file containing annotations (snpeff_annotation.sh, snpeff2tsv.py).

what I tried to add to this pipeline is automation. the entire pipeline runs
through align_annotate.sh, which takes as input a working directory, your
Illumina reads, and a fasta genome. Output files are then saved in the working
directory. 

```
example of align_annotate.sh
```

The pipeline can also be automated across multiple samples using
batch_align_annotate.sh. This takes as input a tab-delimited file that maps
read identifiers (e.g., the sample name provided to or by your sequencer) to the 
names of each sample. This is currently written to work with the file provided
by the Hunstman Cancer Institute Genomics Core facility (which is who does our
genome sequencing). 

```
example of running batch_align_annotate.sh
``` 

*NOTE: currently parallelization isn't implemented here,
but certainly could be at some point. The current implementation performs the
analysis in serial (one sample at a time), but at least is a "set it and forget
it" sort of thing*

## auxiliary programs
this pipeline was written with identifying causative mutations in forward
genetic screens in *C. elegans*, and as such has a few analyses implemented that
we think aid in this.

1. comp-test.py. if your screen yields a lot of hits (like the suppressor
   screen you can find in Labella, et al. Dev Cell 2018 that isolated 48
suppressors of a single gene), it may be worthwhile to perform an *in silico*
complementation test to identify hits in the same gene. this script takes as
input the annotated TSV files (that are output from snpeff2tsv.py in
snpeff_annotation.sh) and outputs another TSV which lists the number of strains
with mutations in a given gene and the strains and their mutation.

```
example for comp_test.py
```

2. bcftools isec. an expectation for smaller screens (e.g., F2 or clonal
   screens) where you aren't screening 100,000+ genomes is that you won't
saturate the screen, and therefore most of your mutations should be independent.
As such, a useful filtering step is to only look at the unique mutations in each
strain. You can use bcftools' *isec* command to do this filtering.

```
#first move all your vcf files into a new folder
mkdir NEW_DIR
cd NEW_DIR
mv *.vcf.gz ISEC/ #or however you need to move your VCFs to ISEC/
cd ISEC/

#isec output goes into a new directory
#-n =1 tells script to output VCFs for mutations that appear in only one VCF
bcftools isec -p OUTPUT_DIR/ -n =1 *.vcf.gz

#rename files based on isec output (0000.vcf, 0001.vcf, etc...)
grep ^ISEC README.txt | bioawk -t '{system("scp ../"$1" "$3")}'
```

