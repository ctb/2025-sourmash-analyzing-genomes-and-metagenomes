tld = '/group/ctbrowngrp2/scratch/annie/2023-swine-sra/results/atlas/'
bins = tld + 'atlas_{acc}/genomes/genomes/'
reads = tld + 'atlas_{{acc}}/{{acc}}/sequence_quality_control/{{acc}}_QC_R{r}.fastq.gz'
assembly = tld + 'atlas_{acc}/{acc}/{acc}_contigs.fasta'

ACCS = ['SRR11125434']
SOURCES = ['reads', 'bins', 'assembly']

rule all:
    input:
        expand("sketches/{acc}_{src}.sig.zip", acc=ACCS, src=SOURCES),
        expand("outputs/{acc}_{src}.x.gtdb-rs220.gather.csv",
               acc=ACCS, src=SOURCES),

rule sketch_reads:
    input:
        reads=expand(reads, r=[1, 2])
    output:
        "sketches/{acc}_reads.sig.zip",
    shell: """
        sourmash scripts singlesketch -p dna,abund,k=21,k=31,k=51 {input} \
           -o {output} --name "{wildcards.acc}_reads"
    """

rule sketch_bins:
    input:
        bindir=bins
    output:
        "sketches/{acc}_bins.sig.zip",
    shell: """
        sourmash scripts singlesketch -p dna,k=21,k=31,k=51 {input}/* \
            -o {output} --name "{wildcards.acc}_bins"
    
    """
        
rule sketch_assembly:
    input:
        assembly=assembly,
    output:
        "sketches/{acc}_assembly.sig.zip",
    shell: """
        sourmash scripts singlesketch -p dna,k=21,k=31,k=51 {input} -o {output} \
            --name "{wildcards.acc}_assembly"
    """
        
rule gather_gtdb:
    input:
        q="sketches/{sketch}.sig.zip",
        db="/group/ctbrowngrp5/sourmash-db/gtdb-rs220/gtdb-rs220-k31.zip",
    output:
        "outputs/{sketch}.x.gtdb-rs220.gather.csv",
    threads: 32
    shell: """
        sourmash scripts fastgather {input.q} {input.db} \
            -o {output} -k 31 -s 1000 -t 0 -c {threads}
    """
