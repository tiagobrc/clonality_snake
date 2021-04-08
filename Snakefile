import pandas as pd

samples_df = pd.read_table('data/runs/primeira/samples.tsv').set_index("samples", drop=False)
SAMPLES = list(samples_df['samples'])
SOURCE = list(samples_df['source'])
DESTINY = list(samples_df['destiny'])
CELL = list(samples_df['cell'])
HEAD = list(samples_df['read'])
EMAIL = list(samples_df['email'])
SERVER = list(samples_df['server'])

PBC_df = pd.read_table('data/barcodes/T/PBC.txt', header=None)
PBC = list(PBC_df[0])

rule all:
    input:
        expand("{destiny}{sample}.zip", destiny=DESTINY, sample=SAMPLES)

#run fastqc for quality control
#trim 2nt in the begining of each read
rule trim_reads:
    input:
        expand("{source}{sample}_S1_R1_001.fastq.gz", source=SOURCE, sample=SAMPLES),
        expand("{source}{sample}_S1_R2_001.fastq.gz", source=SOURCE, sample=SAMPLES)
    output:
        expand("{destiny}{sample}_S1_R1_001.fastq", sample=SAMPLES, destiny=DESTINY),
        expand("{destiny}{sample}_S1_R2_001.fastq", sample=SAMPLES, destiny=DESTINY)
    params:
        QCinput=expand("{source}", source=SOURCE),
        QCfolder=expand("{destiny}", destiny=DESTINY)
    conda:
        "envs/fastx.yaml"
    shell:"""
        zcat {input[0]} | fastx_trimmer -Q33 -f 3 -o {output[0]} && zcat {input[1]} | fastx_trimmer -Q33 -f 3 -o {output[1]}
        fastqc {input[0]}
        fastqc {input[1]}
        mkdir {params.QCfolder}QC
        mv {params.QCinput}*zip {params.QCinput}*html {params.QCfolder}QC/
    """

#run pandaseq
rule pandaseq:
    input:
        expand("{source}{sample}_S1_R1_001.fastq", source=DESTINY, sample=SAMPLES),
        expand("{source}{sample}_S1_R2_001.fastq", source=DESTINY, sample=SAMPLES)
    output:
        expand("{destiny}{sample}.fasta", sample=SAMPLES, destiny=DESTINY)
    conda:
        "envs/pandaseq.yaml"
    shell:
        "pandaseq -f {input[0]} -r {input[1]} -w {output}"

#Split fasta files into plates by each barcodes
rule genfasta:
    input:
        expand("{source}{sample}.fasta", source=DESTINY, sample=SAMPLES),
        expand("data/barcodes/{cell}PBC.txt", cell=CELL, sample=SAMPLES),
        expand("data/barcodes/{cell}BC.txt", cell=CELL, sample=SAMPLES)
    conda:
        "envs/fastx.yaml"
    params:
        fo1=expand("{destiny}split1/", destiny=DESTINY),
        fo2=expand("{destiny}split1/split2/", destiny=DESTINY),
        fo3=expand("{destiny}split1/split2/collapsed/", destiny=DESTINY),
        head=expand("{head}", head=HEAD),
        fo4=expand("{destiny}split1/split2/collapsed/head/", destiny=DESTINY),
        fo5=expand("{destiny}split1/split2/collapsed/head/final/", destiny=DESTINY),
        sample=expand("{sample}", sample=SAMPLES),
        removal=expand("{destiny}", destiny=DESTINY)
    output:
        expand("{destiny}split1/{sample}", destiny=DESTINY, sample=SAMPLES)
    shell:"""
        touch {output}
        cat {input[0]} | fastx_barcode_splitter.pl --bcfile {input[1]} --bol --prefix {output} --suffix .fasta --exact
        rm  {params.fo1}*unmatched*
        find {params.fo1} -name "*fasta" -size -4k -delete
        mkdir -p {params.fo1}split2/
        for f in {params.fo1}*fasta ; do x="${{f##*/}}" ; cat $f | fastx_barcode_splitter.pl --bcfile {input[2]} --eol --prefix {params.fo1}split2/${{x%.*}} --suffix .fasta ; done > {params.removal}BC.stats.txt
        grep -v 'total' {params.removal}BC.stats.txt > {params.removal}tmp && mv {params.removal}tmp {params.removal}BC.stats.txt
        grep -v 'Barcode' {params.removal}BC.stats.txt > {params.removal}tmp && mv {params.removal}tmp {params.removal}BC.stats.txt
        rm {params.fo2}*unmatched*
        find {params.fo2} -name "*fasta" -size -4k -delete
        mkdir -p {params.fo2}collapsed/
        for f in {params.fo2}*fasta ; do x="${{f##*/}}" ; fastx_collapser -i $f -o {params.fo2}collapsed/${{x%.*}}.fasta ; done

        for f in {params.fo3}*fasta ; do x="${{f##*/}}" ; grep ">" $f | sed s/^.*-//g | sort -nr > {params.fo3}${{x%.*}}.collapsed.stats.txt ; done
        ls {params.fo3}*txt | xargs | sed 's/ /\t/g' | sed 's/.collapsed.stats.txt//g' > {params.fo3}col.names
        paste {params.fo3}*txt > {params.fo3}collapsed.stats
        cat {params.fo3}col.names {params.fo3}collapsed.stats > {params.removal}collapsed.stats.txt
        Rscript ./data/scripts/plate_stats.R {params.removal} {params.removal}{params.sample}

        mkdir -p {params.fo3}head/
        for f in {params.fo3}*fasta ; do x="${{f##*/}}" ; head -n{params.head} $f > {params.fo3}head/${{x%.*}}.fasta ; done
        mkdir -p {params.fo4}final/
        awk '/>/{{sub(">","&"FILENAME"_")}}1' {params.fo4}*fasta > {params.fo5}{params.sample}.fasta
        sed -i 's/>.*\//>/g' {params.fo5}{params.sample}.fasta
        sed -i 's/\.fasta//g' {params.fo5}{params.sample}.fasta
    """

rule final:
    input:
        expand("{destiny}split1/{sample}", destiny=DESTINY, sample=SAMPLES),
    output:
        expand("{destiny}{sample}.zip", destiny=DESTINY, sample=SAMPLES)
    params:
        removal=expand("{destiny}", destiny=DESTINY),
        finalFasta=expand("{destiny}split1/split2/collapsed/head/final/{sample}.fasta", destiny=DESTINY, sample=SAMPLES),
        finalRes=expand("{destiny}{sample}.fasta", destiny=DESTINY, sample=SAMPLES),
        QC=expand("{destiny}/QC", destiny=DESTINY),
        sample=expand("{sample}.zip", sample=SAMPLES),
        server=expand("{server}", server=SERVER),
        email=expand("{email}", email=EMAIL),
        subject=expand("{sample}", sample=SAMPLES)
    shell:"""
        mv {params.finalFasta} {params.finalRes}
        zip -r {output} {params.finalRes} {params.QC} {params.removal}*pdf {params.removal}*tsv
        if [ {params.server} == "mucidalab3" ]; then smbclient '//mucidalab3/LMI_shared' -A ~/.smbclient1.conf -c 'cd ./Pipelines/clonality/results/ ; put {output} {params.sample}'; fi
        if [ {params.server} == "victoranas" ]; then smbclient '//victoranas/victoracfs'  -A ~/.smbclient2.conf -c 'cd ./Pipelines/clonality/results/; put {output} {params.sample}'; fi
        mail -s "Your Clonality Run {params.subject} is completed" {params.email} <<< "Your clonality run is completed sucessfuly and it was uploaded to the folder: Pipelines/clonality/results/ of your lab server."
    """
