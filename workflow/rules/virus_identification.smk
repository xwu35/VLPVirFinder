rule genomad:
    """
    Viral contig identification using geNomad
    """
    input: 
        contigs=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta")
    output:
        genomad=directory(os.path.join(dir["output"]["viral_identification"], "genomad", "{sample}")),
        summary=os.path.join(dir["output"]["viral_identification"], "genomad", "{sample}", "{sample}_contigs_1kb_summary", "{sample}_contigs_1kb_virus_summary.tsv"),
        genomad_viruses=os.path.join(dir["output"]["viral_identification"], "genomad", "{sample}", "{sample}_contigs_1kb_summary", "{sample}_contigs_1kb_virus.fna")
    params:
        database=config["genomad"]["database_path"],
        sensitivity=config["genomad"]["sensitivity"], 
        other=config["genomad"]["other"]
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["big_mem"]
    shell:
        """
        {GENOMAD}

        if [ -s "{input.contigs}" ]; then
            genomad end-to-end \
                --cleanup {input.contigs} \
                {output.genomad} \
                {params.database} \
                -t {threads} \
                --sensitivity {params.sensitivity} \
                {params.other}
        else
            echo "File is empty"
            touch {output.summary}
            touch {output.genomad_viruses} # this one is necessary for extracting genomad identified provirus
        fi   
        """

rule virsorter2:
    """
    Viral contig identification using VirSorter2
    """
    input: 
        contigs=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta")
    output:
        virsorter2_outdir=directory(os.path.join(dir["output"]["viral_identification"], "virsorter2-pass1", "{sample}")),
        score=os.path.join(dir["output"]["viral_identification"], "virsorter2-pass1", "{sample}", "final-viral-score.tsv")
    params:
        groups=config["virsorter2"]["groups"],
        length=config["virsorter2"]["length"],
        score=config["virsorter2"]["score"],
        other=config["virsorter2"]["other"],
        steps=config["virsorter2"]["steps"]
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["big_mem"]
    shell:
        """
        {VIRSORTER2}

        if [ -s "{input.contigs}" ]; then
            virsorter run \
                {params.other} \
                -i {input.contigs} \
                -w {output.virsorter2_outdir} \
                --include-groups {params.groups} \
                --min-length {params.length} \
                --min-score {params.score} \
                -j {threads} {params.steps}
        else
            echo "Contig file is empty"
            touch {output.score}
        fi   
        """

rule cenote_taker3:
    """
    Viral contig identification using Cenote-Taker3.
    """
    input: 
        contigs=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta")
    output:
        summary=os.path.join(dir["output"]["viral_identification"], "cenote_taker3", "{sample}", "{sample}_virus_summary.tsv"),
        done=os.path.join(dir["output"]["viral_identification"], "cenote_taker3", "{sample}", ".done")
    params:
        # because snakemake generates the directory with sample name, cenote-taker3 will rename it (e.g. change VLP_Day0 to VLP_Day0_old_OJ5RL)
        # so here use a tmp directory and then move everything to destination directory
        tmp_dir=directory(os.path.join(dir["output"]["viral_identification"], "cenote_taker3_tmp")),
        dst_dir=directory(os.path.join(dir["output"]["viral_identification"], "cenote_taker3")),
        database=config["cenote_taker3"]["database"],
        # because the {} in -db {virion,rdrp,dnarep} makes snakemake thinking they are wildcards, 
        # therefore need to deactivate automatic wildcard expansion in params strings by specifying `lambda wildcards: `
        settings=lambda wildcards: config["cenote_taker3"]["settings"] 
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    shell:
        """
        {CT3}

        if [ -s "{input.contigs}" ]; then
            cenotetaker3 -c {input.contigs} \
                -r {wildcards.sample} \
                -t {threads} \
                -wd {params.tmp_dir} \
                --cenote-dbs {params.database} \
                {params.settings} && # && lets you do something based on whether the previous command completed successfully
            mv {params.tmp_dir}/{wildcards.sample}/* {params.dst_dir}/{wildcards.sample} &&
            rm -r {params.tmp_dir} &&
            touch {output.done}
        else
            echo "Contig file is empty"
            touch {output.summary}
            touch {output.done}
        fi   
        """
