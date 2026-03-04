rule select_viral_contigs:
    """
    Contigs satisfying one of the following criteria were considered viral:
    1) geNomad: all identified contigs
    2) Cenote-taker3: all identified contigs
    3) VirSorter2: all identified contigs
    """
    input:
        genomad=os.path.join(dir["output"]["viral_identification"], "genomad", "{sample}", "{sample}_contigs_1kb_summary", "{sample}_contigs_1kb_virus_summary.tsv"),
        cenote=os.path.join(dir["output"]["viral_identification"], "cenote_taker3", "{sample}", "{sample}_virus_summary.tsv"),
        virsorter2=os.path.join(dir["output"]["viral_identification"], "virsorter2-pass1", "{sample}", "final-viral-score.tsv"),
        done=os.path.join(dir["output"]["viral_identification"], "cenote_taker3", "{sample}", ".done")
    output:
        viral=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_lists", "{sample}_viral_contig_names.tsv"),
        provirus=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_lists", "{sample}_genomad_provirus_contig_names.tsv")
    params:
        script=os.path.join(dir["scripts"], "baldridge_select_viral_contigs.R")
    conda:
        os.path.join(dir["env"], "R.yml")
    shell:
        """
        Rscript {params.script} \
            -g {input.genomad} \
            -c {input.cenote} \
            -v {input.virsorter2} \
            -o {output.viral} \
            -p {output.provirus}
        """

rule extract_identified_viral_contig:
    """
    extract_identified_viral_contig sequences. 
    NOTE: proviruses identified by geNomad were excluded here
    """
    input:
        viral_id=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_lists", "{sample}_viral_contig_names.tsv"),
        contigs=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta")
    output:
        viral_seq=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_sequences", "tmp", "{sample}_identified_viral_contigs_without_genomad_provirus.fa")
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    shell:
        """
        # if viral_id is not empty, extract their sequences
        if [ -s "{input.viral_id}" ]; then
            seqkit grep -i -f {input.viral_id} {input.contigs} > {output.viral_seq}
        else
            echo "viral_id file is empty"
            touch {output.viral_seq}
        fi
        """

rule extract_genomad_provirus_sequences:
    """
    get contig names of provirus identified by geNomad. TODO: double check in termina as PvVT has none
    """
    input:
        provirus_id=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_lists", "{sample}_genomad_provirus_contig_names.tsv"),
        genomad_viruses=os.path.join(dir["output"]["viral_identification"], "genomad", "{sample}", "{sample}_contigs_1kb_summary", "{sample}_contigs_1kb_virus.fna")
    output:
        genomad_provirus=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_sequences", "tmp", "{sample}_genomad_provirus.fa"),
        genomad_provirus_1kb=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_sequences", "tmp", "{sample}_genomad_provirus_1kb.fa")
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    shell:
        """
        # if provirus_id is not empty, extract their sequences
        if [ -s "{input.provirus_id}" ]; then
            seqkit grep -i -f {input.provirus_id} {input.genomad_viruses} > {output.genomad_provirus}
            # only keep trimmed provirus with length >= 1kb 
            seqkit seq -m 1000 {output.genomad_provirus} > {output.genomad_provirus_1kb}
        else
            echo "provirus_id file is empty"
            touch {output.genomad_provirus}
            touch {output.genomad_provirus_1kb}
        fi
        """

rule combine_nonprovirus_and_genomad_provirus:
    """
    combine the non-provirus viral contigs with the genomad proviruss
    """
    input:
        genomad_provirus_1kb=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_sequences", "tmp", "{sample}_genomad_provirus_1kb.fa"),
        viral_seq=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_sequences", "tmp", "{sample}_identified_viral_contigs_without_genomad_provirus.fa")
    output:
        with_provirus=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_sequences", "tmp", "{sample}_identified_viral_contigs_with_genomad_provirus.fa")
    shell:
        """
        cat {input.genomad_provirus_1kb} {input.viral_seq} > {output.with_provirus}
        """

rule checkv_identified:
    """
    run checkV on all identified viral contigs to identify additional proviruses.
    Contigs identified by VirSorter2 and Cenote-taker2 were not checked for provirus, so run checkV on the identified viral contigs and extract the trimmed provirus contigs.
    """
    input:
        with_provirus=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_sequences", "tmp", "{sample}_identified_viral_contigs_with_genomad_provirus.fa")
    output:
        checkv_dir=directory(os.path.join(dir["output"]["viral_identification"],"checkv_identified", "{sample}")),
        checkv_proviruses=os.path.join(dir["output"]["viral_identification"],"checkv_identified", "{sample}", "proviruses.fna"), # this would NOT be a problem since CheckV produces empty provirus files
        checkv_viruses=os.path.join(dir["output"]["viral_identification"],"checkv_identified", "{sample}", "viruses.fna")
    params:
        database_path=config["checkv"]["database_path"]
    threads:
        config["resources"]["small_cpu"]
    resources:
        mem_mb=config["resources"]["med_mem"]
    shell:
        """
        {CHECKV}

        # if with_provirus is not empty, run checkV
        if [ -s "{input.with_provirus}" ]; then
            checkv end_to_end {input.with_provirus} {output.checkv_dir} -t {threads} -d {params.database_path}
        else
            echo "with_provirus file is empty"
            touch {output.checkv_proviruses}
            touch {output.checkv_viruses}
        fi
        """

rule combine_checkv_provirus_1kb_and_virus:
    """
    CheckV might have identified more proviruses, only keep proviruses >= 1kb.
    Combine the kept proviruses and the non-proviruses sequences
    """
    input:
        checkv_proviruses=os.path.join(dir["output"]["viral_identification"],"checkv_identified", "{sample}", "proviruses.fna"), # this would NOT be a problem since CheckV produces empty provirus files
        checkv_viruses=os.path.join(dir["output"]["viral_identification"],"checkv_identified", "{sample}", "viruses.fna")
    output:
        checkv_proviruses_1kb=os.path.join(dir["output"]["viral_identification"],"checkv_identified", "{sample}", "proviruses_1kb.fna"),
        final_viruses=os.path.join(dir["output"]["viral_identification"],"identified_viral_contig_sequences", "{sample}_identified_viral_contigs_final.fa")
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    shell:
        """
        # only keep trimmed provirus with length >= 1kb 
        # deal with the possibility that no provirus identified by checkV
        if [ -s "{input.checkv_proviruses}" ]; then
            seqkit seq -m 1000 {input.checkv_proviruses} > {output.checkv_proviruses_1kb}
        else
            touch {output.checkv_proviruses_1kb}
        fi

        # combine proviruses >= 1kb with viruses
        cat {output.checkv_proviruses_1kb} {input.checkv_viruses} > {output.final_viruses}
        """

rule get_virus_sequence_length:
    """
    Get the length of final viruses
    """
    input:
        final_viruses=os.path.join(dir["output"]["viral_identification"],"identified_viral_contig_sequences", "{sample}_identified_viral_contigs_final.fa")
    output:
        os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_sequences", "{sample}_identified_viral_contigs_final_length.txt")
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    shell:
        """
        # DON'T NEED TO DEAL WITH EMPTY FASTA, seqkit WILL JUST OUTPUT HEADER
        seqkit fx2tab --length --name --header-line {input.final_viruses} > {output}
        """