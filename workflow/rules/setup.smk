#------------ SET UP THE DIRECTORIES
dir = dict()
dir["output"] = dict()

# WORKFLOW DIRs
dir["env"]     = os.path.join(workflow.basedir, "envs")
dir["scripts"] = os.path.join(workflow.basedir, "scripts")

# OUTPUT DIRs
dir["output"]["base"] = RESULTS_DIR
dir["output"]["reads_processing"] = os.path.join(dir["output"]["base"], "reads_processing")
dir["output"]["assembly"] = os.path.join(dir["output"]["base"], "assembly")
dir["output"]["intermediate"] = os.path.join(dir["output"]["reads_processing"], "intermediate")
dir["output"]["viral_identification"] = os.path.join(dir["output"]["base"], "viral_identification")

#------------ SET UP THE OUTPUT
trimmed_reads=[]
contig_status=[]
viral_identification=[]
for sample in SAMPLE:
    trimmed_reads.append(os.path.join(dir["output"]["reads_processing"], "filtered_reads", sample + "_R1.humanfilt.fastq.gz"))
    trimmed_reads.append(os.path.join(dir["output"]["reads_processing"], "filtered_reads", sample + "_R2.humanfilt.fastq.gz"))
    contig_status.append(os.path.join(dir["output"]["assembly"], "contig_evaluation", sample, "report.txt"))
    viral_identification.append(os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_sequences", sample + "_identified_viral_contigs_final_length.txt"))

# raw reads: read counts and fastqc
fastqc_input = [
    os.path.join(dir["scripts"], "anicalc.py"),
    os.path.join(dir["scripts"], "aniclust.py"),
    os.path.join(dir["output"]["reads_processing"], "reads_statistics", "raw_reads", "R1_stats.tsv"),
    os.path.join(dir["output"]["reads_processing"], "reads_statistics", "raw_reads", "R2_stats.tsv"),
    os.path.join(dir["output"]["reads_processing"], "fastqc", "multiqc_raw_reads", "multiqc_report.html")
]

# clean reads: removed adaptors, phix and human contmaination
clean_reads_input = [os.path.join(dir["output"]["reads_processing"], "fastqc", "multiqc_after_trimmomatic", "multiqc_report.html")] + trimmed_reads + [
    os.path.join(dir["output"]["reads_processing"], "reads_statistics", "number_of_reads_removed_at_each_step.txt"),
    os.path.join(dir["output"]["reads_processing"], "reads_statistics", "reads_composition_barplot.svg")
]

# assembly 
assembly_input = contig_status + [os.path.join(dir["output"]["assembly"], "all_contigs_1kb.fasta")] 

# viral identification
identification_input = viral_identification

#---------- DOWNLOAD SCRIPTS
localrules:
    download_checkV_scripts

rule download_checkV_scripts:
    """
    download scripts for calculating ANI and clustering
    """
    output:
        anicalc=os.path.join(dir["scripts"], "anicalc.py"),
        aniclust=os.path.join(dir["scripts"], "aniclust.py")
    params:
        dir=dir["scripts"]
    shell:
        """
        # download the two scripts from checkV which cannot be installed via conda
        # only download them if the files don't exsit
        # script for calculating ANI and AF
        if [[ ! -e {output.anicalc} ]]; then 
            wget -P {params.dir} https://bitbucket.org/berkeleylab/checkv/raw/3f185b5841e8c109848cd0b001df7117fe795c50/scripts/anicalc.py && chmod +x {output.anicalc}
        else
            # in case the files were downloaded, but not executable 
            chmod +x {output.anicalc}
        fi

        # script for clustering contigs into vOTUs
        if [[ ! -e {output.aniclust} ]]; then 
            wget -P {params.dir} https://bitbucket.org/berkeleylab/checkv/raw/3f185b5841e8c109848cd0b001df7117fe795c50/scripts/aniclust.py && chmod +x {output.aniclust}
        else
            # in case the files were downloaded, but not executable
            chmod +x {output.aniclust}
        fi
        """
