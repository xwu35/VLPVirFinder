# VLPVirFinder

VLPVirFinder is a Snakemake workflow designed to identify viral contigs from viromes generated using MDA amplification at the Baldridge lab. 

The quality of sequencing reads is assessed using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.12.1 and [MultiQC](https://github.com/MultiQC/MultiQC) v1.14. Adapter contamination and low-quality reads are removed using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) v0.39. Reads mapped to the Phix and human genomes are removed using the selected mapping tool ([Bowtie2](https://github.com/BenLangmead/bowtie2) v2.4.2 or [Minimap2](https://github.com/lh3/minimap2) v2.28) and [SAMtools](https://github.com/samtools/samtools) v1.21. The identification is carried out using [GeNomad](https://github.com/apcamargo/genomad) v1.11.1 (--sensitivity 7.5 --enable-score-calibration), [VirSorter2](https://github.com/jiarong/VirSorter2) v2.2.4 (--keep-original-seq --include-groups dsDNAphage,ssDNA --min-length 1000 --min-score 0.5) and [Cenote-Taker3](https://github.com/mtisza1/Cenote-Taker3) v3.4.3 (-p F -db {virion,rdrp,dnarep}). Contigs are considered viral if they were identified by at least one of the methods. [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/) is used to evaluate the quality of the identified viral contigs. The host regions from the edges of contigs that were identified as proviruses by geNomad and CheckV are removed, and only proviruses with a length >= 1kb after removing the host regions are kept.

**NOTE**: geNomad, VirSorter2, Cenote-Taker3 and CheckV must be installed before running the workflow. The paths to these tools and their databases need to be updated to match your local setup by editing `VLPVirFinder/config/config.yml`. Automatic installation will be added in the next version.

## Set up environment

### Clone the repository 

```bash
git clone https://github.com/xwu35/VLPVirFinder.git
```

### Install Snakemake

VLPVirFinder is built for Snakemake version 7. Version 8 and above introduce breaking changes and deprecations and have not been tested. It may not function correctly with newer versions. Please install Snakemake version 7 using the script below.

```bash
cd VLPVirFinder

# option 1 using conda
conda env create -n snakemake -f snakemake_env.yml

# option 2 using mamba if it's installed
mamba env create -n snakemake -f snakemake_env.yml
```

### Download snakemake profile

The profile is required to run the workflow on HPC. Skip this step if you already have a SLURM profile in `~/.config/snakemake`.

```bash
# download the profile
git clone https://github.com/xwu35/slurm

# move the profile to the right directory
mv slurm ~/.config/snakemake 
```

## Sample information table

The sample information table should look like this:

| sample   | R1                                    | R2                                    |
|----------|---------------------------------------|---------------------------------------|
| Fp22_10A | Baldridge_10A_Fp22_CD_R1_001.fastq.gz | Baldridge_10A_Fp22_CD_R2_001.fastq.gz | 
| Fp22_3C  | Baldridge_3C_Fp22_HHC_R1_001.fastq.gz | Baldridge_3C_Fp22_HHC_R2_001.fastq.gz | 
| Fp22_7B  | Baldridge_7B_Fp22_CD_R1_001.fastq.gz  | Baldridge_7B_Fp22_CD_R2_001.fastq.gz  | 
| Fp22_8G  | Baldridge_8G_Fp22_CD_R1_001.fastq.gz  | Baldridge_8G_Fp22_CD_R2_001.fastq.gz  | 

## Usage

VLPVirFinder supports two mapping software options (bowtie2 and minimap2) and two assembler options (megahit and metaspades), with minimap2 and metaspades used by default. Detailed usage information can be viewed using the -h or --help flags `python VLPVirFinder.py -h`.

Do not run the analysis on the login node. Submit it as a sbatch job. See `run_vlpvirfinder.sh` for an example, or check the HTCF usage guide here (https://github.com/xwu35/baldridge_lab/blob/main/HTCF.md). 

A dry-run can be performed to check which rules will be executed and which files will be produced by specifying `--dryrun`.

```bash
conda activate snakemake

python VLPVirFinder.py \
    --reads_dir /path/to/dir/containing/sequences \
    --sample_info /path/to/sample_info_table \
    --output_dir output_dir 
```

### Specific steps

Specific steps can be run using the `--step` flag. 

- **fastqc**: QC on raw reads
- **preprocess**: run fastqc and trim the reads
- **assemble**: run fastqc, preprocess, and assembly
- **identification**: run all the steps (fastqc, preprocess, assemble and viral contigs identification)

VLPVirFinder runs all steps by default.

## Output description

|                     Filename                                                                            |                     Description                    |
|---------------------------------------------------------------------------------------------------------|----------------------------------------------------|
| `reads_processing/fastqc`                                                                               | Quality control results                            |
| `reads_processing/filtered_reads`                                                                       | Trimmed reads used for assembly                    |
| `reads_statistics/number_of_reads_removed_at_each_step.txt`                                             | Read counts                                        |
| `assembly/all_contigs_1kb.fasta`                                                                        | Sequences of all contigs >= 1kb                    | 
| `viral_identification/identified_viral_contig_sequences/{sample_name}_identified_viral_contigs_final.fa`| Sequences of identified viral contigs              |
 	


