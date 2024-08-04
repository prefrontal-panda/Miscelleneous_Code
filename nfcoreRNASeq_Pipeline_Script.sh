#!/bin/bash

# Job partition
#SBATCH -p comp
# Job name
#SBATCH --job-name="nf-core/rnaseq"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=74GB
#SBATCH --account=xm41
#SBATCH --time=168:00:00

# Initialising mamba
. "/usr/local/miniforge3/23.3.1/conda/etc/profile.d/conda.sh"
. "/usr/local/miniforge3/23.3.1/conda/etc/profile.d/mamba.sh"
mamba init

# Activate Environment
mamba activate nextflow

# Loading modules
module load java/openjdk-21.0.1
module load star/2.7.10b

# RSEM path
# Path variable should be something like: PATH variable=/scratch/.../debbie_rnaseq/nfcore_trial/RSEM
export PATH=$PATH:/scratch/xm41/debbie_rnaseq/nfcore_trial/RSEM #"{where the RSEM is saved}"

# Files
# Make sure in your sample sheet, you specify the path to the files
assembly_file=/scratch/xm41/debbie_rnaseq/Rat_7.2.111_dna.sm.primary_assembly.fa
gtf_file=/scratch/xm41/debbie_rnaseq/Rattus_norvegicus.mRatBN7.2.111.gtf
sample_in=/scratch/xm41/debbie_rnaseq/PTE_samplesheet.csv #/scratch/xm41/debbie_rnaseq/GAERS_samplesheet.csv	

# Run
which nextflow
echo "Before nextflow"
nextflow run nf-core/rnaseq -r 3.14.0 \
	--input $sample_in \
	--outdir PTE_nfOut \
	--fasta $assembly_file \
	--gtf $gtf_file \
	--skip_dupradar \
	--skip_markduplicates \
	--save_reference \
	-resume