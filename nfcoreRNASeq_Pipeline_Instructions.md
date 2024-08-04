# This file is to go through the steps of transcriptomics analysis using bash
# For this, we will be using the nf-core/rnaseq pipeline
# This is following the tutorial over at: https://github.com/mucosal-immunology-lab/RNAseq_NFCORE 

# First, log in to server and change to your directory
# e.g. cd /scratch/.../debbie_rnaseq

# Then, download the reads and gtf files
# DNA reads obtained from: https://ftp.ensembl.org/pub/release-112/fasta/rattus_norvegicus/dna/
wget -r --no-parent -A 'Rattus_norvegicus.mRatBN7.2.dna_sm.primary_assembly.*.fa.gz' https://ftp.ensembl.org/pub/release-112/fasta/rattus_norvegicus/dna/
# Combining all the files
cat *.fa.gz > Rat_7.2.111_dna.sm.primary_assembly.fa.gz
gunzip Rat_7.2.111_dna.sm.primary_assembly.fa.gz

# GTF file obtained from: https://ftp.ensembl.org/pub/release-111/gtf/rattus_norvegicus/
wget https://ftp.ensembl.org/pub/release-111/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.111.gtf.gz
guunzip Rattus_norvegicus.mRatBN7.2.111.gtf.gz

# Now we prepare to run the pipeline
# First, we load mamba
module avail miniforge3 # Checking for the module availability
module load miniforge3/23.3.1 # Loading
# Creating environment
mamba create -n nextflow nextflow \
    salmon=1.10.0 fq fastqc umi_tools \
    trim-galore bbmap sortmerna samtools \
    picard stringtie bedtools rseqc \
    qualimap preseq multiqc subread \
    ucsc-bedgraphtobigwig ucsc-bedclip \
    bioconductor-deseq2 bioconductor-tximport \
    bioconductor-tximeta r-pheatmap -c bioconda -c conda-forge -c default \
    --skip_dupradar 

# Activate environment
mamba activate nextflow
# Deactivate environment
# mamba deactivate
# Environment location: /home/debbiec/.../tools/multiomics/debbie.chong/miniconda/conda/envs/nextflow

# now we need to download and RSEM for the estimation of gene and isoform expression levels from RNA-Seq data.
git clone https://github.com/deweylab/RSEM # download
cd RSEM; make # compiling

# Run in a separate smux session
smux n --time=3-05:00:00 --mem=74GB --ntasks=1 --cpuspertask=12 -J nf-core/rnaseq

# Load JAVA AND STAR before running the pipeline
module load java/openjdk-17.0.2
module load star/2.7.10b

# Testing the pipeline before running the actual thing
nextflow run nf-core/rnaseq -r 3.14.0 -profile test --outdir test -resume --skip_dupradar --skip_markduplicates

# Run the file
bash data.pipeline.sh