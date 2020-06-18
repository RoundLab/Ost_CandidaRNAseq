#!/bin/bash

# Add SBATCH directives for run with slurm scheduler here.

###### User Defined Variables ############################
# Kallisto reference index directory (index made in script)
REFSEQDIR=
mkdir -p $REFSEQDIR
# Host (mouse) transcriptome bowtie2 reference index directory (must be present)
HOSTREFSEQDIR=
# Raw fastq sequences directory (seqs must be present with file names as: *_R[12]_001.fastq.gz)
RAWSEQDIR=
# Scratch directory for temporary files
SCRATCH=
mkdir -p $SCRATCH
# Working directory to copy final alignments to
WRKDIR=
# Number of processes to run in parallel
NumProc=
#########################################################

# Load required software / document versions
module load cutadapt/1.14
module load fastqc/0.11.4
module load trim_galore/0.4.4
module load bowtie2/2-2.2.9
module load kallisto/0.45.0

# Step 1: Create reference index with kallisto
# 	 Note: we provide archived version link to explicitly record version, as "current" stable download links update with each version.
cd ${REFSEQDIR}
wget http://candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/archive/C_albicans_SC5314_version_A22-s07-m01-r99_default_coding.fasta.gz
kallisto index -i C_albicans_SC5314_version_A22-s07-m01-r99_default_coding C_albicans_SC5314_version_A22-s07-m01-r99_default_coding.fasta.gz

# Step 2: Adapter and Quality Trimming
mkdir -p ${SCRATCH}/trim_out
cd ${RAWSEQDIR}

mkdir -p ${SCRATCH}/trim_out
cd ${RAWSEQDIR}

ls -1 *R1_001.fastq.gz | cut -f 1-7 -d _ | parallel -j ${NumProc} 'trim_galore --paired --retain_unpaired --fastqc --length 20 -q 20 -o /scratch/general/nfs1/u0210816/17599R_CalbRNAseq/trim_out {}_R1_001.fastq.gz {}_R2_001.fastq.gz'

# Step 3: Filter mouse (host) reads

mkdir -p ${SCRATCH}/MmFilt

cd ${SCRATCH}/trim_out/

for f in *_R1_001_val_1.fq.gz
 do bowtie2 -x ${HOSTREFSEQDIR} -p ${NumProc} --very-fast -1 ${f} -2 ${f%_R1_001_val_1.fq.gz}_R2_001_val_2.fq.gz -S ../MmFilt/${f%_R1_001_val_1.fq.gz}.sam --un-conc-gz ../MmFilt/${f%%_*}_trimmed_unAlignPairs.fq.gz
done

# remove mouse aligned reads if not of interest
# rm MmFilt/*.sam

# Step 4: Run Kallisto on non-host reads with bootstrapping
mkdir -p ${SCRATCH}/kallisto_outputs_MmFilt
cd ${SCRATCH}/MmFilt

for QCREADS in *.fq.1.gz
 do kallisto quant -t ${NumProc} -b 100 -i ${REFSEQDIR}/C_albicans_SC5314_version_A22-s07-m01-r99_default_coding -o ../kallisto_outputs_MmFilt/${QCREADS%_trimmed_unAlignPairs.fq.1.gz} ${QCREADS} ${QCREADS%.fq.1.gz}.fq.2.gz
done

cp ${SCRATCH}/kallisto_outputs_MmFilt ${WRKDIR}