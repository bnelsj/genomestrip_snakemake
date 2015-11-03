module purge
module load modules modules-init modules-gs/prod modules-eichler/prod
module load java/7u17 R/2.15.0 tabix/0.2.6 samtools/1.2 genomestrip

export classpath="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar"
export PROJECT_NAME=testing
export JOB_QUEUE=all.q
export RUN_DIR=`pwd`
export JOB_LOG_DIR=genomestrip_log
export JOB_RUNNER=Drmaa
export REF_DIR=/net/eichler/vol2/eee_shared/assemblies/human_1kg_v37/genomestrip_metadata/1000G_phase1
export REF_FASTA=human_g1k_v37.fasta
export CONFIG_FILE=genstrip_parameters.txt
export INPUT_BAM_FILES=hgdp_bam_file.list
export METADATA_DIR=genomestrip_metadata
