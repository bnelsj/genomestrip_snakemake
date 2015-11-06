import os

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

shell.prefix("source %s/config.sh; " % SNAKEMAKE_DIR)

if config == {}:
    configfile: "%s/config.json" % SNAKEMAKE_DIR

BAM_LIST = config["BAM_LIST"]
REFERENCE = config["REFERENCE"]
GENSTRIP_PARAMETERS = config["GENSTRIP_PARAMETERS"]
METADATA_DIR = config["METADATA_DIR"]
LOG_DIR = config["LOG_DIR"]
REF_DIR = config["REF_INFO"][REFERENCE]["REF_DIR"]
REF_PREFIX = config["REF_INFO"][REFERENCE]["REF_PREFIX"]

VCF_FILE = "results/gs_cnv.genotypes.vcf.gz"

if not os.path.exists("log"):
    os.makedirs("log")

localrules: all

rule all:
    input: "GS_PREPROCESS_FINISHED"

rule genstrip_genotyper:
    input: VCF_FILE, BAM_LIST
    output: "gs_genotypes.vcf"
    params: sge_opts = "-l mfree=8G -N gs_del_discovery"
    run:
        cmd = """java -Xmx4g -cp $classpath \
                org.broadinstitute.gatk.queue.QCommandLine \
                -S $SV_DIR/qscript/SVQScript.q \
                -S $SV_DIR/qscript/SVGenotyper.q \
                -cp $classpath \
                -gatk $SV_DIR/lib/gatk/GenomeAnalysisTK.jar \
                -configFile {SNAKEMAKE_DIR}/{GENSTRIP_PARAMETERS} \
                -R {REF_DIR}/{REF_PREFIX}.fasta \
                -genderMapFile {GENDER_MAP} \
                -vcf {input[0]} \
                -I {input[1]} \
                -O {output} \
                -runDirectory {SNAKEMAKE_DIR} \
                -jobLogDir gs_genotyper_log \
                -jobRunner Drmaa \
                -gatkJobRunner Drmaa \
                -jobNative "-l mfree=8G" \
                -jobNative "-q all.q" \
                -jobNative " -V -cwd" \
                -jobNative \"-w n\" \
                -jobNative \"-l L3cache=30M\" \
                -md gs_md \
                -bamFilesAreDisjoint false \
                -parallelJobs 50 \
                -run"""
        print(cmd)
        shell(cmd)

rule genstrip_del_discovery:
    input: BAM_LIST, "GS_CNV_DISCOVERY_FINISHED"
    output: "svdiscovery.dels.vcf"
    params: sge_opts = "-l mfree=8G -N gs_del_discovery"
    run:
        cmd = """java -Xmx4g -cp $classpath \
                org.broadinstitute.gatk.queue.QCommandLine \
                -S $SV_DIR/qscript/SVQScript.q \
                -S $SV_DIR/qscript/SVDiscovery.q \
                -cp $classpath \
                -gatk $SV_DIR/lib/gatk/GenomeAnalysisTK.jar \
                -configFile {SNAKEMAKE_DIR}/{GENSTRIP_PARAMETERS} \
                -R {REF_DIR}/{REF_PREFIX}.fasta \
                -ploidyMapFile {REF_DIR}/{REF_PREFIX}.ploidymap.txt \
                -I {input[0]} \
                -O {output} \
                -minimumSize 100 \
                -maximumSize 100000 \
                -runDirectory {SNAKEMAKE_DIR} \
                -jobLogDir gs_del_discovery_log \
                -jobRunner Drmaa \
                -gatkJobRunner Drmaa \
                -jobNative "-l mfree=8G" \
                -jobNative "-q all.q" \
                -jobNative " -V -cwd" \
                -jobNative \"-w n\" \
                -jobNative \"-l L3cache=30M\" \
                -md gs_md \
                -bamFilesAreDisjoint false \
                -run"""
        print(cmd)
        shell(cmd)

rule genstrip_cnv_discovery:
    input: BAM_LIST, "GS_PREPROCESS_FINISHED"
    output: "GS_CNV_DISCOVERY_FINISHED"
    params: sge_opts = "-l mfree=8G -N gs_cnv_discovery"
    run:
        cmd = """java -Xmx4g -cp $classpath \
                org.broadinstitute.gatk.queue.QCommandLine \
                -S $SV_DIR/qscript/discovery/cnv/CNVDiscoveryPipeline.q \
                -S $SV_DIR/qscript/SVQScript.q \
                -cp $classpath \
                -gatk $SV_DIR/lib/gatk/GenomeAnalysisTK.jar \
                -configFile {SNAKEMAKE_DIR}/{GENSTRIP_PARAMETERS} \
                -R {REF_DIR}/{REF_PREFIX}.fasta \
                -I {input[0]} \
                -md gs_md \
                -runDirectory {SNAKEMAKE_DIR} \
                -jobLogDir gs_cnv_discovery_log \
                -intervalList {REF_DIR}/{REF_PREFIX}.interval.list \
                -tilingWindowSize 1000 \
                -tilingWindowOverlap 500 \
                -maximumReferenceGapLength 1000 \
                -boundaryPrecision 100 \
                -minimumRefinedLength 500 \
                -jobRunner Drmaa \
                -gatkJobRunner Drmaa \
                -jobNative "-l mfree=20G" \
                -jobNative "-V -cwd" \
                -jobNative "-q all.q" \
                -jobNative \"-w n\" \
                -jobNative \"-l L3cache=30M\" \
                -run"""

        print(cmd)
        shell(cmd)
        shell("touch {output}")


rule genstrip_preprocess:
    input: BAM_LIST
    output: "GS_PREPROCESS_FINISHED"
    params: sge_opts = "-l mfree=8G -N gs_preproc"
    run:
        cmd = """java -Xmx4g -cp $classpath \
                org.broadinstitute.gatk.queue.QCommandLine \
                -S $SV_DIR/qscript/SVQScript.q \
                -S $SV_DIR/qscript/SVPreprocess.q \
                -cp $classpath \
                -gatk $SV_DIR/lib/gatk/GenomeAnalysisTK.jar \
                -configFile {SNAKEMAKE_DIR}/{GENSTRIP_PARAMETERS} \
                -R {REF_DIR}/{REF_PREFIX}.fasta \
                -ploidyMapFile {REF_DIR}/{REF_PREFIX}.ploidymap.txt \
                -I {input} \
                -jobLogDir gs_preprocess_log \
                -md gs_md \
                -bamFilesAreDisjoint false \
                -jobRunner Drmaa \
                -gatkJobRunner Drmaa \
                -jobQueue all.q \
                -jobNative \"-l mfree=8G\" \
                -jobNative \"-q all.q\" \
                -jobNative \"-V -cwd\" \
                -jobNative \"-w n\" \
                -jobNative \"-l L3cache=30M\" \
                -run """

        print(cmd)
        shell(cmd)
        shell("touch {output}")
