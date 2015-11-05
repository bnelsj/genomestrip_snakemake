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


if not os.path.exists("log"):
    os.makedirs("log")

localrules: all

rule all:
    input: "GS_PREPROCESS_FINISHED"

rule genstrip_preprocess:
    input: BAM_LIST
    output: "GS_PREPROCESS_FINISHED"
    params: sge_opts = "-l mfree=8G -N gs_preproc"
    run:
        shell(
                "java -Xmx4g -cp $classpath "
                "org.broadinstitute.gatk.queue.QCommandLine "
                "-S $SV_DIR/qscript/SVQScript.q "
                "-S $SV_DIR/qscript/SVPreprocess.q "
                "-cp $classpath "
                "-gatk $SV_DIR/lib/gatk/GenomeAnalysisTK.jar "
                "-configFile {SNAKEMAKE_DIR}/{GENSTRIP_PARAMETERS} "
                "-R {REF_DIR}/{REF_PREFIX}.fasta "
                "-genomeMaskFile {REF_DIR}/{REF_PREFIX}.svmask.fasta "
                "-ploidyMapFile {REF_DIR}/{REF_PREFIX}.ploidymap.txt "
                "-I {input} "
                "-jobLogDir $TMPDIR/gs_preprocess_log "
                "-md $TMPDIR/gs_md "
                "-bamFilesAreDisjoint false "
                "-jobRunner Drmaa "
                "-gatkJobRunner Drmaa "
                "-jobQueue all.q "
                "-jobNative \"-l mfree=8G\" "
                "-jobNative \"-q all.q\" "
                "-jobNative \"-V -cwd\" "
                "-jobNative \"-l L3cache=30M\" "
                "-run "
                )
        shell("rsync $TMPDIR .")
        shell("touch {output}")
