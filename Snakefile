import os
import pysam

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

shell.prefix("source %s/config.sh; " % SNAKEMAKE_DIR)

if config == {}:
    configfile: "%s/config.json" % SNAKEMAKE_DIR

BAM_LIST = config["BAM_LIST"]
REFERENCE = config["REFERENCE"]
GENDER_MAP = config["GENDER_MAP"]
GENSTRIP_PARAMETERS = config["GENSTRIP_PARAMETERS"]
METADATA_DIR = config["METADATA_DIR"]
LOG_DIR = config["LOG_DIR"]
REF_DIR = config["REF_INFO"][REFERENCE]["REF_DIR"]
REF_PREFIX = config["REF_INFO"][REFERENCE]["REF_PREFIX"]

VCF_FILE = config["VCF_FILE"]
VCF_GT = config["VCF_GT"]

if not os.path.exists("log"):
    os.makedirs("log")

VCF_FILES = {"del": "del_discovery/svdiscovery.dels.vcf", "cnv": "cnv_discovery/results/gs_cnv.genotypes.vcf.gz", "other": VCF_FILE}

localrules: all

rule all:
    input: expand("gs_{vcf}.gt.vcf", vcf = VCF_GT)

rule genstrip_genotyper:
    input: lambda wildcards: VCF_FILES[wildcards.vcf], 
           BAM_LIST
    output: "gs_{vcf}.gt.vcf"
    params: sge_opts = "-l mfree=8G -N gs_genotyper", type = "{vcf}"
    run:
        cmd = """java -Xmx4g -cp $classpath \
                org.broadinstitute.gatk.queue.QCommandLine \
                -S $SV_DIR/qscript/SVQScript.q \
                -S $SV_DIR/qscript/SVGenotyper.q \
                -cp $classpath \
                -gatk $SV_DIR/lib/gatk/GenomeAnalysisTK.jar \
                -configFile {SNAKEMAKE_DIR}/{GENSTRIP_PARAMETERS} \
                -R {REF_DIR}/{REF_PREFIX}.fasta \
                -vcf {input[0]} \
                -I {input[1]} \
                -O {output} \
                -runDirectory {SNAKEMAKE_DIR}/gt_{params.type} \
                -jobLogDir gs_genotyper_log \
                -jobRunner Drmaa \
                -gatkJobRunner Drmaa \
                -jobNative "-l mfree=8G" \
                -jobNative "-q all.q" \
                -jobNative " -V -cwd" \
                -jobNative \"-w n\" \
                -md gs_md \
                -bamFilesAreDisjoint true \
                -parallelJobs 50 \
                -run"""
        print(cmd)
        shell(cmd)

rule genstrip_del_discovery:
    input: BAM_LIST, "cnv_discovery/results/gs_cnv.genotypes.vcf.gz"
    output: "del_discovery/svdiscovery.dels.vcf"
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
                -runDirectory {SNAKEMAKE_DIR}/del_discovery \
                -jobLogDir del_discovery_log \
                -jobRunner Drmaa \
                -gatkJobRunner Drmaa \
                -jobNative "-l mfree=8G" \
                -jobNative "-q all.q" \
                -jobNative " -V -cwd" \
                -jobNative \"-w n\" \
                -md gs_md \
                -rmd {REF_DIR} \
                -bamFilesAreDisjoint true \
                -run"""
        print(cmd)
        shell(cmd)

rule genstrip_cnv_discovery:
    input: BAM_LIST, "GS_PREPROCESS_FINISHED"
    output: "cnv_discovery/results/gs_cnv.genotypes.vcf.gz"
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
                -runDirectory {SNAKEMAKE_DIR}/cnv_discovery \
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
                -run"""

        print(cmd)
        shell(cmd)
        shell("touch {output}")


rule genstrip_preprocess:
    input: BAM_LIST, "SAMPLES_PASS"
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
                -I {input[0]} \
                -jobLogDir gs_preprocess_log \
                -md gs_md \
                -bamFilesAreDisjoint true \
                -reduceInsertSizeDistributions true \
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

rule check_LB_tag:
    input: BAM_LIST
    output: "SAMPLES_PASS"
    params: sge_opts=""
    run:
        bams = []
        with open(input[0], "r") as manifest:
            for line in manifest:
                bams.append(line.rstrip())
        for file in bams:
            samfile = pysam.AlignmentFile(file)
            header_lines = samfile.text.split("\n")
            read_group = None
            for line in header_lines:
                if line.startswith("@RG"):
                    read_group = line.split("\t")
            if read_group is None:
                sys.exit("File %s missing @RG header line" % file)
            tags = map(lambda x: x.split(":")[0], read_group)

            for tag in ["ID", "SM", "LB"]:
                if not tag in tags:
                    sys.exit("File %s missing %s tag" % (file, tag))
        shell("touch {output[0]}")
