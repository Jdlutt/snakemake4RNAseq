# 指定配置文件
configfile: "config.yaml"

# 读取sample.txt第二列的样本信息列表
SAMPLES = []

with open('samples.txt', 'r') as file:
    for line in file:
        columns = line.strip().split('\t')
        SAMPLES.append(columns[1])

# 定义一个函数实现特异性文库参数的转换
def convert_hisat2_to_featurecounts(hisat2_strandness):
    if hisat2_strandness == "RF":
        return 2  # FeatureCounts: -s 2 (反向链)
    elif hisat2_strandness == "FR":
        return 1  # FeatureCounts: -s 1 (正向链)
    elif hisat2_strandness == "FF":
        return 1  # FeatureCounts: -s 1 (正向链)
    else:
        raise ValueError("Invalid HISAT2 strandness value")

# 定义合并表达矩阵的输出文件路径
def get_counts_file_path(MATRIX_HOME):
    return f"{MATRIX_HOME}.counts.matrix"

def get_tpm_file_path(MATRIX_HOME):
    return f"{MATRIX_HOME}.TMM.TPM.matrix"

rule all:
  input:
    expand("14.data/clean/{sample}_{i}.fastp.fastq.gz",sample = SAMPLES,i = [1,2]),
    expand("genome.hisat2.{i}.ht2",i = range(1,9)),
    expand("21.mapping/{sample}.bam",sample = SAMPLES),
    expand("22.quantification/{sample}.count",sample = SAMPLES),
    counts=get_counts_file_path(config["MATRIX_HOME"]),
    tpm=get_tpm_file_path(config["MATRIX_HOME"]),
    de_result="23.DE_analysis/de_result.tsv"

# 使用fastp对数据进行过滤和质控
rule fastp:
  input:
    r1=config["r1_RAWDATA"],
    r2=config["r2_RAWDATA"],
  output:
    r1="14.data/clean/{sample}_1.fastp.fastq.gz",
    r2="14.data/clean/{sample}_2.fastp.fastq.gz",
    html="14.data/clean/{sample}.fastp.html",
    json="14.data/clean/{sample}.fastp.json"
  threads:5
  log:"14.data/clean/{sample}.fastp.log"
  shell:
    """
      fastp -w {threads} \
      -i {input.r1} \
      -I {input.r2} \
      -o {output.r1} \
      -O {output.r2} \
      -h {output.html} \
      -j {output.json} \
      1>{log} 2>&1
    """

# 构建参考基因组index
rule build_hisat2_index:
  input:
    genome=config["genome"]
  output:
    idx=expand("genome.hisat2.{i}.ht2",i = range(1,9))
  log:"12.ref/hisat2-build.log"
  shell:
    """
      hisat2-build {input.genome} genome.hisat2 1>{log} 2>&1
    """

# 使用hisat2比对到参考基因组
rule align_hisat2:
  input:
    r1="14.data/clean/{sample}_1.fastp.fastq.gz",
    r2="14.data/clean/{sample}_2.fastp.fastq.gz",
    idx=expand("genome.hisat2.{i}.ht2",i = range(1,9))
  output:
    sam=temp("21.mapping/{sample}.sam")
  log:"21.mapping/{sample}.log"
  threads:4
  params:
    libType=config["libType"]
  shell:
    """
      hisat2 -p {threads} \
      -x genome.hisat2 \
      -1 {input.r1} \
      -2 {input.r2} \
      --new-summary \
      --rna-strandness {params.libType} \
      -S {output.sam} \
      2>{log}
    """

# 将sam文件转换为bam文件
rule sam2bam:
  input:
    sam="21.mapping/{sample}.sam",
  output:
    bam="21.mapping/{sample}.bam",
    index="21.mapping/{sample}.bam.bai"
  threads:4
  shell:
    """
      samtools sort -@ {threads} -o {output.bam} {input.sam}
      samtools index {output.bam}
    """

# 使用run-featurecounts表达定量
rule quanti_featurecounts:
  input:
    bam="21.mapping/{sample}.bam",
    gtf=config["gtf"]
  output:
    counts="22.quantification/{sample}.count"
  params:
    out_prefix="22.quantification/{sample}",
    libType=convert_hisat2_to_featurecounts(config["libType"])
  shell:
    """
      Rscript 11.software/RunFeatureCounts/run-featurecounts.R \
      -b {input.bam} \
      -g {input.gtf} \
      -s {params.libType} \
      -o {params.out_prefix}
    """

# 合并表达矩阵
rule abundance_estimate_to_matrix:
  input:
    counts_files=expand("22.quantification/{sample}.count",sample=SAMPLES)
  output:
    counts=get_counts_file_path(config["MATRIX_HOME"]),
    tpm=get_tpm_file_path(config["MATRIX_HOME"])
  params:
    MATRIX_HOME=config["MATRIX_HOME"],
    ids=config["ids"]
  shell:
    """
      ls {input.counts_files} > quant_files.txt
      perl 11.software/RunFeatureCounts/abundance_estimates_to_matrix.pl \
      --est_method featureCounts \
      --quant_files quant_files.txt \
      --out_prefix {params.MATRIX_HOME}

      sed -i '1s/^/{params.ids}/' {output.counts}
      sed -i '1s/^/{params.ids}/' {output.tpm}
    """

# 差异表达分析
rule DE_analysis:
  input:
    counts=get_counts_file_path(config["MATRIX_HOME"]),
    samples_file=config["samples_file"],
    contrasts_file=config["contrasts_file"]
  output:
    de_result="23.DE_analysis/de_result.tsv"
  params:
    TRINITY_HOME=config["TRINITY_HOME"],
    DE_METHOD=config["DE_METHOD"]
  shell:
    """
      perl {params.TRINITY_HOME}/Analysis/DifferentialExpression/run_DE_analysis.pl \
      --matrix {input.counts} \
      --method {params.DE_METHOD} \
      --samples_file {input.samples_file} \
      --contrasts {input.contrasts_file} \
      --output 23.DE_analysis

      awk 'FNR==1 && NR!=1{{next}}{{print}}' 23.DE_analysis/*DE_results > {output.de_result}

      snakemake --dag | dot -Tsvg > dag.svg #生成流程图
    """

