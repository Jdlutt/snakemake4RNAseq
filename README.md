# snakemake4RNAseq
Use snakemake to automate RNAseq analysis processes

## The software needed to complete the process 
fastp
https://github.com/OpenGene/fastp

hisat2
https://github.com/DaehwanKimLab/hisat2

samtools
https://github.com/samtools/samtools

**The above three software packages have been packaged into the environment**

RunFeatureCounts
http://git.genek.cn:3333/zhxd2/RunFeatureCounts.git

```bash
git clone http://git.genek.cn:3333/zhxd2/RunFeatureCounts.git
```

trinityrnaseq
https://github.com/trinityrnaseq/trinityrnaseq/releases

```bash
https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.15.1/trinityrnaseq-v2.15.1.FULL.tar.gz
tar -zxvf trinityrnaseq-v2.15.1.FULL.tar.gz
```

**You need to install the above two software manually**

## mamba installation and environment configuration

https://github.com/conda-forge/miniforge#mambaforge

Configure env
```bash
cd mambaforge/bin
./mamba init
source ~/.bashrc
```

copy env.yaml file to your workspace,then create the env
```bash
mamba env create -f env.yaml -n snakemake4rnaseq
```

## R-packages installation
make sure your R-packages including
* BiocManager
* argparser
* Rsubread
* limma
* DESeq2

download R-packages
```r
install.packages("BiocManager")
install.packages("argparser")
BiocManager::install("Rsubread")
BiocManager::install("limma")
BiocManager::install("DESeq2")
```

