# snakemake4RNAseq
Use snakemake to automate RNAseq analysis processes

## The software needed to complete the process 
fastp
https://github.com/OpenGene/fastp

hisat2
https://github.com/DaehwanKimLab/hisat2

samtools
https://github.com/samtools/samtools

RunFeatureCounts
http://git.genek.cn:3333/zhxd2/RunFeatureCounts.git

trinityrnaseq
https://github.com/trinityrnaseq/trinityrnaseq

All software is packaged into the mamba environment or software directory

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

