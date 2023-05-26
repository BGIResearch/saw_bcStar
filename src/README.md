# Introduction
bcStar implement these steps in RNA-seq pipeline:
1. CID mapping
2. adapter sequences filter
3. align to reference

# Install
## Platform & Environment
* centos-7.0+
* gcc-5.2.0
* GNU Make 3.82

## Depended libraries
* hdf5-1.12.1
* deflate
* folly
* pthread
* zlib

## Compile
```shell
git clone git@github.com:BGIResearch/SAW_Mapping.git
cd SAW_Mapping
make
```
## Run
```text
Program: bcStar
Version: 1.0.4
Contact: GongChun<gongchun@genomics.cn>
Usage: bcStar  [options]... --genomeDir REFERENCE   --readFilesIn R1.fq R2.fq --bcPara <barcodeMapping parameter file>
Parameters of barcode mapping should be written in a file as the value of parameter --bcPara
e.g. echo --in=input.gz >bcp.txt, and run bcStar [other options] --bcPara bcp.txt
To list parameters of barcode mapping, run bcStar --helpBcPara
To list all parameters, run bcStar --help
```
## Results
* standard output : stat of barcode mapping
* *.bc.txt : barcoded reads count in each coordinate
* other files : output files of STAR