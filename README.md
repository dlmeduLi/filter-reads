# Alignment reads: parent-of-origin filtering

Filter out parent-of-origin alignment reads according to SNP databases.

## Install

`filter_reads.py` requires `pysam` library.

You can install `pysam`:

```bash
sudo yum install python-pip
pip install pysam
```

## Usage

`filter_reads.py in.sam snp1.csv snp2.csv`

SNP data files `*.csv` have the following format:
```csv
chr,pos,ref,alt
1,3000248,G,T
1,3000289,T,G
1,3000353,C,T
1,3000355,T,C
1,3000444,T,A
1,3000906,G,T
```
chr: chrome name without `chr` prefix
pos: SNP position
ref: Reference base 
alt: Alternative base

## *.mup files

*.mup files illustrate the underline algorithms the script uses, they can be opened with [https://www.mindmup.com](https://www.mindmup.com) an on-line mind storm editor.


