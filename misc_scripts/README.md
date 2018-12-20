Two indels were genotyped semi-manually, as alignment issues in repetitive regions caused them to be marked at several different adjacent positions in the bam cigar strings.

### Step 1: extract data from region

    # 809653 CEBPA deletion
    region=19:33301989-33301990
    samtools view -hb ${bam} ${region} > 809653_CEBPA.bam
    samtools index 809653_CEBPA.bam

    # 721214 & 548327 NPM1 insertions
    region=5:171410538-171410546
    samtools view -hb ${bam} ${region} > 721214_NPM1.bam
    samtools index 721214_NPM1.bam

### Step2: run python script

     python fetch_manual_ins_del.py 809653_CEBPA.bam CEBPA 809_barcodes.tsv 809653 > 809653.barcodes.tsv


### Tool Usage
```{shell}
python3 fetch_manual_ins_del.py -h
usage: fetch_manual_ins_del.py [-h] bam_file gene barcodes upn

Parse CB barcodes from Single cell rna seq data

positional arguments:
  bam_file    BAM file
  gene        NPM1 or CEBPA
  barcodes    list of good barcodes file
  upn         upn/sample name: will be used as prefix for out_file

optional arguments:
  -h, --help  show this help message and exit
```
