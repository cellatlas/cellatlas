# Universal preprocessing of single-cell genomics data

The cellatlas tool uses seqspec and kallisto bustools to facilitate universal preprocessing of single-cell genomics data.


## Installation
The cellatlas command-line tool can be installed with pip:

```bash
pip install git+https://github.com/cellatlas/cellatlas.git
```

and can be run with 
```bash
cellatlas build \
-o out \
-m modality \
-s spec.yaml \
-fa http://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
-g http://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz \
-fb feature_barcodes.txt \
fastqs/R1.fastq.gz fastqs/R2.fastq.gz",
```
- `-o` is the output folder
- `-m` is the modality (rna/atac etc)
- `-s` is the `seqspec` for the supplied FASTQs
- `-fa` is either a link to the genome FASTA or the file itself
- `-g` is either a link to the genome annotation (GTF) or the file itself
- `-fb` is optional and is the feature barcode file for tag/protein/crispr assays

FASTQs are supplied at the end of the command in any order.
