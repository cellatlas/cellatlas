---
title: Preprocess reads
date: 2024-07-07
authors:
  - title: A. Sina Booeshaghi
---

Preprocessing data is the procedure where

1. Sequencing reads are aligned to a reference
2. Barcodes errors are corrected
3. UMIs/reads are counted

The goal is to produce a count matrix, where rows are cells or samples and columns are biological features such as genes, proteins, or genomic regions.

There are many tools that perform single-cell RNA-sequencing preprocessing. For this tutorial we will use the `cellatlas` infrastructure which couples the `kb-python` (which uses `kallisto` and `bustools`) and `seqspec` to perform quantification of any modality. Alignment is performed on the sequencing reads with `kallisto` and barcode correction and UMI counting is performed with `bustools`. With `cellatlas` we can effectively preprocess all kinds of single-cell data.

We need to install `kb-python`, `cellatlas`, `seqspec`, and a few helper tools:

```bash
# install with pip
pip install cellatlas kb-python seqspec gget # cellatlas has dependencies that must be installed

# cellatlas
cellatlas --version

# seqspec
seqspec --version

# verify installation
kb --version
```

`seqspec` works seamlessly with `kb-python` to preprocess sequencing reads from any modality.

Throughout this tutorial we will use the `dogmaseq-dig` dataset which is a multimodal assay (RNA/ATAC/PROTEIN/TAG).

To understand how each tool works, please review their manuscript:

- seqspec
- cellatlas
- gget
- kallisto
- bustools
- BUS file

# RNA data

Examples include (list specific papers corresponding to specs)

- scRNAseq
- snRNAseq

**Inputs**

- Fastqs (must be represented as regions in the `seqspec`)
  - `fastqs/rna_R1_SRR18677638.fastq.gz`
  - `fastqs/rna_R2_SRR18677638.fastq.gz`
- Links to references
  - Human Genome: `http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`
  - Human Annotation: `http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz`
- validated seqspec (and associated files)
  - `dogmaseq-dig/spec.yaml`
  - `dogmaseq-dig/RNA-737K-arc-v1.txt`

**Outputs:**

- Reference files
  - `kallisto` alignment transcriptome: `rna/transcriptome.fa`
  - `kallisto` alignment reference `rna/index.idx` (built from `transcriptome.fa`)
  - `bustools` transcripts to genes map: `rna/t2g.txt`
- QC files
  - `cellatlas` run info: `cellatlas_info.json`
  - `kb_python` run info: `kb_info.json`
  - `kallisto` alignment info: `run_info.json`
  - `bustools` quantification info: `inspect.json`
- Count matrix
  - AnnData file`rna/counts_unfiltered/adata.h5ad` built from
    - Barcodes: `rna/counts_unfiltered/cells_x_genes.barcodes.txt`
    - Genes: `rna/counts_unfiltered/cells_x_genes.genes.txt`
    - Matrix: `rna/counts_unfiltered/cells_x_genes.mtx`
- Alignment files
  - Alignment BUS file: `rna/output.bus`
  - Error-corrected & sorted BUS file: `rna/output.unfiltered.bus`
  - Equivalence class to transcripts map: `rna/matrix.ec`
  - Transcripts: `rna/transcripts.txt`

## Preprocessing

First we download the links to the references

```bash
$ gget ref -w dna,gtf homo_sapiens
Sun Oct 15 05:02:38 2023 INFO Fetching reference information for homo_sapiens from Ensembl release: 110.
{
    "homo_sapiens": {
        "genome_dna": {
            "ftp": "http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
            "ensembl_release": 110,
            "release_date": "2023-04-21",
            "release_time": "17:28",
            "bytes": "841M"
        },
        "annotation_gtf": {
            "ftp": "http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz",
            "ensembl_release": 110,
            "release_date": "2023-04-24",
            "release_time": "11:54",
            "bytes": "52M"
        }
    }
}
```

We can preprocess our data with either `cellatlas`, `kb-python`, or `kallisto` and `bustools`. While each tool is distinct, they all work together to reduce the number of steps to preprocess data. You can think of them as existing in a hierarchy with `cellatlas` living at the top with the least complex usage, and `kallisto` and `bustools` living at the bottom.

### `cellatlas`

```bash
cellatlas build \
-o rna -m rna -s spec.yaml \
-fa http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
-g http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz -fb '' \
fastqs/rna_R1_SRR18677638.fastq.gz fastqs/rna_R2_SRR18677638.fastq.gz
```

### `kb-python`

`cellatlas` generates the following `kb-python` code to preprocess the data.

```bash
kb ref \
-i rna/index.idx \
-g rna/t2g.txt \
-f1 rna/transcriptome.fa \
http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

kb count \
-i rna/index.idx \
-g rna/t2g.txt \
-x 0,0,16:0,16,28:1,0,102 \ # this can also be generated with seqspec
-w RNA-737K-arc-v1.txt \
-o rna --h5ad -t 2 \
fastqs/rna_R1_SRR18677638.fastq.gz fastqs/rna_R2_SRR18677638.fastq.gz
```

### `kallisto` and `bustools`

`kb-python` generates the following `kallisto` and `bustools` code to run.

```bash
# reference building kb index uses custom code to download and slice the genome with the relevant features (extracted from the gtf)
kallisto index -i rna/index.idx rna/transcriptome.fa

# alignment with kallisto
kallisto bus -i rna/index.idx -o rna -x 0,0,16:0,16,28:1,0,102 -t 2 fastqs/rna_R1_SRR18677637.fastq.gz fastqs/rna_R2_SRR18677637.fastq.gz fastqs/rna_R1_SRR18677638.fastq.gz fastqs/rna_R2_SRR18677638.fastq.gz

# error correction and counting with bustools
bustools inspect rna/output.bus
bustools sort -o rna/tmp/output.s.bus -T rna/tmp -t 2 -m 4G rna/output.bus
bustools inspect rna/tmp/output.s.bus
bustools inspect -o rna/inspect.json -w RNA-737K-arc-v1.txt rna/tmp/output.s.bus
bustools correct -o rna/tmp/output.s.c.bus -w RNA-737K-arc-v1.txt rna/tmp/output.s.bus
bustools inspect rna/tmp/output.s.c.bus
bustools sort -o rna/output.unfiltered.bus -T rna/tmp -t 2 -m 4G rna/tmp/output.s.c.bus
bustools inspect rna/output.unfiltered.bus
bustools count -o rna/counts_unfiltered/cells_x_genes -g rna/t2g.txt -e rna/matrix.ec -t rna/transcripts.txt --genecounts rna/output.unfiltered.bus
```

Note that in the `kb-python` and `kallisto` `bustools` examples above the `-x` and the `-w` arguments can also be generated directly with `seqspec`:

```bash
$ seqspec index \
-t kb -m rna \
-r fastqs/rna_R1_SRR18677638.fastq.gz,rna_R2_SRR18677638.fastq.gz \
spec.yaml
0,0,16:0,16,28:1,0,102

$ seqspec onlist -m rna -r barcode spec.yaml
/path/to/the/relevant/onlist/dogmaseq-dig/RNA-737K-arc-v1.txt
```

# TAG data

Examples include

- Cell Hashing
- Other tags?

**Inputs**

- Fastqs (must be represented as regions in the `seqspec`)
  - `fastqs/tag_R1_SRR18677640.fastq.gz`
  - `fastqs/tag_R2_SRR18677640.fastq.gz`
- Links to references
  - Human Genome: `http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`
  - Human Annotation: `http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz`
  - Tag barcode file: `tag_0419_feature_barcodes.txt`
- validated seqspec (and associated files)
  - `spec.yaml`
  - `RNA-737K-arc-v1.txt`

**Outputs:**

- Reference files
  - `kallisto` alignment "transcriptome" (this is the kITE barcode mismatch file): `tag/transcriptome.fa`
  - `kallisto` alignment reference `tag/index.idx` (built from `transcriptome.fa`)
  - `bustools` barcode mismatch to barcode map: `tag/t2g.txt`
- QC files
  - `cellatlas` run info: `cellatlas_info.json`
  - `kb_python` run info: `kb_info.json`
  - `kallisto` alignment info: `run_info.json`
  - `bustools` quantification info: `inspect.json`
- Count matrix
  - AnnData file `tag/counts_unfiltered/adata.h5ad` built from
    - Barcodes: `tag/counts_unfiltered/cells_x_features.barcodes.txt`
    - Tags: `tag/counts_unfiltered/cells_x_features.genes.txt`
    - Matrix: `tag/counts_unfiltered/cells_x_features.mtx`
- Alignment files
  - Alignment BUS file: `tag/output.bus`
  - Error-corrected & sorted BUS file: `tag/output.unfiltered.bus`
  - Equivalence class to transcripts map: `tag/matrix.ec`
  - Transcripts: `tag/transcripts.txt`

## Preprocessing

Insert introduction here

### `cellatlas`

```bash
cellatlas build \
-o tag \
-m tag \
-s spec.yaml \
-fa http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
-g http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz \
-fb tag_0419_feature_barcodes.txt \
fastqs/tag_R1_SRR18677640.fastq.gz fastqs/tag_R2_SRR18677640.fastq.gz
```

### `kb-python`

`cellatlas` generates the following `kb-python` code to preprocess the data.

```bash
# build alignment reference
kb ref \
--workflow kite \
-i tag/index.idx \
-g tag/t2g.txt \
-f1 tag/transcriptome.fa \
tag_0419_feature_barcodes.txt

# perform alignment, error correction, and counting
kb count \
--workflow kite \
-i tag/index.idx \
-g tag/t2g.txt \
-x 0,0,16:0,16,28:1,0,15 \
-w RNA-737K-arc-v1.txt \
-o tag --h5ad -t 2 \
fastqs/tag_R1_SRR18677640.fastq.gz fastqs/tag_R2_SRR18677640.fastq.gz
```

### `kallisto` and `bustools`

`kb-python` generates the following `kallisto` and `bustools` code to run.

```bash
# reference building kb index uses custom code to create the "barcode mismatch" kITE index from the tag barcode file (the naming convention of "transciptome.fa" is standard despite the file containing synthetic barcodes to be indexed)
kallisto index -i tag/index.idx tag/transcriptome.fa

# alignment with kallisto
kallisto bus -i tag/index.idx -o tag -x 0,0,16:0,16,28:1,0,15 -t 2 fastqs/tag_R1_SRR18677640.fastq.gz fastqs/tag_R2_SRR18677640.fastq.gz

# error correction and counting with bustools
bustools inspect tag/output.bus
bustools sort -o tag/tmp/output.s.bus -T tag/tmp -t 2 -m 4G tag/output.bus
bustools inspect tag/tmp/output.s.bus
bustools inspect -o tag/inspect.json -w RNA-737K-arc-v1.txt tag/tmp/output.s.bus
bustools correct -o tag/tmp/output.s.c.bus -w RNA-737K-arc-v1.txt tag/tmp/output.s.bus
bustools inspect tag/tmp/output.s.c.bus
bustools sort -o tag/output.unfiltered.bus -T tag/tmp -t 2 -m 4G tag/tmp/output.s.c.bus
bustools inspect tag/output.unfiltered.bus
bustools count -o tag/counts_unfiltered/cells_x_features -g tag/t2g.txt -e tag/matrix.ec -t tag/transcripts.txt --genecounts tag/output.unfiltered.bus
```

Note that in the `kb-python` and `kallisto` `bustools` examples above the `-x` and the `-w` arguments can also be generated directly with `seqspec`:

```bash
$ seqspec index \
-t kb \
-m tag \
-r fastqs/tag_R1_SRR18677640.fastq.gz,fastqs/tag_R2_SRR18677640.fastq.gz \
spec.yaml
0,0,16:0,16,28:1,0,15

$ seqspec onlist -m tag -r barcode spec.yaml
/path/to/the/relevant/onlist/dogmaseq-dig/RNA-737K-arc-v1.txt
```

# Protein data

Example assays include:

- CITEseq
- MultiSeq
- ClickTags

**Inputs**

- Fastqs (must be represented as regions in the `seqspec`)
  - `fastqs/protein_R1_SRR18677644.fastq.gz`
  - `fastqs/protein_R2_SRR18677644.fastq.gz`
- Links to references
  - Human Genome: `http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`
  - Human Annotation: `http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz`
  - Protein barcode file: `protein_feature_barcodes.txt`
- Validated seqspec (and associated files)
  - `spec.yaml`
  - `RNA-737K-arc-v1.txt`

**Outputs:**

- Reference files
  - `kallisto` alignment "transcriptome" (this is the kITE barcode mismatch file): `protein/transcriptome.fa`
  - `kallisto` alignment reference `tag/index.idx` (built from `transcriptome.fa`)
  - `bustools` barcode mismatch to barcode map: `protein/t2g.txt`
- QC files
  - `cellatlas` run info: `cellatlas_info.json`
  - `kb_python` run info: `kb_info.json`
  - `kallisto` alignment info: `run_info.json`
  - `bustools` quantification info: `inspect.json`
- Count matrix
  - AnnData file `protein/counts_unfiltered/adata.h5ad` built from
    - Barcodes: `protein/counts_unfiltered/cells_x_features.barcodes.txt`
    - Tags: `protein/counts_unfiltered/cells_x_features.genes.txt`
    - Matrix: `protein/counts_unfiltered/cells_x_features.mtx`
- Alignment files
  - Alignment BUS file: `protein/output.bus`
  - Error-corrected & sorted BUS file: `protein/output.unfiltered.bus`
  - Equivalence class to transcripts map: `protein/matrix.ec`
  - Transcripts: `protein/transcripts.txt`

## Preprocessing

TODO : Insert introduction here

### `cellatlas`

```bash
cellatlas build \
-o protein \
-m protein \
-s spec.yaml \
-fa http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
-g http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz \
-fb protein_feature_barcodes.txt \
fastqs/protein_R1_SRR18677644.fastq.gz fastqs/protein_R2_SRR18677644.fastq.gz
```

### `kb-python`

`cellatlas` generates the following `kb-python` code to preprocess the data.

```bash
# build alignment reference
kb ref \
--workflow kite \
-i protein/index.idx \
-g protein/t2g.txt \
-f1 protein/transcriptome.fa \
protein_feature_barcodes.txt

# perform alignment, error correction, and counting
kb count \
--workflow kite \
-i protein/index.idx \
-g protein/t2g.txt \
-x 0,0,16:0,16,28:1,0,15 \
-w RNA-737K-arc-v1.txt \
-o protein --h5ad -t 2 \
fastqs/protein_R1_SRR18677644.fastq.gz fastqs/protein_R2_SRR18677644.fastq.gz
```

### `kallisto` and `bustools`

`kb-python` generates the following `kallisto` and `bustools` code to run.

```bash
# reference building kb index uses custom code to create the "barcode mismatch" kITE index from the tag barcode file (the naming convention of "transciptome.fa" is standard despite the file containing synthetic barcodes to be indexed)
kallisto index -i protein/index.idx protein/transcriptome.fa

# alignment with kallisto
kallisto bus -i protein/index.idx -o protein -x 0,0,16:0,16,28:1,0,15 -t 2 fastqs/protein_R1_SRR18677644.fastq.gz fastqs/protein_R2_SRR18677644.fastq.gz

# error correction and counting with bustools
bustools inspect tag/output.bus
bustools sort -o tag/tmp/output.s.bus -T tag/tmp -t 2 -m 4G tag/output.bus
bustools inspect tag/tmp/output.s.bus
bustools inspect -o tag/inspect.json -w RNA-737K-arc-v1.txt tag/tmp/output.s.bus
bustools correct -o tag/tmp/output.s.c.bus -w RNA-737K-arc-v1.txt tag/tmp/output.s.bus
bustools inspect tag/tmp/output.s.c.bus
bustools sort -o tag/output.unfiltered.bus -T tag/tmp -t 2 -m 4G tag/tmp/output.s.c.bus
bustools inspect tag/output.unfiltered.bus
bustools count -o tag/counts_unfiltered/cells_x_features -g tag/t2g.txt -e tag/matrix.ec -t tag/transcripts.txt --genecounts tag/output.unfiltered.bus
```

Note that in the `kb-python` and `kallisto` `bustools` examples above the `-x` and the `-w` arguments can also be generated directly with `seqspec`:

```bash
$ seqspec index \
-t kb \
-m protein \
-r fastqs/protein_R1_SRR18677644.fastq.gz fastqs/protein_R2_SRR18677644.fastq.gz \
spec.yaml
0,0,16:0,16,28:1,0,15

$ seqspec onlist -m protein -r barcode spec.yaml
/path/to/the/relevant/onlist/dogmaseq-dig/RNA-737K-arc-v1.txt
```

# CRISPR data

Example assays include

- PerturbSeq
- TAPseq

Note that single-cell CRISPR guide RNAs can be quantified in the same way as TAG and PROTEIN data. Simply supply the guide RNA barcode file as the "feature barcodes" file.

# ATAC data

Examples include

- snATACseq
- 10xMultiome

**Inputs**

- Fastqs (must be represented as regions in the `seqspec`)
  - `fastqs/atac_R1_SRR18677642.fastq.gz`
  - `fastqs/atac_R2_SRR18677642.fastq.gz`
  - `fastqs/atac_R3_SRR18677642.fastq.gz`
- Links to references
  - Human Genome: `http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`
  - Human Annotation: `http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz`
- Validated seqspec (and associated files)
  - `spec.yaml`
  - `ATA-737K-arc-v1_rc.txt`

**Outputs:**

- Reference files
  - `kallisto` alignment "transcriptome" (this is the set of putative "regions/peaks" found by peak calling): `atac/peaks.fa`
  - `kallisto` alignment reference `atac/index.idx` (built from `peaks.fa`)
  - `bustools` peak to peak map: `atac/t2g.txt`
- QC files
  - `cellatlas` run info: `cellatlas_info.json`
  - `kb_python` run info: `kb_info.json`
  - `kallisto` alignment info: `run_info.json`
  - `bustools` quantification info: `inspect.json`
- Count matrix
  - AnnData file `atac/counts_unfiltered/adata.h5ad` built from
    - Barcodes: `atac/counts_mult/cells_x_genes.barcodes.txt`
    - Regions: `atac/counts_mult/cells_x_genes.genes.txt`
    - Matrix: `atac/counts_mult/cells_x_genes.mtx`
- Alignment files
  - Alignment BUS file: `atac/output.bus`
  - Error-corrected & sorted BUS file: `atac/output.unfiltered.bus`
  - Equivalence class to transcripts map: `atac/matrix.ec`
  - Transcripts: `atac/transcripts.txt`

## Preprocessing

todo: Insert intro here

### `cellatlas`

```bash
cellatlas build \
-o atac \
-m atac \
-s spec.yaml \
-fa http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
-g http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz \
-fb '' \
fastqs/atac_R1_SRR18677642.fastq.gz fastqs/atac_R2_SRR18677642.fastq.gz fastqs/atac_R3_SRR18677642.fastq.gz
```

### `kb-python`

`cellatlas` generates the following `kb-python` code to preprocess the data.

```bash
# build alignment reference
minimap2 -d atac/ref.mmi reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
zcat reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | fold -w 80 > atac/genome.fa
minimap2 -o atac/genome.sam -a -x sr -t 32 atac/ref.mmi fastqs/atac_R1_SRR18677641.fastq.gz fastqs/atac_R3_SRR18677641.fastq.gz
samtools view -@ 8 -o atac/genome.u.bam -b atac/genome.sam
samtools sort -@ 8 -o atac/genome.bam -n -m 8G atac/genome.u.bam
Genrich -t atac/genome.bam -o atac/genome.bed -f atac/genome_peaks.log -v
cat atac/genome.bed | bedtools sort | bedtools merge > atac/peaks.bed
bedtools getfasta -fi atac/genome.fa -bed atac/peaks.bed -fo atac/peaks.fa
cat atac/peaks.fa | awk '{if($1~/>/)print $1\"\t\"$1\"\t\"$1}' > atac/t2g.txt
sed -i 's/>//g' atac/t2g.txt
kallisto index -i atac/index.idx atac/peaks.fa

# perform alignment, error correction, and counting
kb count -i atac/index.idx -g atac/t2g.txt -x 1,8,24:-1,-1,-1:0,0,52,2,0,52 -w ATA-737K-arc-v1.txt -o atac --h5ad -t 2 fastqs/atac_R1_SRR18677642.fastq.gz fastqs/atac_R2_SRR18677642.fastq.gz fastqs/atac_R3_SRR18677642.fastq.gz

mkdir -p atac/counts_mult

# generate read counts instead of UMI counts
bustools count -o atac/counts_mult/cells_x_genes -g atac/t2g.txt -e atac/matrix.ec -t atac/transcripts.txt --genecounts --cm atac/output.unfiltered.bus
```

### `kallisto` and `bustools`

`kb-python` generates the following `kallisto` and `bustools` code to run.

```bash
# Note that the reference is built in a slightly different way using tools like Genrich and minimap2
kallisto index -i atac/index.idx atac/peaks.fa

# alignment with kallisto
kallisto bus -i atac/index.idx -o atac -x 1,8,24:-1,-1,-1:0,0,52,2,0,52 -t 2 fastqs/atac_R1_SRR18677642.fastq.gz fastqs/atac_R2_SRR18677642.fastq.gz fastqs/atac_R3_SRR18677642.fastq.gz

# error correction and counting with bustools
bustools sort -o atac/tmp/output.s.bus -T atac/tmp -t 2 -m 4G atac/output.bus
bustools inspect atac/tmp/output.s.bus
bustools inspect -o atac/inspect.json -w ATA-737K-arc-v1_rc.txt atac/tmp/output.s.bus
bustools correct -o atac/tmp/output.s.c.bus -w ATA-737K-arc-v1_rc.txt atac/tmp/output.s.bus
bustools inspect atac/tmp/output.s.c.bus
bustools sort -o atac/output.unfiltered.bus -T atac/tmp -t 2 -m 4G atac/tmp/output.s.c.bus
bustools inspect atac/output.unfiltered.bus
bustools count -o atac/counts_unfiltered/cells_x_genes -g atac/t2g.txt -e atac/matrix.ec -t atac/transcripts.txt --genecounts atac/output.unfiltered.bus

# this last step is manually added to generate read counts instead of UMI counts
bustools count -o atac/counts_mult/cells_x_genes -g atac/t2g.txt -e atac/matrix.ec -t atac/transcripts.txt --genecounts --cm atac/output.unfiltered.bus
```

Note that in the `kb-python` and `kallisto` `bustools` examples above the `-x` and the `-w` arguments can also be generated directly with `seqspec`:

```bash
$ seqspec index \
-t kb \
-m atac \
-r fastqs/atac_R1_SRR18677642.fastq.gz,fastqs/atac_R2_SRR18677642.fastq.gz,fastqs/atac_R3_SRR18677642.fastq.gz \
spec.yaml
1,8,24:-1,-1,-1:0,0,52,2,0,52

$ seqspec onlist -m atac -r barcode spec.yaml
/path/to/the/relevant/onlist/dogmaseq-dig/ATA-737K-arc-v1_rc.txt
```
