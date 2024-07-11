---
title: Understand the data
date: 2024-07-07
authors:
  - name: A. Sina Booeshaghi
---

# Understand the data

Before we download, QC, and process the data we must first understand the read structure. NB: the term "read" is used in this context to refer to an entry in a FASTQ file (which contains a header, string of nucleotides, spacer, and quality score string.) While sometimes we have to work directly with the BCL files produced by Illumina sequencers, we will assume we are dealing with FASTQ files. In order to perform QC, we must understand the annotation of the nucleotides in the sequencing reads contained in the FASTQ files.

A sequencing library contains molecules with sequencing primers and assay-specific constructs which can contain, for example, barcodes and biological features. A `seqspec` read specification provides an annotation for the sequencing reads and tells us, nominally, where different features are expected in each read. The `seqspec` for our reads can be found here. We will use the `seqspec` command line tool to understand the structure of the data.

Install `seqspec`

```bash
# install seqspec with pip
pip install git+https://github.com/IGVF/seqspec.git

# verify the instal
seqspec --version
```

and then view the specification

```bash
# print out the specification
seqspec print spec.yaml
```

This should print out

TODO: put image here

TODO: Add the picture of the modalities

As is immediately evident from the `seqspec`, we should expect nine FASTQ files, each annotated with modality specific features like barcodes, UMIs and spacer sequences. Note that the number next to the feature on the very right is the length of the feature in number of nucleotides.

We do not recommend acquiring, QC'ing, or analyzing sequencing data without first understanding the structure of the reads- `seqspec` is one way to achieve this.
