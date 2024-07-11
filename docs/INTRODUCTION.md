---
title: Introduction
date: 2024-07-07
authors:
  - name: A. Sina Booeshaghi
---

# You get sequencing data, now what?

This document is intended to be a guide for researchers who wish to generate, analyze, or process genomic sequencing data. It is structure as follows:

- [[#Understand the data]]
- [[#Download the data]]
- [[#Quality control the data]]
- [[#Preprocess the data]]

---

Genomics experiments encompass the design, production, and analysis of sequencing data. Broadly speaking they occur in three steps

1. Design of sequencing readout
2. Experimentation and generation of sequencing library
3. Sequencing and computational analysis

Each step is complex and comprises multiple substeps. Here we discuss the what to do when we get sequencing data in the form of FASTQ files. We first start with understanding the data and metadata, then we transfer the data to our computer, install software for QC, perform QC, preprocess, and perform analysis.

We will make the following three assumptions about the data:

1. Single-cell data was generated with short read Illumina sequencer
2. The data was made available on the machine or on a server (with either FTP or USB access)
3. A `seqspec` was made available for the data

To guide this tutorial, we will use multi-modal (RNA, Protein, Tag, ATAC) data generated with the [DOGMAseq assay](- https://doi.org/10.1186/s13059-022-02698-8) (with DIG permeabilization). A sample dataset derived from the original data can be found here: https://github.com/IGVF/seqspec/tree/main/specs/dogmaseq-dig. This folder contains

- 1 million sequencing reads for each of the four modalities
- A `seqspec` specification that annotates the read structure
- Barcode onlists for the RNA cell barcodes and the ATAC cell barcodes
- Reference protein and tag barcodes for the protein and tag modalities

To download the example data, simply run the following in the command line

```bash
# download the repo
git clone https://github.com/IGVF/seqspec.git

# go into the DOGMAseq-dig folder
cd seqspec/specs/dogmaseq-dig

# list out the contents of the folder
ls -lht
```
