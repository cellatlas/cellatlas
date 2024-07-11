---
title: QC Reads
date: 2024-07-07
authors:
  - name: A. Sina Booeshaghi
---

Now that we have data on our machine we can QC it.

# QC reads with `fastqc`

Data quality control must be assessed so that we are confident in our results. We start by first using the `fastqc` command line tool. This tool will check the reads in FASTQ files against 10 metrics.

NB: The following is summarized from the [`fastq` documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/).

1. Per base sequence quality
   - The Per Base Sequence Quality tool in FastQC visualizes the quality of base calls in FastQ files using a BoxWhisker plot at each position, where the central red line indicates the median quality, the yellow box shows the inter-quartile range (25-75%), the whiskers mark the 10% and 90% points, and a blue line represents the mean; the y-axis divides the quality into green (very good), orange (reasonable), and red (poor) regions, and issues warnings or failures based on lower quartile and median values falling below set thresholds, commonly triggered by long sequencing runs, transient issues, or low coverage in varying-length reads.
2. Per tile sequence quality
   - The Per Tile Sequence Quality tool in FastQC provides a color-coded plot for Illumina libraries retaining original sequence identifiers, showing the deviation of quality scores from the average on a per-tile basis within a flowcell; blue indicates average or above-average quality, while hotter colors signify tiles with below-average quality, with warnings and failures triggered if any tile's mean Phred score deviates by more than 2 or 5 units, respectively, from the base's mean across all tiles, commonly due to transient issues like bubbles or more permanent problems like flowcell smudges.
3. Per sequence quality scores
   - The Per Sequence Quality Scores tool in FastQC identifies sequences with universally low quality values in a run, issuing a warning if the most frequently observed mean quality falls below 27 (0.2% error rate) and a failure if below 20 (1% error rate); this is especially useful for diagnosing systematic issues in part of a sequencing run and the results should be cross-referenced with per-tile qualities when a bi-modal or complex distribution is observed.
4. Per base sequence content
   - The Per Base Sequence Content tool in FastQC plots the proportion of each of the four DNA bases at each position in a file, expecting parallel lines for a balanced genome, and issues a warning or failure if the difference between A and T or G and C exceeds 10% or 20% respectively; while intrinsic biases from libraries primed by random hexamers or fragmented by transposases can't be corrected, the module's warnings may also flag overrepresented sequences, biased fragmentation in the first 12 base pairs, and libraries with inherently skewed base composition, such as sodium bisulphite-treated or aggressively adapter-trimmed libraries.
5. Per sequence GC content
   - The Per Sequence GC Content tool in FastQC measures the GC content across each sequence in a file and compares it to a modeled normal distribution, issuing a warning if deviations from this normal distribution account for more than 15% of the reads and a failure if more than 30%; anomalies in the distribution could indicate library contamination, biased subsets, or systematic biases that won't be flagged if they result in a consistently shifted normal distribution.
6. Per base N content
   - The Per Base N Content tool in FastQC plots the percentage of 'N' base calls at each position in a sequence, triggering a warning if any position has more than 5% 'N' content and a failure if over 20%; higher proportions of 'N' usually indicate poor base calling due to general loss of quality or biased sequence composition, and should be cross-referenced with other quality modules for comprehensive assessment.
7. Sequence Length Distribution
   - The Sequence Length Distribution module in FastQC generates a graph to show the distribution of sequence fragment sizes, issuing a warning if sequences are not of uniform length and an error if any sequence has zero length, with the understanding that variable read lengths may be normal for some sequencing platforms.
8. Sequence Duplication Levels
   - The Duplicate Sequences module in FastQC quantifies sequence duplication levels by analyzing the first 100,000 sequences in a file, truncating reads over 75bp to 50bp for the analysis, and plots the relative number of sequences with different degrees of duplication, issuing a warning if more than 20% of sequences are non-unique and an error if more than 50%, while acknowledging that different library types and scenarios—such as RNA-Seq over-sequencing of highly expressed transcripts—may naturally yield high duplication levels.
9. Overrepresented sequences
   - The Overrepresented Sequences module lists sequences comprising more than 0.1% of the total, analyzing only the first 100,000 sequences to conserve memory, truncates reads over 75bp to 50bp for analysis, and cross-references overrepresented sequences with a database of common contaminants, issuing a warning if any sequence surpasses the 0.1% threshold and an error if exceeding 1%, while noting the possibility of false positives or misses, especially in specialized libraries like small RNA libraries.
10. Adapter Content
    - The Adapter Content module specifically searches for a set of pre-defined Kmers to assess the proportion of your library containing these adapter sequences, generating a cumulative percentage plot that tracks adapter presence at each read position, and issues a warning if any adapter sequence is present in more than 5% of reads and an error if it exceeds 10%, commonly triggering in libraries with short insert sizes relative to read lengths, indicating the need for adapter trimming.

First, we need to install `fastqc`. There is extensive documentation for how to do this on multiple operating systems. Here is the [official documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt) and here is a [simplified version](https://pbertinblog.wordpress.com/fastqc-installation/). If you are on a mac (and use Homebrew), then you can install `fastqc` by simply running the following:

```bash
brew install fastqc
```

Once we have `fastqc` installed, we can evaluate our sequencing data.

```bash
# make a folder for the results
mkdir fqc

# run fastqc on all sequencing reads
fastqc -o fqc/ fastqs/*
```

This will produce two HTML reports in the `fqc` folder that we can open and look at. One HTML file is made for each FASTQ file.

# QC reads with `seqspec`

Now that we have confidence that our sequencing data passes an initial round of QC, we can get more fine grain and QC relevant regions of our reads (e.g. barcodes and UMIs).

Some questions to consider when QCing reads

- For sequences corresponding to cellular/nuclei/sample barcodes
  - What fraction of barcodes are 0-distance from the onlist? 1-dist? 2?
  - For the barcodes that are on the onlist, are there overrepresented sequences?
    - Do we expect uniformity?
    - Are the barcodes in the expected position within the read?
- For sequences corresponding to UMIs
  - Are there overrepresented sequences? Is their distribution more or less uniform?
  - Can we compute the expected number of UMIs on the bead that was captured?
- For the sequences corresponding to biological features
  - Are there overrepresented sequences?
- Generally, are linkers and spacers where we expect them to be?

## QC Barcodes

Often single-cells/nuclei are labeled with nucleotide barcodes that arise from an "onlist", or a list of permissible barcodes. We want to extract a representative sample of them from our reads to QC them. To do this, we will use `seqspec` to get the position of the barcodes, and `seqkit` to extract the sequences, and python scripts to perform some simple analysis.

To do this, we will need to install `seqkit`, an all-purpose command line tool for manipulating FASTQ data. Here is official documentation on how to do this on the command line. If you are on a mac (and use HomeBrew), then seqkit can simply be installed by running

```bash
brew install seqkit
```

Extract the first 1,000 barcodes (indexed by `seqspec`) with `seqkit` and look at their distribution:

```bash
zcat fastqs/tag_R1_SRR18677640.fastq.gz | seqkit subseq -r $(seqspec index -t seqkit -m tag -s barcode -r tag_R1_SRR18677640.fastq.gz spec.yaml)  | awk '(NR-2)%4==0' | head -1000 | sort | uniq -c | sort -nr | head
```

## QC UMIs

Extract the first 1,000 UMIs (indexed by `seqspec`) with `seqkit` and check the uniformity of the distribution. UMIs are designed to be random so we expect a uniform distribution. In other words, we expect maximum entropy. We compute the fraction of the maximum entropy that the distribution of the first 1,000 UMIs give us (closer to 1 indicates more uniformity, closer to 0 indicates less uniformity).

```bash
# todo
```

Next, we extract the tag from the FASTQs and look at its distribution. The process is the same as before.

```bash
zcat fastqs/tag_R2_SRR18677640.fastq.gz | seqkit subseq -r $(seqspec index -t seqkit -m tag -s tag -r tag_R2_SRR18677640.fastq.gz spec.yaml)  | awk '(NR-2)%4==0' | head -1000 | sort | uniq -c | sort -nr | head # look at the first 1000 barcodes
```
