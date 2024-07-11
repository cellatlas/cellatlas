---
title: Download the data
date: 2024-07-07
authors:
  - name: A. Sina Booeshaghi
---

# Download the data

As illustrated by the `seqspec`, there are many FASTQ files which (if we did not already have them) would need to be downloaded. Often FASTQ files located on an remote FTP server, remote server, Dropbox, or local sequencing machine. In either case, it is straightforward to download the data.

## Remote FTP server

If FASTQ files are stored on a remote FTP server, you will need to have access credentials (username and password). This will allow you to log on to the FTP server, navigate to the FASTQ files, and download them to your computer.

```bash
# example taken from https://ena-docs.readthedocs.io/en/latest/retrieval/file-download.html
ftp ftp.sra.ebi.ac.uk
Name: anonymous # provided username
Password:       # provided password
ftp> cd vol1/fastq/ERR164/ERR164407
ftp> get ERR164407.fastq.gz
```

The `get` command is used to specify the file to download to your machine.

Oftentimes direct links to to sequencing data are provided and we can use the command line tools `wget` or `curl` to download them without needing a username and password.

```bash
curl -Ls -O R1.fastq.gz ftp://link.to.fastq/R1.fastq.gz
wget
```

## Remote server

If FASTQ files are stored on a remote server that requires credentials to log on, you will need a user name and password to download the data to your own computer. With an account, there are many ways to download the data- one such was is with the command line tool `scp`.

```bash
scp username@remoteserver /path/to/fastq/files/*fastq.gz
```

## AWS s3 bucket

If FASTQ files are stored in an AWS s3 bucket, you will need to [install the AWS command line tool](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) so that the files can be transferred to your machine. Additionally, you may need access credentials to view and download data in s3 buckets. Information on how to get access credentials can be found [here](https://docs.aws.amazon.com/powershell/latest/userguide/pstools-appendix-sign-up.html)

```bash
aws configure # supply the Access Key and Secret
aws s3 cp s3://path/to/files/ . --recursive --exclude "*" --include "*.fastq.gz
```

Info on using glob ("\*") with the cp command: https://stackoverflow.com/questions/38834708/how-can-i-use-wildcards-to-cp-a-group-of-files-with-the-aws-cli

## Dropbox

If FASTQ files are stored on dropbox, we will need to obtain the direct file download link and use `wget` or `curl` to download them. Here is [documentation](https://www.dropboxforum.com/t5/Dropbox-API-Support-Feedback/Downloading-Dropbox-files-using-curl-or-wget/td-p/485857) from Dropbox discussing how to do this. A direct download link can usually be obtained by clicking "share" next to the file. Click "Share File" and finally "Copy link". `wget` can be used with that link

```bash
wget https://www.dropbox.com/s/ew2jket9lisdf4oor/example.fastq.gz
```

## Local sequencing machine

Illumina machines often directly transfer FASTQ files to a connected computer. Obtaining FASTQ files from this local machine requires connected a USB storage device and manually transferring the files over. Illumina [documentation](https://support-docs.illumina.com/SHARE/NetworkSecurity/Content/SHARE/NetworkSecurity/DataOutputStorage.htm) explains when transfers can be initiated after sequencing.

## SRA

todo: prefetch and fastq-dump
