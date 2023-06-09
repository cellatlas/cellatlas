import json
import os
import sys
from datetime import datetime
from seqspec.seqspec_onlist import run_onlist
from seqspec.seqspec_index import run_index
from seqspec.utils import load_spec
from seqspec.seqspec_find import run_find_by_type


def setup_build_args(parser):
    subparser = parser.add_parser(
        "build",
        description="build cellatlas for one observation",
        help="build cellatlas for one observation",
    )
    subparser.add_argument("fastqs", nargs="+", help="List of FASTQ files")
    subparser.add_argument(
        "-o",
        metavar="OUT",
        help=("Path to output file"),
        type=str,
        default=None,
    )
    subparser.add_argument(
        "-s",
        metavar="SEQSPEC",
        help=("Path to seqspec file"),
        type=str,
        default=None,
    )
    subparser.add_argument(
        "-fa",
        metavar="FASTA",
        help=("Path to genome fasta file"),
        type=str,
        default=None,
    )
    subparser.add_argument(
        "-fb",
        metavar="Feature barcodes",
        help=("Path to feature barcode file"),
        type=str,
        default=None,
    )
    subparser.add_argument(
        "-g",
        metavar="GTF",
        help=("Path to genome gtf file"),
        type=str,
        default=None,
    )
    subparser.add_argument(
        "-m",
        metavar="Modality",
        help=("Modality of the data"),
        type=str,
        default=None,
    )
    return subparser


def validate_build_args(parser, args):
    fastqs = [f for f in args.fastqs]
    seqspec = args.s
    fasta = args.fa
    gtf = args.g
    feature_barcodes = args.fb
    modality = args.m
    output = args.o

    if not os.path.isdir(output):
        os.makedirs(output)

    return run_build(modality, fastqs, seqspec, fasta, gtf, feature_barcodes, output)


def run_build(modality, fastqs, seqspec, fasta, gtf, feature_barcodes, output):
    call_time = datetime.now().strftime("%a %b %d %H:%M:%S %Y %Z")

    ref = run_build_ref(modality, fastqs, seqspec, fasta, gtf, feature_barcodes, output)
    count = run_build_count(modality, fastqs, seqspec, output)

    cmds = [
        {"ref": ref},
        {"count": count},
    ]

    run_json = {
        "call": " ".join(sys.argv),
        "start_time": call_time,
        "fastqs": [{"file": f, "source": ""} for f in fastqs],
        "seqspec": seqspec,
        "genome_fasta": fasta,
        "genome_gtf": gtf,
        "commands": cmds,
    }
    with open(os.path.join(output, "cellatlas_info.json"), "w") as f:
        print(json.dumps(run_json, indent=4), file=f)
    return


def run_build_ref(modality, fastqs, seqspec_fn, fasta, gtf, feature_barcodes, output):
    spec = load_spec(seqspec_fn)
    MOD2FEATURE = {
        "TAG": "tags",
        "PROTEIN": "tags",
        "ATAC": "gDNA",
        "RNA": "cDNA",
        "CRISPR": "gRNA",
    }

    # search the modality, and the feature type associated with the modality to get the fastq file names from the region_id
    rgns = run_find_by_type(spec, modality, MOD2FEATURE.get(modality.upper(), ""))
    relevant_fqs = [rgn.parent_id for rgn in rgns]
    # get the paths from fastqs that match relevant_fqs
    fqs = [f for f in fastqs if os.path.basename(f) in relevant_fqs]
    # add spec to ref build, find the fastq files that are gDNA then select those from fastqs parameter
    REF = {
        "TAG": build_kb_ref_kite,
        "PROTEIN": build_kb_ref_kite,
        "CRISPR": build_kb_ref_kite,
        "ATAC": build_kb_ref_snATAK,
        "RNA": build_kb_ref_standard,
    }
    ref = REF[modality.upper()](fqs, fasta, gtf, feature_barcodes, output)

    return ref


def run_build_count(modality, fastqs, seqspec_fn, output):
    seqspec = load_spec(seqspec_fn)
    x_string = run_index(
        seqspec, modality, [os.path.basename(i) for i in fastqs], fmt="kb"
    )
    onlist = run_onlist(seqspec, modality, "barcode")
    # get onlist path relative to seqspec_fn path
    onlist = os.path.join(os.path.dirname(seqspec_fn), onlist)

    COUNT = {
        "TAG": build_kb_count_kite,
        "PROTEIN": build_kb_count_kite,
        "CRISPR": build_kb_count_kite,
        "ATAC": build_kb_count_snATAK,
        "RNA": build_kb_count_standard,
    }

    count = COUNT[modality.upper()](fastqs, x_string, onlist, output)

    return count


def build_kb_ref_standard(fastqs, fasta, gtf, feature_barcodes, output):
    cmd = ["kb ref"]
    cmd.append(f"-i {os.path.join(output, 'index.idx')}")
    cmd.append(f"-g {os.path.join(output, 't2g.txt')}")
    cmd.append(f"-f1 {os.path.join(output, 'transcriptome.fa')}")
    cmd.append(fasta)
    cmd.append(gtf)
    return [" ".join(cmd)]


def build_kb_ref_kite(fastqs, fasta, gtf, feature_barcodes, output):
    cmd = ["kb ref --workflow kite"]
    cmd.append(f"-i {os.path.join(output, 'index.idx')}")
    cmd.append(f"-g {os.path.join(output, 't2g.txt')}")
    cmd.append(f"-f1 {os.path.join(output, 'transcriptome.fa')}")
    cmd.append(feature_barcodes)
    return [" ".join(cmd)]


def build_kb_ref_snATAK(fastqs, fasta, gtf, feature_barcodes, output):
    fqs = " ".join(fastqs)

    cmd = [
        f"minimap2 -d {os.path.join(output, 'ref.mmi')} {fasta}",
        f"zcat {fasta} | fold -w 80 > {os.path.join(output, 'genome.fa')}",
        f"minimap2 -o {os.path.join(output, 'genome.sam')} -a -x sr -t 32 {os.path.join(output, 'ref.mmi')} {fqs}",
        f"samtools view -@ 8 -o {os.path.join(output, 'genome.u.bam')} -b {os.path.join(output, 'genome.sam')}",
        f"samtools sort -@ 8 -o {os.path.join(output,'genome.bam')} -n -m 8G {os.path.join(output, 'genome.u.bam')}",
        f"Genrich -t {os.path.join(output,'genome.bam')} -o {os.path.join(output,'genome.bed')} -f {os.path.join(output,'genome_peaks.log')} -v",
        f"cat {os.path.join(output,'genome.bed')} | bedtools sort | bedtools merge > {os.path.join(output, 'peaks.bed')}",
        f"bedtools getfasta -fi {os.path.join(output, 'genome.fa')} -bed {os.path.join(output, 'peaks.bed')} -fo {os.path.join(output, 'peaks.fa')}",
        f"cat {os.path.join(output, 'peaks.fa')} | awk '{{if($1~/>/)print $1\"\t\"$1\"\t\"$1}}' > {os.path.join(output, 't2g.txt')}",
        f"sed -i 's/>//g' {os.path.join(output, 't2g.txt')}",
        f"kallisto index -i {os.path.join(output, 'index.idx')} {os.path.join(output, 'peaks.fa')}",
    ]
    # cmd = ["kb ref --workflow kite"]
    # cmd.append(f"-i {os.path.join(output, 'index.idx')}")
    # cmd.append(f"-g {os.path.join(output, 't2g.txt')}")
    # cmd.append(f"-f1 {os.path.join(output, 'transcriptome.fa')}")
    # cmd.append(feature_barcodes)
    return cmd


def build_kb_count_standard(fastqs, x_string, onlist, output):
    # make technology string with seqspec
    # get whitelist from seqspec
    cmd = [
        "kb count",
        f"-i {os.path.join(output, 'index.idx')}",
        f"-g {os.path.join(output, 't2g.txt')}",
        f"-x {x_string}",
        f"-w {onlist}",
        f"-o {output}",
        "--h5ad",
        "-t 2",
        " ".join(fastqs),
    ]
    return [" ".join(cmd)]


def build_kb_count_kite(fastqs, x_string, onlist, output):
    # make technology string with seqspec
    # get whitelist from seqspec
    cmd = ["kb count --workflow kite"]
    cmd.append(f"-i {os.path.join(output, 'index.idx')}")
    cmd.append(f"-g {os.path.join(output, 't2g.txt')}")
    cmd.append(f"-x {x_string}")
    cmd.append(f"-w {onlist}")
    cmd.append(f"-o {output}")
    cmd.append("--h5ad")
    cmd.append("-t 2")
    cmd.append(" ".join(fastqs))
    return [" ".join(cmd)]


def build_kb_count_snATAK(fastqs, x_string, onlist, output):
    # make technology string with seqspec
    # get whitelist from seqspec
    cmd = build_kb_count_standard(fastqs, x_string, onlist, output)
    cmd.append(f"mkdir -p {os.path.join(output, 'counts_mult')}")
    cmd.append(
        f"bustools count -o {os.path.join(output, 'counts_mult/cells_x_genes')} -g {os.path.join(output, 't2g.txt')} -e {os.path.join(output, 'matrix.ec')} -t {os.path.join(output, 'transcripts.txt')} --genecounts --cm {os.path.join(output, 'output.unfiltered.bus')}"
    )
    return cmd
