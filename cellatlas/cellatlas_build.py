import json
import os
import sys
from datetime import datetime
from seqspec.seqspec_onlist import run_onlist
from seqspec.seqspec_index import run_index
from seqspec.utils import load_spec


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
        "-f",
        metavar="FASTA",
        help=("Path to genome fasta file"),
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
    fastqs = [os.path.abspath(f) for f in args.fastqs]
    seqspec = os.path.abspath(args.s)
    fasta = os.path.abspath(args.f)
    gtf = os.path.abspath(args.g)
    modality = args.m
    output = os.path.abspath(args.o)
    if not os.path.isdir(output):
        os.makedirs(output)
    return run_build(modality, fastqs, seqspec, fasta, gtf, output)


def run_build(modality, fastqs, seqspec_fn, fasta, gtf, output):
    call_time = datetime.now().strftime("%a %b %d %H:%M:%S %Y %Z")

    seqspec = load_spec(seqspec_fn)
    x_string = run_index(
        seqspec, modality, [os.path.basename(i) for i in fastqs], fmt="kb"
    )
    onlist = os.path.abspath(run_onlist(seqspec, modality, "barcode"))

    cmds = [
        {"kb_ref": build_kb_ref(fasta, gtf)},
        {"kb_count": build_kb_count(fastqs, x_string, onlist)},
    ]

    run_json = {
        "call": " ".join(sys.argv),
        "start_time": call_time,
        "fastqs": [{"file": f, "source": ""} for f in fastqs],
        "seqspec": seqspec_fn,
        "genome_fasta": fasta,
        "genome_gtf": gtf,
        "commands": cmds,
    }
    with open(os.path.join(output, "cellatlas_info.json"), "w") as f:
        print(json.dumps(run_json, indent=4), file=f)
    return


def build_kb_ref(genome_fasta, genome_gtf):
    s = f"kb ref -i index.idx -g t2g.txt -f1 transcriptome.fa {genome_fasta} {genome_gtf}"
    return s


def build_kb_count(fastqs, x_string, onlist):
    # make technology string with seqspec
    # get whitelist from seqspec
    s = f"kb count -i index.idx -g t2g.txt -x {x_string} -w {onlist} -o output --h5ad -t 2 {' '.join(fastqs)}"
    return s
