import json
import os
import sys
from datetime import datetime
from typing import List
from seqspec.seqspec_onlist import run_onlist
from seqspec.seqspec_index import run_index
from seqspec.utils import load_spec, region_ids_in_spec
from seqspec.seqspec_find import run_find_by_type
from cellatlas.UniformData import UniformData

MOD2FEATURE = {
    "TAG": "tags",
    "PROTEIN": "protein",
    "ATAC": "gDNA",
    "RNA": "cDNA",
    "CRISPR": "gRNA",
}


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
        help=("Path(s) to output file"),
        type=str,
        nargs="+",
        default=None,
    )
    # -s is a list of seqspecs
    subparser.add_argument(
        "-s",
        metavar="SEQSPEC",
        help="Path to seqspec file(s)",
        type=str,
        nargs="+",
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
        nargs="+",
    )
    subparser.add_argument(
        "--joint", help="Joint quantification", action="store_true", default=False
    )

    return subparser


def validate_build_args(parser, args):
    fastqs = [f for f in args.fastqs]
    fasta = args.fa
    gtf = args.g
    feature_barcodes = args.fb

    modalities = args.m
    seqspec_fns = args.s
    outputs = args.o
    joint = args.joint

    len(set(fastqs)) == len(fastqs) or parser.error("FASTQs must be unique")

    # Case 1, O,M,S: (1, 1, 1)
    if len(set(modalities)) == 1 and len(outputs) == 1 and len(seqspec_fns) == 1:
        udatas = [
            UniformData(
                seqspec_fns[0],
                modalities[0],
                fastqs,
                fasta,
                gtf,
                feature_barcodes,
                outputs[0],
            )
        ]
        multimodal = False
        joint = False
    # Case 2, O,M,S: (>1, >1, 1)
    elif len(set(modalities)) > 1 and len(outputs) > 1 and len(seqspec_fns) == 1:
        udatas = []
        for o, m in zip(outputs, modalities):
            udatas.append(
                UniformData(
                    seqspec_fns[0],
                    m,
                    fastqs,
                    fasta,
                    gtf,
                    feature_barcodes,
                    o,
                )
            )
        multimodal = True
        joint = False
    # Case 2, O,M,S: (1, >1, >1)
    elif len(set(modalities)) == 1 and len(outputs) > 1 and len(seqspec_fns) > 1:
        udatas = []
        for o, s in zip(outputs, seqspec_fns):
            udatas.append(
                UniformData(
                    s,
                    modalities[0],
                    fastqs,
                    fasta,
                    gtf,
                    feature_barcodes,
                    o,
                )
            )
        multimodal = False
        joint = False
    # Case 2, O,M,S: (1, 1, >1)
    elif len(set(modalities)) == 1 and len(outputs) == 1 and len(seqspec_fns) > 1:
        udatas = []
        for s in seqspec_fns:
            udatas.append(
                UniformData(
                    s,
                    modalities[0],
                    fastqs,
                    fasta,
                    gtf,
                    feature_barcodes,
                    outputs[0],
                )
            )
        multimodal = False
        joint = True
    else:
        raise Exception("Invalid input.")

    # add all of the feature fastqs to all of the udatas
    all_feature_fastqs = []
    for udata in udatas:
        all_feature_fastqs.append(udata.spec_feature_fastqs)
    for udata in udatas:
        udata.all_feature_fastqs = all_feature_fastqs

    for output in outputs:
        if not os.path.isdir(output):
            os.makedirs(output)

    # return run_build(
    #     modalities, fastqs, seqspecs, fasta, gtf, feature_barcodes, outputs
    # )
    return run_build(udatas, multimodal, joint)


def run_build(udatas: List[UniformData], multimodal=False, joint=False):
    # single (len(udata) == 1)
    # multimodal same seqspec
    # multimodal different seqspecs
    # unimodal different seqspecs not joint
    # unimodal different seqspecs joint

    if not joint:
        for udata in udatas:
            run_build_separate(udata)
    else:
        run_build_joint(udatas[0])
    return


def run_build_separate(udata: UniformData):
    call_time = datetime.now().strftime("%a %b %d %H:%M:%S %Y %Z")
    # push the selection of which fastqs to pass into here

    if udata.modality.upper() == "ATAC":
        ref = run_build_ref_joint(
            udata.modality,
            udata.all_feature_fastqs,
            udata.fasta,
            udata.gtf,
            udata.feature_barcodes,
            udata.output,
        )
    else:
        ref = run_build_ref(
            udata.modality,
            udata.spec_feature_fastqs,
            udata.fasta,
            udata.gtf,
            udata.feature_barcodes,
            udata.output,
        )

    count = run_build_count(
        udata.modality,
        udata.spec_all_fastqs,
        udata.x_string,
        udata.onlist_fn,
        udata.output,
    )

    cmds = [
        {"ref": ref},
        {"count": count},
    ]

    run_json = {
        "call": " ".join(sys.argv),
        "start_time": call_time,
        "fastqs": [{"file": f, "source": ""} for f in udata.all_fastqs],
        "seqspec": udata.seqspec_fn,
        "genome_fasta": udata.fasta,
        "genome_gtf": udata.gtf,
        "commands": cmds,
    }

    with open(os.path.join(udata.output, "cellatlas_info.json"), "w") as f:
        print(json.dumps(run_json, indent=4), file=f)
    return


def run_build_joint(udata: UniformData):
    call_time = datetime.now().strftime("%a %b %d %H:%M:%S %Y %Z")

    ref = run_build_ref_joint(
        udata.modality,
        udata.all_feature_fastqs,
        udata.fasta,
        udata.gtf,
        udata.feature_barcodes,
        udata.output,
    )

    # the counting is the same as count_separate (to count_joint), only the fastqs change
    count = run_build_count(
        udata.modality, udata.all_fastqs, udata.x_string, udata.onlist_fn, udata.output
    )

    cmds = [
        {"ref": ref},
        {"count": count},
    ]
    # since the seqspecs have to yield the same technology string
    # we only list the first one
    # though in principle they are not the same
    run_json = {
        "call": " ".join(sys.argv),
        "start_time": call_time,
        "fastqs": [{"file": f, "source": ""} for f in udata.all_fastqs],
        "seqspec": udata.seqspec_fn,
        "genome_fasta": udata.fasta,
        "genome_gtf": udata.gtf,
        "commands": cmds,
    }

    with open(os.path.join(udata.output, "cellatlas_info.json"), "w") as f:
        print(json.dumps(run_json, indent=4), file=f)
    return


def run_build_ref(
    modality: str,
    fastqs: List[str],
    fasta: str,
    gtf: str,
    feature_barcodes: str,
    output: str,
):
    REF = {
        "TAG": build_kb_ref_kite,
        "PROTEIN": build_kb_ref_kite,
        "CRISPR": build_kb_ref_kite,
        "ATAC": build_kb_ref_snATAK,
        "RNA": build_kb_ref_standard,
    }

    ref = REF[modality.upper()](fastqs, fasta, gtf, feature_barcodes, output)

    return ref


def run_build_ref_joint(
    modality: str,
    fastqs: List[List[str]],
    fasta: str,
    gtf: str,
    feature_barcodes: str,
    output: str,
):
    REF = {
        "TAG": build_kb_ref_kite_joint,
        "PROTEIN": build_kb_ref_kite_joint,
        "CRISPR": build_kb_ref_kite_joint,
        "ATAC": build_kb_ref_snATAK_joint,
        "RNA": build_kb_ref_standard_joint,
    }

    # joint ref requires passing all of the feature specific fastqs across all of the udatas
    ref = REF[modality.upper()](
        fastqs,
        fasta,
        gtf,
        feature_barcodes,
        output,
    )

    return ref


def run_build_count(
    modality: str, fastqs: List[str], x_string: str, onlist: str, output: str
):
    COUNT = {
        "TAG": build_kb_count_kite,
        "PROTEIN": build_kb_count_kite,
        "CRISPR": build_kb_count_kite,
        "ATAC": build_kb_count_snATAK,
        "RNA": build_kb_count_standard,
    }

    count = COUNT[modality.upper()](fastqs, x_string, onlist, output)

    return count


def build_kb_ref_standard(
    fastqs: List[str], fasta: str, gtf: str, feature_barcodes: str, output: str
):
    cmd = ["kb ref"]
    cmd.append(f"-i {os.path.join(output, 'index.idx')}")
    cmd.append(f"-g {os.path.join(output, 't2g.txt')}")
    cmd.append(f"-f1 {os.path.join(output, 'transcriptome.fa')}")
    cmd.append(fasta)
    cmd.append(gtf)
    return [" ".join(cmd)]


def build_kb_ref_standard_joint(
    fastqs: List[List[str]], fasta: str, gtf: str, feature_barcodes: str, output: str
):
    cmd = ["kb ref"]
    cmd.append(f"-i {os.path.join(output, 'index.idx')}")
    cmd.append(f"-g {os.path.join(output, 't2g.txt')}")
    cmd.append(f"-f1 {os.path.join(output, 'transcriptome.fa')}")
    cmd.append(fasta)
    cmd.append(gtf)
    return [" ".join(cmd)]


def build_kb_ref_kite(
    fastqs: List[str], fasta: str, gtf: str, feature_barcodes: str, output: str
):
    cmd = ["kb ref --workflow kite"]
    cmd.append(f"-i {os.path.join(output, 'index.idx')}")
    cmd.append(f"-g {os.path.join(output, 't2g.txt')}")
    cmd.append(f"-f1 {os.path.join(output, 'transcriptome.fa')}")
    cmd.append(feature_barcodes)
    return [" ".join(cmd)]


def build_kb_ref_kite_joint(
    fastqs: List[List[str]], fasta: str, gtf: str, feature_barcodes: str, output: str
):
    cmd = ["kb ref --workflow kite"]
    cmd.append(f"-i {os.path.join(output, 'index.idx')}")
    cmd.append(f"-g {os.path.join(output, 't2g.txt')}")
    cmd.append(f"-f1 {os.path.join(output, 'transcriptome.fa')}")
    cmd.append(feature_barcodes)
    return [" ".join(cmd)]


def build_kb_ref_snATAK(
    fastqs: List[str], fasta: str, gtf: str, feature_barcodes: str, output: str
):
    # build minimap ref
    cmd = [f"minimap2 -d {os.path.join(output, 'ref.mmi')} {fasta}"]

    # get sample specific peaks
    cmd += get_peaks(fastqs, fasta, gtf, feature_barcodes, output)

    # merge the peaks
    cmd += [
        f"cat {os.path.join(output,'peaks.*.bed')} | bedtools sort | bedtools merge > {os.path.join(output, 'peaks.bed')}",
    ]

    # Build the peak fasta and index with kallisto
    cmd += [
        f"zcat {fasta} | fold -w 80 > {os.path.join(output, 'genome.fa')}",
        f"bedtools getfasta -fi {os.path.join(output, 'genome.fa')} -bed {os.path.join(output, 'peaks.bed')} -fo {os.path.join(output, 'peaks.fa')}",
        f"cat {os.path.join(output, 'peaks.fa')} | awk '{{if($1~/>/)print $1\"\t\"$1\"\t\"$1}}' > {os.path.join(output, 't2g.txt')}",
        f"sed -i 's/>//g' {os.path.join(output, 't2g.txt')}",
        f"kallisto index -i {os.path.join(output, 'index.idx')} {os.path.join(output, 'peaks.fa')}",
    ]
    return cmd


def build_kb_ref_snATAK_joint(
    fastqs: List[List[str]], fasta: str, gtf: str, feature_barcodes: str, output: str
):
    # build minimap ref
    cmd = [f"minimap2 -d {os.path.join(output, 'ref.mmi')} {fasta}"]

    # fastqs is a list of lists of fastq paths for each sample
    for idx, fqs in enumerate(fastqs):
        cmd += get_peaks(fqs, fasta, gtf, feature_barcodes, output, idx)

    # merge the peaks
    cmd += [
        f"cat {os.path.join(output,'peaks.*.bed')} | bedtools sort | bedtools merge > {os.path.join(output, 'peaks.bed')}",
    ]

    # Build the peak fasta and index with kallisto
    cmd += [
        f"zcat {fasta} | fold -w 80 > {os.path.join(output, 'genome.fa')}",
        f"bedtools getfasta -fi {os.path.join(output, 'genome.fa')} -bed {os.path.join(output, 'peaks.bed')} -fo {os.path.join(output, 'peaks.fa')}",
        f"cat {os.path.join(output, 'peaks.fa')} | awk '{{if($1~/>/)print $1\"\t\"$1\"\t\"$1}}' > {os.path.join(output, 't2g.txt')}",
        f"sed -i 's/>//g' {os.path.join(output, 't2g.txt')}",
        f"kallisto index -i {os.path.join(output, 'index.idx')} {os.path.join(output, 'peaks.fa')}",
    ]

    return cmd


def get_peaks(
    fastqs: List[str],
    fasta: str,
    gtf: str,
    feature_barcodes: str,
    output: str,
    sample_index: int = 0,
):
    fqs = " ".join(fastqs)

    cmd = [
        f"minimap2 -o {os.path.join(output, 'genome.' + str(sample_index) + '.sam')} -a -x sr -t 32 {os.path.join(output, 'ref.mmi')} {fqs}",
        f"samtools view -@ 8 -o {os.path.join(output, 'genome.u.' + str(sample_index) + '.bam')} -b {os.path.join(output, 'genome.' + str(sample_index) + '.sam')}",
        f"samtools sort -@ 8 -o {os.path.join(output,'genome.' + str(sample_index) + '.bam')} -n -m 8G {os.path.join(output, 'genome.u.' + str(sample_index) + '.bam')}",
        f"Genrich -t {os.path.join(output,'genome.' + str(sample_index) + '.bam')} -o {os.path.join(output,'genome.' + str(sample_index) + '.bed')} -f {os.path.join(output,'genome_peaks.' + str(sample_index) + '.log')} -v",
        f"cat {os.path.join(output,'genome.' + str(sample_index) + '.bed')} | bedtools sort | bedtools merge > {os.path.join(output, 'peaks.' + str(sample_index) + '.bed')}",
    ]
    return cmd


def build_kb_count_standard(
    fastqs: List[str], x_string: str, onlist_fn: str, output: str
):
    # make technology string with seqspec
    # get whitelist from seqspec
    cmd = [
        "kb count",
        f"-i {os.path.join(output, 'index.idx')}",
        f"-g {os.path.join(output, 't2g.txt')}",
        f"-x {x_string}",
        f"-w {onlist_fn}",
        f"-o {output}",
        "--h5ad",
        "-t 2",
        " ".join(fastqs),
    ]
    return [" ".join(cmd)]


def build_kb_count_kite(fastqs: List[str], x_string: str, onlist_fn: str, output: str):
    # make technology string with seqspec
    # get whitelist from seqspec
    cmd = ["kb count --workflow kite"]
    cmd.append(f"-i {os.path.join(output, 'index.idx')}")
    cmd.append(f"-g {os.path.join(output, 't2g.txt')}")
    cmd.append(f"-x {x_string}")
    cmd.append(f"-w {onlist_fn}")
    cmd.append(f"-o {output}")
    cmd.append("--h5ad")
    cmd.append("-t 2")
    cmd.append(" ".join(fastqs))
    return [" ".join(cmd)]


def build_kb_count_snATAK(
    fastqs: List[str], x_string: str, onlist_fn: str, output: str
):
    # make technology string with seqspec
    # get whitelist from seqspec
    cmd = build_kb_count_standard(fastqs, x_string, onlist_fn, output)
    cmd.append(f"mkdir -p {os.path.join(output, 'counts_mult')}")
    cmd.append(
        f"bustools count -o {os.path.join(output, 'counts_mult/cells_x_genes')} -g {os.path.join(output, 't2g.txt')} -e {os.path.join(output, 'matrix.ec')} -t {os.path.join(output, 'transcripts.txt')} --genecounts --cm {os.path.join(output, 'output.unfiltered.bus')}"
    )
    return cmd
