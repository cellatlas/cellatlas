import json
import os
import sys
from datetime import datetime
from typing import List
from seqspec.seqspec_onlist import run_onlist
from seqspec.seqspec_index import run_index
from seqspec.utils import load_spec, region_ids_in_spec
from seqspec.seqspec_find import run_find_by_type

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
    # subparser.add_argument(
    #     "--joint", help="Joint quantification", action="store_true", default=False
    # )

    return subparser


def validate_build_args(parser, args):
    fastqs = [f for f in args.fastqs]
    fasta = args.fa
    gtf = args.g
    feature_barcodes = args.fb

    modalities = args.m
    seqspecs = args.s
    outputs = args.o

    run_build = run_build_independent

    len(set(fastqs)) == len(fastqs) or parser.error("FASTQs must be unique")

    # Case 1, O,M,S: (1, 1, 1)
    if len(set(modalities)) == 1 and len(outputs) == 1 and len(seqspecs) == 1:
        run_build = run_build_independent
    # Case 2, O,M,S: (>1, >1, 1)
    elif len(set(modalities)) > 1 and len(outputs) > 1 and len(seqspecs) == 1:
        run_build = run_build_independent
    # Case 2, O,M,S: (1, >1, >1)
    elif len(set(modalities)) == 1 and len(outputs) > 1 and len(seqspecs) > 1:
        run_build = run_build_independent
    # Case 2, O,M,S: (1, 1, >1)
    elif len(set(modalities)) == 1 and len(outputs) == 1 and len(seqspecs) > 1:
        run_build = run_build_joint
    else:
        raise Exception("Invalid input.")

    for output in outputs:
        if not os.path.isdir(output):
            os.makedirs(output)

    return run_build(
        modalities, fastqs, seqspecs, fasta, gtf, feature_barcodes, outputs
    )


def run_build_independent(
    modalities: List[str],
    fastqs: List[str],
    seqspecs: List[str],
    fasta: str,
    gtf: str,
    feature_barcodes: str,
    outputs: List[str],
):
    if len(set(modalities)) == len(outputs) == len(seqspecs) == 1:
        return run_build_single(
            modalities[0], fastqs, seqspecs[0], fasta, gtf, feature_barcodes, outputs[0]
        )
    # Case 2, O,M,S: (>1, >1, 1)
    elif len(set(modalities)) == len(outputs) > 1 and len(seqspecs) == 1:
        return run_build_independent_multiple_mm(
            modalities, fastqs, seqspecs[0], fasta, gtf, feature_barcodes, outputs
        )
    # Case 2, O,M,S: (1, >1, >1)
    elif len(set(modalities)) == 1 and len(outputs) == len(seqspecs) > 1:
        return run_build_independent_multiple_sm(
            modalities[0], fastqs, seqspecs, fasta, gtf, feature_barcodes, outputs
        )
    return


def run_build_independent_multiple_mm(
    modalities: List[str],
    fastqs: List[str],
    seqspec: str,
    fasta: str,
    gtf: str,
    feature_barcodes: str,
    outputs: List[str],
):
    spec = load_spec(seqspec)
    for m, o in zip(modalities, outputs):
        rids = [os.path.basename(f) for f in fastqs]
        found = region_ids_in_spec(spec, m, rids)
        run_build_single(m, found, seqspec, fasta, gtf, feature_barcodes, o)
    return


def run_build_independent_multiple_sm(
    modality: str,
    fastqs: List[str],
    seqspecs: List[str],
    fasta: str,
    gtf: str,
    feature_barcodes: str,
    outputs: List[str],
):
    for spec, o in zip(seqspecs, outputs):
        found = region_ids_in_spec(
            spec, modality, [os.path.basename(f) for f in fastqs]
        )
        run_build_single(modality, found, spec, fasta, gtf, feature_barcodes, o)
    return


# this will only ever run if the modalities are the same and joint is set, is automatically run_build_single
def run_build_joint(
    modalities: List[str],
    fastqs: List[str],
    seqspecs: List[str],
    fasta: str,
    gtf: str,
    feature_barcodes: str,
    outputs: List[str],
):
    # only one modality
    # only one output
    # must have multiple seqspecs and the technology string must be the same for both ie bc umi must be in the same location for all

    run_build_joint_single(
        modalities[0], fastqs, seqspecs, fasta, gtf, feature_barcodes, outputs[0]
    )
    return


def run_build_joint_single(
    modality: str,
    fastqs: List[str],
    seqspecs: List[str],
    fasta: str,
    gtf: str,
    feature_barcodes: str,
    output: str,
):
    call_time = datetime.now().strftime("%a %b %d %H:%M:%S %Y %Z")

    ref = run_build_ref_joint(
        modality, fastqs, seqspecs, fasta, gtf, feature_barcodes, output
    )
    count = run_build_count_joint(modality, fastqs, seqspecs, output)

    cmds = [
        {"ref": ref},
        {"count": count},
    ]

    run_json = {
        "call": " ".join(sys.argv),
        "start_time": call_time,
        "fastqs": [{"file": f, "source": ""} for f in fastqs],
        "seqspec": ",".join(seqspecs),
        "genome_fasta": fasta,
        "genome_gtf": gtf,
        "commands": cmds,
    }

    with open(os.path.join(output, "cellatlas_info.json"), "w") as f:
        print(json.dumps(run_json, indent=4), file=f)
    return


def run_build_single(
    modality: str,
    fastqs: List[str],
    seqspec: str,
    fasta: str,
    gtf: str,
    feature_barcodes: str,
    output: str,
):
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


def run_build_ref_joint(
    modality: str,
    fastqs: List[str],
    seqspecs: List[str],
    fasta: str,
    gtf: str,
    feature_barcodes: str,
    output: str,
):
    specs = [load_spec(ss) for ss in seqspecs]

    REF = {
        "TAG": build_kb_ref_kite_joint,
        "PROTEIN": build_kb_ref_kite_joint,
        "CRISPR": build_kb_ref_kite_joint,
        "ATAC": build_kb_ref_snATAK_joint,
        "RNA": build_kb_ref_standard_joint,
    }

    fqs = []
    for spec in specs:
        ss_rgns = run_find_by_type(
            spec, modality, MOD2FEATURE.get(modality.upper(), "")
        )
        # We want all of the fastq files that are relevant to the modality and that are provided
        relevant_fqs = [rgn.parent_id for rgn in ss_rgns]
        # get the paths from fastqs that match relevant_fqs
        fqs += [f for f in fastqs if os.path.basename(f) in relevant_fqs]

    ref = REF[modality.upper()](fqs, fasta, gtf, feature_barcodes, output)

    #
    return ref


def run_build_ref(modality, fastqs, seqspec_fn, fasta, gtf, feature_barcodes, output):
    spec = load_spec(seqspec_fn)

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


# this needs to be fixed
def run_build_count_joint(modality, fastqs, seqspec_fns, output):
    # For now assume that the technology strings match
    # In the future TODO check that the technology strings match
    # For now assume that the onlist files match
    # In the future TODO check that the onlist files match

    x_strings = []
    joined_fastqs = []
    for seqspec_fn in seqspec_fns:
        seqspec = load_spec(seqspec_fn)
        # search the modality, and the feature type associated with the modality to get the fastq file names from the region_id
        rgns = run_find_by_type(
            seqspec, modality, MOD2FEATURE.get(modality.upper(), "")
        )
        relevant_fqs = [rgn.parent_id for rgn in rgns]
        # get the paths from fastqs that match relevant_fqs
        fqs = [f for f in fastqs if os.path.basename(f) in relevant_fqs]
        print(fqs)
        print(modality)
        x_strings += run_index(seqspec, modality, fqs, fmt="kb")
        print(x_strings)
        joined_fastqs += fqs

    x_string = x_strings[0]  # assumes the technology strings are the same

    onlist = run_onlist(
        seqspec, modality, "barcode"
    )  # assumes that the onlists are the same
    # get onlist path relative to seqspec_fn path
    onlist = os.path.join(os.path.dirname(seqspec), onlist)

    COUNT = {
        "TAG": build_kb_count_kite,
        "PROTEIN": build_kb_count_kite,
        "CRISPR": build_kb_count_kite,
        "ATAC": build_kb_count_snATAK,
        "RNA": build_kb_count_standard,
    }

    count = COUNT[modality.upper()](joined_fastqs, x_string, onlist, output)

    return count


def build_kb_ref_standard(fastqs, fasta, gtf, feature_barcodes, output):
    cmd = ["kb ref"]
    cmd.append(f"-i {os.path.join(output, 'index.idx')}")
    cmd.append(f"-g {os.path.join(output, 't2g.txt')}")
    cmd.append(f"-f1 {os.path.join(output, 'transcriptome.fa')}")
    cmd.append(fasta)
    cmd.append(gtf)
    return [" ".join(cmd)]


def build_kb_ref_standard_joint(fastqs, fasta, gtf, feature_barcodes, output):
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


def build_kb_ref_kite_joint(fastqs, asta, gtf, feature_barcodes, output):
    cmd = ["kb ref --workflow kite"]
    cmd.append(f"-i {os.path.join(output, 'index.idx')}")
    cmd.append(f"-g {os.path.join(output, 't2g.txt')}")
    cmd.append(f"-f1 {os.path.join(output, 'transcriptome.fa')}")
    cmd.append(feature_barcodes)
    return [" ".join(cmd)]


def build_kb_ref_snATAK(fastqs, fasta, gtf, feature_barcodes, output):
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


def build_kb_ref_snATAK_joint(fastqs, fasta, gtf, feature_barcodes, output):
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


def get_peaks(fastqs, fasta, gtf, feature_barcodes, output, sample_index=0):
    fqs = " ".join(fastqs)

    cmd = [
        f"minimap2 -o {os.path.join(output, 'genome.' + str(sample_index) + '.sam')} -a -x sr -t 32 {os.path.join(output, 'ref.mmi')} {fqs}",
        f"samtools view -@ 8 -o {os.path.join(output, 'genome.u.' + str(sample_index) + '.bam')} -b {os.path.join(output, 'genome.' + str(sample_index) + '.sam')}",
        f"samtools sort -@ 8 -o {os.path.join(output,'genome.' + str(sample_index) + '.bam')} -n -m 8G {os.path.join(output, 'genome.u.' + str(sample_index) + '.bam')}",
        f"Genrich -t {os.path.join(output,'genome.' + str(sample_index) + '.bam')} -o {os.path.join(output,'genome.' + str(sample_index) + '.bed')} -f {os.path.join(output,'genome_peaks.' + str(sample_index) + '.log')} -v",
        f"cat {os.path.join(output,'genome.' + str(sample_index) + '.bed')} | bedtools sort | bedtools merge > {os.path.join(output, 'peaks.' + str(sample_index) + '.bed')}",
    ]
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
