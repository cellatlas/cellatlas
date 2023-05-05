from typing import List
from seqspec.seqspec_index import run_index, format_kallisto_bus
from seqspec.Assay import Assay


def create_kb_string_from_spec(spec: Assay, modality: str, regions: List[str]):
    indices = []
    for r in regions:
        index = run_index(spec, modality, r)
        indices.append({r: index})
    x = format_kallisto_bus(indices)
    return x
