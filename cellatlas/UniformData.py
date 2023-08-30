from seqspec.utils import load_spec
from seqspec.seqspec_find import run_find_by_type
from seqspec.utils import region_ids_in_spec
from seqspec.seqspec_index import run_index
from seqspec.seqspec_onlist import run_onlist
import os
from typing import List

MOD2FEATURE = {
    "TAG": "tags",
    "PROTEIN": "protein",
    "ATAC": "gDNA",
    "RNA": "cDNA",
    "CRISPR": "gRNA",
}


class UniformData:
    def __init__(
        self,
        seqspec_fn: str,
        modality: str,
        fastqs: List[str],
        fasta: str,
        gtf: str,
        feature_barcodes: str,
        output: str,
        all_feature_fastqs: List[List[str]] = [[]],
        x_string: str = "",
        onlist_fn: str = "",
    ) -> None:
        self.seqspec_fn = seqspec_fn
        self.seqspec = load_spec(seqspec_fn)
        self.modality = modality
        self.output = output
        self.all_fastqs = fastqs
        self.all_feature_fastqs = all_feature_fastqs
        rids_in_spec = region_ids_in_spec(
            self.seqspec, self.modality, [os.path.basename(i) for i in self.all_fastqs]
        )
        self.spec_all_fastqs = [
            f for f in self.all_fastqs if os.path.basename(f) in rids_in_spec
        ]

        # filter the fastqs to feature fastqs (the type being feature associated with the modality)
        rgns = run_find_by_type(
            self.seqspec, self.modality, MOD2FEATURE.get(self.modality.upper(), "")
        )
        relevant_fqs = [rgn.parent_id for rgn in rgns]
        fqs = [f for f in self.all_fastqs if os.path.basename(f) in relevant_fqs]
        self.spec_feature_fastqs = fqs

        # note the use of rids_in_spec here, which is the same as the rids_in_spec above
        self.x_string = run_index(self.seqspec, self.modality, rids_in_spec, fmt="kb")

        self.fasta = fasta
        self.gtf = gtf
        self.feature_barcodes = feature_barcodes

        onlist = run_onlist(self.seqspec, self.modality, "barcode")
        # get onlist path relative to seqspec_fn path
        self.onlist_fn = os.path.join(os.path.dirname(self.seqspec_fn), onlist)
