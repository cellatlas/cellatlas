---
title: QC Matrix
date: 2024-07-07
authors:
  - name: A. Sina Booeshaghi
---

Preprocessing produces count matrices that can also be QC'ed. Below are important count matrix level QC measures that can be useful for assessing the quality of the data.

We will use the [`mx`](https://github.com/cellatlas/mx/) tool to generate these matrix level QC metrics. First install `mx`:

```bash
pip install git+https://github.com/cellatlas/mx/
```

# Matrix-level metrics

By way of example, here are whole matrix level measures that can be used to QC.

- **ncells**: The number of barcodes quantified in the matrix
  - Too low and there was an issue with barcode diversity or error-correction to an incorrect onlist
- **ngenes**: The number of genes quantified in the matrix
  - This will often be equal to the number of genes used to construct the alignment index
- **nvals**: The number of entries in the matrix that are non-zero
  - Often a measure for the diversity of barcode-umi pairs in the library
- **density**: `nvals/(ncells*ngenes)`
  - Can be thought of as an efficiency measure of umi capture
- **avg_per_cell**: The average number of counts per cell
- **avg_per_gene**: The average number of counts per gene
- **min_cell**: The number of counts for the cell with the fewest counts
- **max_cell**: The number of counts for the cell with the most counts
- **total_count**: The sum of the counts in the count matrix
- **overdispersion**: `sigma^2 = mu + alpha*mu`, computed on the genes

## Command-line

Then navigate to where you have your matrix files (note the matrix orientation here is `cells x genes`) and run the `mx inspect` command:

```bash
mx inspect -o inspect.json -gi cells_x_genes.genes.txt -a all -bi cells_x_genes.barcodes.txt cells_x_genes.mtx
```

Your json file will look like something like this

```bash
$ cat inspect.json
{
    ncells: 1000,
    ngenes: 20000,
    nvals: 50000,
    density: 0.0025,
    avg_per_cell: 3.4,
    avg_per_gene: 5.6,
    min_cell: 0,
    max_cell: 100,
    total_count: 123123,
    overdispersion: 0.123
}
```

## Python

If you have your count matrix loaded as an anndata you can generate a these metrics from within a python environment

```python
import anndata
from mx.mx_inspect import mx_inspect_rows, mx_inspect_cols, mx_inspect

adata = anndata.read_h5ad("adata.h5ad")

metrics = mx_inspect(adata.X.copy())
print(metrics)
```

# Cell and gene-level metrics

We can also compute cell-wise and gene-wise metrics that are useful in more fine-grained QC metrics. For all cells or all genes we compute

- **counts_min**: the minimum number of counts across
- **counts_max**: maximum number of counts
- **counts_sum**: sum of counts
- **counts_mean**: `counts_sum / size of axis`
- **counts_nnzero**: number of non-zero entries <= size of axis
- **counts_variance**: variance of counts across axis

## Command line

To generate these, run the `mx inspect` command with the `-a rows` for cell metrics or `-a cols` for gene metadata

```bash
# cells
mx inspect -o inspect.cells.txt -gi cells_x_genes.genes.txt -a rows -bi cells_x_genes.barcodes.txt cells_x_genes.mtx

# genes
mx inspect -o inspect.genes.txt -gi cells_x_genes.genes.txt -a cols -bi cells_x_genes.barcodes.txt cells_x_genes.mtx
```

## Python

If you have your count matrix loaded as an anndata you can generate a these metrics from within a python environment

```python
import anndata
from mx.mx_inspect import mx_inspect_rows, mx_inspect_cols, mx_inspect

adata = anndata.read_h5ad("adata.h5ad")

adata.obs = mx_inspect_rows(adata.X.copy(), adata.obs.index.values).copy()
adata.var = mx_inspect_cols(adata.X.copy(), adata.var.index.values).copy()
```

# Visualizing metrics

These measures make visualizing QC measures straightforward

- knee plot on cells and genes
- features detected vs counts
- genes mean vs variance

Plots can be generated from within python

## Command line

The `mx` tool has the facility to generate the knee plot directly from the command line

```bash
mx plot -o knee.png -m knee -gi cells_x_genes.genes.txt -a cols -bi cells_x_genes.barcodes.txt cells_x_genes.mtx
```

## Python

You can also make the knee plot from within a python environment

```python
import matplotlib.pyplot as plt
from mx.mx_plot import mx_knee_plot
import anndata

adata = anndata.read_h5ad("adata.h5ad")

fig, ax = plt.subplots(figsize=(10,10))
ax = mx_plot_knee(mtx, ax)

fig.savefig(knee.png, dpi=300)
fig.show()
```
