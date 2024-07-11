---
title: Preprocess Matrix
date: 2024-07-07
authors:
  - name: A. Sina Booeshaghi
---

Short description of the steps we are taking (insert image from cell atlas paper)

1. Filter matrix
2. Normalize matrix
3. Assign celltypes or cell categories

# Filter matrix

## Command line

```bash
mx filter -bi barcodes.txt -c 2 2 -bo filter.txt -o filtered_matrix.mtx matrix.mtx
```

## Python

TODO: Verify

```python
import matplotlib.pyplot as plt
from mx.filter import mx_filter
import anndata

adata = anndata.read_h5ad("adata.h5ad")
filtered_matrix = mx_filter(adata.X.copy(), adata.obs.index.values)
```

# Normalize counts

# Assign cell types / categories

## Compiling markers

Finding DEG lists online

Extracting celltype-marker gene pairs

The EC file format

## Assigning cell types/categories
