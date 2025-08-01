# Bam2cell


[![Latest PyPI Version][pb]][pypi] [![PyPI Downloads][db]][pypi] [![tests][gb]][yml] 

[gb]: https://github.com/davidrm-bio/bam2cell/actions/workflows/release.yml/badge.svg
[yml]: https://github.com/davidrm-bio/bam2cell/actions/workflows/release.yml
[pb]: https://img.shields.io/pypi/v/bam2cell.svg
[pypi]: https://pypi.org/project/bam2cell/

[db]: https://img.shields.io/pypi/dm/bam2cell?label=pypi%20downloads


A package to split a BAM file based on cell type annotation in 
AnnData objects.


# Usage and Examples

There are two modes **sequential** and **parallel**. The sequential mode will process cell types one by one but is more
disk space friendly, the parallel is more disk space hungry but faster, since it process all cell types at the same time.

A minimal example is shown here:

```python
import bam2cell
import anndata as ad

adata = ad.read_h5ad("data/adata.h5ad")

generator = bam2cell.GenerateCellTypeBAM(adata, 
                                         annot_key="annotation",
                                         output_path="data/",
                                         input_bam="data/AllCellsSorted_toy.bam",
                                         tmp_path="data/",
                                         workers=8,
                                         )
generator.process_all_parallel()  # Case 1 - Process all cell types at the same time
generator.process_cts_sequential() # Case 2 - Process cell types one by one

```

For a more advanced usage, you can use the function `bam2cell`, which allow to process an AnnData with multiple samples.

```python
import bam2cell
import anndata as ad
import pandas as pd

adata = ad.read_h5ad("data/adata.h5ad")
artificial_batch = ["batch1"] * 100 + ["batch2"] * 91
adata.obs["batch"] = pd.Categorical(artificial_batch)
adata.obs["bam_path"] = "data/AllCellsSorted_toy.bam"

bam2cell.bam2cell(adata,
                  annot_key="annotation",
                  input_bam=None,  # Only when we have 1 batch in the AnnData
                  output_path="data/",  
                  tmp_path="data/",
                  bam_key="bam_path",  # For each barcode we have the path to the BAM file
                  batch_key="batch",  
                  mode="parallel",
                  suffix=None,  # Suffix in the barcode to be removed (e.g., BC-1-suffix --> BC-1)
                  prefix=None,  # Prefix in the BC to be removed (e.g., prefix-BC-1 --> BC-1) 
                  workers=8
                  )

```

# Installation

You need to have Python 3.10 or newer installed on your system. There are several alternative options 
to install `bam2cell`:

1. Install the latest release of `bam2cell` from [PyPI](https://pypi.org/project/bam2cell/):
```bash
pip install bam2cell  
```

2. Install the latest development version:
 ```bash
pip install git+https://github.com/davidrm-bio/bam2cell.git@main
```

