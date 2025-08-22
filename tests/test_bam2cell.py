import os
import shutil
from time import time
import anndata as ad
import pandas as pd

import bam2cell

# test_bam2cell_multiple uses parallel
# def test_bam2cell_parallel():
#     adata = ad.read_h5ad("data/adata.h5ad")
#
#     start_time = time()
#     generator = bam2cell.GenerateCellTypeBAM(adata, "annotation", "data/", "data/AllCellsSorted_toy.bam", tmp_path='data/', workers=10)
#     generator.process_all_parallel()
#     end_time = time()
#     print ("Elapse time: ", round((end_time - start_time) / 60, 2), ' min')  # ~6 min
#     files = os.listdir("data/")
#     files_tmp = os.listdir("data/" + [f for f in files if 'BAM_' in f][0])
#     files_tmp = [f for f in files_tmp if not f.startswith('.')]
#     expected_files = ['T_cells_sorted.bam', 'T_cells_sorted.bam.bai',
#                       'Monocytes_sorted.bam', 'Monocytes_sorted.bam.bai',
#                       'pDC_sorted.bam',  'pDC_sorted.bam.bai',
#                       'NK_sorted.bam.bai', 'NK_sorted.bam']
#
#     assert  all(file in files for file in expected_files), "All expected files not found"
#     assert len(files_tmp) == 0, "tmp folder not empty"
#
#     # Clean-Up
#     shutil.rmtree("data/" + [f for f in files if 'BAM_' in f][0])
#     for f in expected_files:
#         os.remove("data/" + f)
#     return


def test_bam2cell_sequential():
    adata = ad.read_h5ad("data/adata.h5ad")

    start_time = time()
    generator = bam2cell.GenerateCellTypeBAM(adata, "annotation", "data/", "data/AllCellsSorted_toy.bam",
                                             tmp_path='data/', workers=10)
    generator.process_cts_sequential()
    end_time = time()
    print("Elapse time: ", round((end_time - start_time) / 60, 2), ' min') #  29.2  min
    files = os.listdir("data/")
    files_tmp = os.listdir("data/" + [f for f in files if 'BAM_' in f][0])
    files_tmp = [f for f in files_tmp if not f.startswith('.')]
    expected_files = ['T_cells_sorted.bam', 'T_cells_sorted.bam.bai',
                      'Monocytes_sorted.bam', 'Monocytes_sorted.bam.bai',
                      'pDC_sorted.bam', 'pDC_sorted.bam.bai',
                      'NK_sorted.bam.bai', 'NK_sorted.bam']

    assert all(file in files for file in expected_files), "All expected files not found"
    assert len(files_tmp) == 0, "tmp folder not empty"

    # Clean-Up
    shutil.rmtree("data/" + [f for f in files if 'BAM_' in f][0])
    for f in expected_files:
        os.remove("data/" + f)
    return


def test_bam2cell_multiple():
    adata = ad.read_h5ad("data/adata.h5ad")
    artificial_batch = ["batch1"] * 100 + ["batch2"] * 91
    adata.obs["batch"] = pd.Categorical(artificial_batch)
    adata.obs["bam_path"] = "data/AllCellsSorted_toy.bam"

    start_time = time()
    bam2cell.bam2cell(adata,
                      annot_key="annotation",
                      output_path="data/",
                      tmp_path="data/",
                      bam_key="bam_path",
                      batch_key="batch",
                      mode="parallel",
                      workers=8
                      )
    end_time = time()
    print("Elapse time: ", round((end_time - start_time) / 60, 2), ' min') # ~ 12.93 min

    expected_files = ['T_cells_sorted.bam', 'T_cells_sorted.bam.bai',
                      'Monocytes_sorted.bam', 'Monocytes_sorted.bam.bai',
                      'pDC_sorted.bam', 'pDC_sorted.bam.bai',
                      'NK_sorted.bam.bai', 'NK_sorted.bam']

    # Check that a folder for each batch has been created
    files = os.listdir("data/")
    for batch in adata.obs["batch"].unique():
        assert batch in files, 'batch folder not in output path'
        inner_folder = os.listdir("data/" + batch)
        assert all(file in inner_folder for file in expected_files)

    # Check that the tmp folders are empty
    files_tmp = []
    for f in [f for f in files if 'BAM_' in f]:
        files_tmp = files_tmp + os.listdir("data/" + f)
    files_tmp = [f for f in files_tmp if not f.startswith('.')]
    assert len(files_tmp) == 0

    # Clean-Up
    for f in [f for f in files if 'BAM_' in f]:
        shutil.rmtree("data/" + f)
    for f in adata.obs["batch"].unique():
        shutil.rmtree("data/" + f)
    return

