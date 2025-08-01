# Set-Up
import os
import shutil
import uuid
import multiprocessing
from pathlib import Path
from tqdm import tqdm
from typing import Literal

import anndata as ad
import pysam
import pysam.samtools

from ._utils import timing, logger, clean_bcs

DEFAULT_TMP_PATH = "/tmp"


class GenerateCellTypeBAM:
    """Allow to split a BAM based on cell types.

    This class allow to split a BAM files for a sample into separate files based on the cell type
    annotation. If an AnnData contains several samples, it should be subsetted to be run per sample to
    avoid misalignment due to shared barcodes.

    ..note::
        The barcodes in the AnnData should match the barcodes in the BAM files.

    Parameters
    ----------
    adata
        Processed AnnData with cell type annotation.
    annot_key
        Name of `obs` column with cell type annotation. Any categorical column is accepted.
    output_path
        Path to save the BAM files.
    input_bam
        Path to the input BAM file for the sample.
    tmp_path
        Path to save tmp files. Default is '/tmp'.
    workers
        Number of threads to use during the splitting process.

    Examples
    --------
    >>> import bam2cell
    >>> import anndata
    >>> adata = anndata.read_h5ad('data/adata.h5ad')
    >>> generator = bam2cell.GenerateCellTypeBAM(adata, "annotation", "data/", "data/AllCellsSorted_toy.bam", workers=20)
    >>> generator.process_all_parallel()  # Option 1 - Can lead to high disk usage
    >>> generator.process_cts_sequential()  # Option 2 - Disk usage friendly

    """

    def __init__(self,
                 adata: ad.AnnData,
                 annot_key: str,
                 output_path: str | Path,
                 input_bam: str | Path,
                 tmp_path: str | Path = "/tmp",
                 workers: int = 20,
                 ) -> None:

        self.adata = adata
        self.annot_key = annot_key
        self.categories = list(adata.obs[annot_key].unique())
        self.output_path = self._convert_path(output_path)
        self.input_bam = self._convert_path(input_bam)
        self.workers = workers
        self.tmp_path = self._convert_path(tmp_path)
        return

    def _create_tmp_folder(self):
        tmpdir_path = Path(self.tmp_path) / f"BAM_{uuid.uuid4().hex}"
        tmpdir_path.mkdir(parents=True, exist_ok=False)

        for ct in self.categories:
            new_folder = tmpdir_path / ct
            new_folder.mkdir(parents=True, exist_ok=False)
        self.tmp_path = tmpdir_path
        return None

    def _get_bcs_for_ct(self, ct: str, idx: int) -> dict:
        return {bc: idx for bc in self.adata.obs.index[self.adata.obs[self.annot_key] == ct]}

    @staticmethod
    def _convert_path(path: str | Path) -> Path:
        if isinstance(path, str):
            return Path(path)
        else:
            return path

    @staticmethod
    def _index_bam(bam_path: str | Path) -> None:
        pysam.index(bam_path)
        return None

    @staticmethod
    def process_contig(args) -> None:
        in_bam, out_bam, contig, db, workers = args

        with pysam.AlignmentFile(in_bam, "rb", threads=workers) as in_bam:
            with pysam.AlignmentFile(out_bam, "wb", header=in_bam.header, threads=workers) as out_bam:
                for read in in_bam.fetch(contig):
                    if read.has_tag("CB") and read.get_tag("CB") in db:
                        out_bam.write(read)
        return None

    def _merge_sort_and_index_bam(self, ct: str) -> None:
        files = [os.path.join(self.tmp_path, ct, f) for f in os.listdir(self.tmp_path / ct) if
                 '.bam' in f and '.bai' not in f]

        pysam.merge(
            "-f",
            "-@",
            str(self.workers),
            os.path.join(self.output_path, f"{ct}.bam"),
            *files
        )

        pysam.sort(
            "-@",
            str(self.workers),
            "-o",
            os.path.join(self.output_path, f"{ct}_sorted.bam"),
            os.path.join(self.output_path, f"{ct}.bam"),
        )

        self._index_bam(os.path.join(self.output_path, f"{ct}_sorted.bam"))
        os.remove(os.path.join(self.output_path, f"{ct}.bam"))
        shutil.rmtree(self.tmp_path / ct)

    @timing
    def process_all_parallel(self) -> None:
        """Split the BAM file parallelizing everything.

        :return:
        """

        self._create_tmp_folder()

        with pysam.AlignmentFile(str(self.input_bam), "rb") as bam:
            contigs = bam.references

        args = []
        for idx, ct in enumerate(self.categories):
            db_bcs = self._get_bcs_for_ct(ct, idx)
            for contig in contigs:
                out_bam = self.tmp_path / ct / f"{contig}.bam"
                args.append((self.input_bam, out_bam, contig, db_bcs, self.workers))

        with multiprocessing.Pool(self.workers) as pool:
            for _ in tqdm(pool.imap_unordered(self.process_contig, args), total=len(args), desc="Contigs processed"):
                pass

        logger.info('Merging and indexing Bam files')
        with multiprocessing.Pool(self.workers) as pool:
            for _ in tqdm(pool.imap_unordered(self._merge_sort_and_index_bam, self.categories),
                          total=len(self.categories), desc="CellTypes processed"):
                pass
        return None

    @timing
    def process_cts_sequential(self) -> None:
        """Split the BAM file sequentially, parallelizing for each cell type the process.

        :return:
        """

        self._create_tmp_folder()

        with pysam.AlignmentFile(str(self.input_bam), "rb") as bam:
            contigs = bam.references

        for idx, ct in tqdm(enumerate(self.categories), desc='CellType processed', total=len(self.categories)):
            db_bcs = self._get_bcs_for_ct(ct, idx)
            args = []
            for contig in contigs:
                out_bam = self.tmp_path / ct / f"{contig}.bam"
                args.append((self.input_bam, out_bam, contig, db_bcs, self.workers))

            with multiprocessing.Pool(self.workers) as pool:
                for _ in pool.imap_unordered(self.process_contig, args):
                    pass
            self._merge_sort_and_index_bam(ct)

        return None


def bam2cell(
        adata: ad.AnnData,
        annot_key: str,
        output_path: str | Path,
        input_bam: str | Path = None,
        bam_key: str = None,
        batch_key: str = None,
        tmp_path: str | Path = None,
        mode: Literal["sequential", "parallel"] = "sequential",
        suffix: str = None,
        prefix: str = None,
        workers: int = 5,
) -> None:
    """Split BAM files based on AnnData annotation.

    There are two modes: `sequential` and `parallel`. The sequential mode will process each cell type one by one,
    parallelizing the analysis and is more disk space friendly, however, it takes longer. The parallel mode will
    process all the cell types at the same time, even though it can be much faster it can lead to high disk space usage
    due to the tmp files generated.

    :param adata: AnnData object.
    :param annot_key: Key in `obs` with cell type annotation. Any other type of categorical column can also be used.
    :param output_path: Path to save the BAM files for each cell type.
    :param input_bam: Path to the input BAM file. Specify only if the AnnData contains only one sample.
    :param bam_key: Key in `obs` with the path to the BAM file for each barcode.
    :param batch_key: Key in `obs` with sample information.
    :param tmp_path: Tmp folder to save intermediate files. If set to None, the files will be saved in '/tmp'
    :param mode: Mode to run the splitting. If set to sequential, each cell type or category with be processed sequentially
                 parallelizing the process. This mode is more disk space friendly. If set to parallel, all the cell type will
                 be processed simultaneously but is more disk space hungry.
    :param suffix: Suffix in the barcodes that need to be removed.
    :param prefix: Prefix in the barcodes that need to be removed.
    :param workers: Number of threads to use for the parallelization.
    :return: Returns None. If batch_key is set, a subfolder will be created for each batch in `output_path`. For each
             cell type or group a BAM file will be saved in the `output_folder`.


    See Also
    --------
        :func:`bam2cell.GenerateCellTypeBAM`: class to generate the BAM files for each cell type.


    Examples
    --------
    >>> import bam2cell
    >>> import anndata
    >>> import pandas as pd
    >>> adata = ad.read_h5ad("data/adata.h5ad")
    >>> artificial_batch = ["batch1"] * 100 + ["batch2"] * 91
    >>> adata.obs["batch"] = pd.Categorical(artificial_batch)
    >>> adata.obs["bam_path"] = "data/AllCellsSorted_toy.bam"
    >>> bam2cell.bam2cell(adata,
    ...                  annot_key="annotation",
    ...                  output_path="data/",
    ...                  tmp_path="data/",
    ...                  bam_key="bam_path",
    ...                  batch_key="batch",
    ...                  mode="parallel",
    ...                  workers=8
    ...                  )

    """

    tmp_path = DEFAULT_TMP_PATH if tmp_path is None else tmp_path

    adata_copy = clean_bcs(adata, suffix=suffix, prefix=prefix)

    # Case 1 - Path to BAM files are in the object
    if bam_key is not None:
        assert batch_key is not None, "Specify the obs key with batch information"
        assert input_bam is None, "Specify either the bam_key or input_bam not both"

        for batch in adata_copy.obs[batch_key].unique():
            logger.info(f"Generating BAM files for {batch}")
            sdata = adata_copy[adata_copy.obs[batch_key] == batch].copy()
            input_bam = sdata.obs[bam_key].to_list()[0]
            batch_path = Path(output_path) / batch

            os.makedirs(batch_path, exist_ok=True)  # Create a subfolder for each sample

            generator = GenerateCellTypeBAM(sdata,
                                            annot_key=annot_key,
                                            output_path=batch_path,
                                            input_bam=input_bam,
                                            tmp_path=tmp_path,
                                            workers=workers)
            if mode == "sequential":
                generator.process_cts_sequential()
            elif mode == "parallel":
                generator.process_all_parallel()
            else:
                raise Exception(f"{mode} is not a valid key, use 'sequential' or 'parallel'")

    # Case 2 - Path to BAM files is provided
    elif input_bam is not None:

        os.makedirs(output_path, exist_ok=True)  # Make sure the output path exist

        generator = GenerateCellTypeBAM(adata_copy,
                                        annot_key=annot_key,
                                        output_path=output_path,
                                        input_bam=input_bam,
                                        tmp_path=tmp_path,
                                        workers=workers)
        if mode == "sequential":
            generator.process_cts_sequential()
        elif mode == "parallel":
            generator.process_all_parallel()
        else:
            raise Exception(f"{mode} is not a valid key, use 'sequential' or 'parallel'")

   # Case 3 - Neither of both are specified or both are specified.
    else:
        raise Exception("Specify either input_bam or bam_key")
