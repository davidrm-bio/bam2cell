import sys
import time
import logging
import warnings
import functools
import anndata as ad

warnings.filterwarnings("ignore")

logger = logging.getLogger("bam2cell")


def _setup_logger() -> None:
    """Logger settings.

    :return:
    """
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter("%(asctime)s - %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.propagate = False
    logger.setLevel(logging.INFO)
    return None


def info(msg: str):
    """Produce an info message.

    :param msg:
    :return:
    """
    logger.info(msg)


_setup_logger()


def timing(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        logger = logging.getLogger(func.__module__)
        start_time = time.time()
        result = func(*args, **kwargs)
        duration = (time.time() - start_time) / 60
        metric = 'min'
        if duration > 60:  # More than 60 minutes --> show as hours
            duration = duration/60
            metric = 'h'
        logger.info(f"{func.__qualname__} completed in {duration:.2f} {metric}")
        return result

    return wrapper


def clean_bcs(
        adata: ad.AnnData,
        suffix: str = None,
        prefix: str = None
) -> ad.AnnData:
    """Remove prefix or suffix from the barcodes of an AnnData object.

    :param adata: AnnData object.
    :param suffix: Suffix in the barcodes
    :param prefix: Prefix in the barcodes
    :return: Returns AnnData object with the suffix or prefix removed.

    Example
    -------

    >>> import anndata
    >>> adata = anndata.read_h5ad('data/adata.h5ad')
    >>> # Case 1 - Suffix present in the barcodes
    >>> adata.obs_names = [bc + '-batch1' for bc in adata.obs_names]  # Add a suffix
    >>> adata.obs_names[:2]
    Index(['AAACGAACAATCGTCA-1-batch1', 'AAAGGGCAGATTGCGG-1-batch1'], dtype='object')
    >>> adata = clean_bcs(adata, suffix='-batch1')
    >>> adata.obs_names[:2]
    Index(['AAACGAACAATCGTCA-1', 'AAAGGGCAGATTGCGG-1'], dtype='object')
    >>> # Case 2 - Suffix present in the barcodes
    >>> adata.obs_names = ['batch1-' + bc  for bc in adata.obs_names]  # Add a prefix
    >>> adata.obs_names[:2]
    Index(['batch1-AAACGAACAATCGTCA-1', 'batch1-AAAGGGCAGATTGCGG-1'], dtype='object')
    >>> adata = clean_bcs(adata, prefix='batch1-')
    >>> adata.obs_names[:2]
    Index(['AAACGAACAATCGTCA-1', 'AAAGGGCAGATTGCGG-1'], dtype='object')

    """
    adata_copy = adata.copy()
    if suffix is not None:
        adata_copy.obs_names = adata_copy.obs_names.str.split(suffix).str[0]
    elif prefix is not None:
        adata_copy.obs_names = adata_copy.obs_names.str.split(prefix).str[-1]
    elif suffix is not None and prefix is not None:
        adata_copy.obs_names = adata_copy.obs_names.str.split(suffix).str[0].str.split(prefix).str[-1]
    else:
        pass
    return adata_copy