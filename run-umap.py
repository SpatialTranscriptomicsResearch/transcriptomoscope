#!/usr/bin/env python3

from enum import Enum
from functools import partial
import logging
from logging import INFO

import numpy as np


logging.basicConfig(level=INFO)
LOG = logging.getLogger().log


class Dred(Enum):
    """
    Enum of available dimensionality reduction methods.
    """
    UMAP = 'UMAP'
    PCA = 'PCA'
    TSNE = 't-SNE'
    def __str__(self):
        return self.value

def get_dred_fnc(dred, opts):
    """
    Returns the dimensionality reduction function corresponding to the given
    `dred` method.

    Parameters
    ----------
    dred : Dred
    opts : dict
        Dictionary with options for the given dimensionality reduction method.

    Returns
    -------
    dred_fnc: array-like -> array-like
    """
    from umap import UMAP
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE

    def do_umap(x, *args, **kwargs):
        return UMAP(*args, **kwargs).fit_transform(x)
    def do_pca(x, *args, **kwargs):
        return PCA(*args, **kwargs).fit_transform(x)
    def do_tsne(x, perplexity, initial_dims, *args, **kwargs):
        if x.shape[1] > initial_dims:
            LOG(INFO, f'reducing initial dimensionality to {initial_dims}')
            pca_kwargs = kwargs
            pca_kwargs.pop("n_components")
            y = PCA(*args, **pca_kwargs, n_components=initial_dims).fit_transform(x)
        else:
            y = x

        local_perplexity = min(y.shape[0] / 3.5, perplexity)
        if not local_perplexity == perplexity:
            LOG(INFO, f'lowering perplexity to {local_perplexity}')
        return TSNE(*args, perplexity=local_perplexity, **kwargs).fit_transform(y)

    if dred == Dred.UMAP:
        return partial(do_umap, n_components=opts.dim, n_neighbors=opts.neighbors)
    if dred == Dred.PCA:
        return partial(do_pca, n_components=opts.dim)
    if dred == Dred.TSNE:
        return partial(do_tsne, n_components=opts.dim, perplexity=opts.perplexity, initial_dims=opts.initial_dims)
    raise ValueError(f"Method {dred:s} not recognized.")


class Cluster(Enum):
    """
    Enum of available clustering methods.
    """
    AGGLOMERATIVE = 'Agglomerative'
    GMIX = 'GaussianMixture'
    KMEANS = 'k-means'
    def __str__(self):
        return self.value

def get_cluster_fnc(cluster, nclusters, opts):
    """
    Returns the clustering function corresponding to the given `cluster`
    method.

    Parameters
    ----------
    cluster : Cluster
    opts : dict
        Dictionary with options for the given clustering method.

    Returns
    -------
    cluster_fnc: array-like -> array-like
    """
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.cluster import KMeans
    from sklearn.mixture import GaussianMixture

    def cache(f):
        from tempfile import TemporaryDirectory
        def wrapper(*args, **kwargs):
            with TemporaryDirectory() as td:
                return f(*args, cache_dir=td, **kwargs)
        return wrapper

    @cache
    def do_agglomerative(x, nclusters, cache_dir, *args, **kwargs):
        clusterer = AgglomerativeClustering(memory=cache_dir)
        def do_clustering(n):
            clusterer.set_params(n_clusters=n)
            return clusterer.fit_predict(x)
        return np.stack(map(do_clustering, nclusters), axis=1) + 1

    def do_gmix(x, nclusters, *args, **kwargs):
        clusterer = GaussianMixture(*args, **kwargs)
        return np.stack([
            clusterer.set_params(n_components=n).fit(x).predict(x)
            for n in nclusters
        ], axis=1) + 1

    def do_kmeans(x, nclusters, *args, **kwargs):
        return np.stack([
            KMeans(*args, n_clusters=n, **kwargs).fit_predict(x)
            for n in nclusters
        ], axis=1) + 1

    if cluster == Cluster.AGGLOMERATIVE:
        return partial(do_agglomerative, nclusters=nclusters)
    if cluster == Cluster.GMIX:
        return partial(do_gmix, nclusters=nclusters)
    if cluster == Cluster.KMEANS:
        return partial(do_kmeans, nclusters=nclusters)
    raise ValueError(f"Method {cluster:s} not recognized.")

def hyphen_range(s):
    """ Takes a range in form of "a-b" and generate a list of numbers between a and b inclusive.
    Also accepts comma separated ranges like "a-b,c-d,f" will build a list which will include
    Numbers from a to b, a to d and f"""
    s = "".join(s.split()) #removes white space
    r = set()
    for x in s.split(','):
        t = x.split('-')
        if len(t) not in [1,2]: raise SyntaxError("hash_range is given its arguement as "+s+" which seems not correctly formated.")
        r.add(int(t[0])) if len(t)==1 else r.update(set(range(int(t[0]),int(t[1])+1)))
    l = list(r)
    l.sort()
    return l

def get_run_fnc(opts):
    """Returns the transformation that will be applied to the data"""
    if opts.cluster is not None:
        LOG(INFO, f'Using clustering method {opts.cluster:s}')
        return get_cluster_fnc(
            opts.cluster,
            hyphen_range(opts.nclusters),
            opts)
    LOG(INFO, f'Using dimensionality reduction method {opts.dim_red:s}')
    return get_dred_fnc(opts.dim_red, opts)


import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--neighbors", metavar="N", type=int, default=256, help="number of neighbors for UMAP")
parser.add_argument("-d", "--dim", metavar="N", type=int, default=3, help="number of dimensions to project to")
parser.add_argument("paths", metavar="path", nargs="+", type=str, help="input matrix, rows are samples")
parser.add_argument("-o", "--out", metavar="path", type=str, default="", help="prefix for generated output files")
parser.add_argument("-t", "--transpose", action="store_true", help="input matrix, rows are samples")
parser.add_argument("-i", "--individually", action="store_true", help="also run separately for each file")
parser.add_argument("-a", "--absfreq", action="store_true", help="use absolute frequency; default is to use relative frequency")
parser.add_argument("-p", "--pseudocnt", metavar="FP", type=float, default=0.0, help="pseudo count to add to counts")
parser.add_argument("-f", "--filter", metavar="FP", type=float, default=0.0, help="filter samples with sums less than this value")
parser.add_argument("--dim-red", type=Dred, choices=list(Dred), default=Dred.UMAP, help="method for dimensionality reduction")
parser.add_argument("--perplexity", type=float, default=50, help="perplexity parameter of t-SNE dimensionality reduction")
parser.add_argument("--initial_dims", type=int, default=50, help="option for t-SNE: the number of dimensions that should be retained in the initial PCA step (default: 50)")
parser.add_argument("--cluster", type=Cluster, choices=list(Cluster), help="clustering method. if set, clusters the data instead of doing dimensionality reduction on it.")
parser.add_argument("--nclusters", default="1-12", type=str, help="comma-separated list of number ranges of clusters to use when --cluster is set, e.g. \"2,5-7,12\". [default = 1-12]")

args = parser.parse_args()

run = get_run_fnc(args)

import pandas as pd

matrices = []
colnames = set()
idx = 0
for path in args.paths:
    LOG(INFO, f'Reading path {path:s}')

    matrix = pd.DataFrame(pd.read_csv(path, index_col=0, delimiter='\t')) + args.pseudocnt

    if args.transpose:
        matrix = matrix.transpose()

    rs = matrix.sum(axis=1)

    if not args.absfreq:
        matrix = (matrix.T / rs).T

    if args.filter > 0.0:
        matrix = matrix.loc[rs >= args.filter, :]

    matrices = matrices + [matrix]

    if args.individually:
        LOG(INFO, 'Running individual analysis '
            f'(samples={matrix.shape[0]:d}, features={matrix.shape[1]:d})')
        embedding = run(matrix)
        df = pd.DataFrame(embedding, index=matrix.index)
        df.to_csv("%sindividual_%03d.txt" % (args.out, idx), sep="\t")
        idx = idx + 1

    colnames = colnames.union(matrix.columns.values)

def union_matrix(matrix, colnames):
    m = pd.DataFrame(0.0, index=matrix.index, columns=colnames)
    cols = colnames.intersection(matrix.columns.values)
    cols = list(cols)
    m.loc[:,cols] = matrix.loc[:,cols]
    return m

matrices = list(map(lambda matrix: union_matrix(matrix, colnames), matrices))

matrix = pd.concat(matrices)
LOG(INFO, 'Performing joint analysis '
    f'(samples={matrix.shape[0]:d}, features={matrix.shape[1]:d})')
embedding = run(matrix)

pd.DataFrame(embedding).to_csv("%sjoint.txt" % args.out, sep="\t")

idx = 0
begin = 0
for matrix in matrices:
    nrow, ncol = matrix.shape
    end = begin + nrow

    df = pd.DataFrame(embedding[range(begin, end), :], index=matrix.index)
    df.to_csv("%sjoint_%03d.txt" % (args.out, idx), sep="\t")

    idx = idx + 1
    begin = end
