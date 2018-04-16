#!/usr/bin/env python3

from enum import Enum
from functools import partial

import numpy as np


class Dred(Enum):
    UMAP = 'UMAP'
    PCA = 'PCA'
    TSNE = 't-SNE'
    def __str__(self):
        return self.value

def get_dred_fnc(dred, opts):
    from umap import UMAP
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE

    def do_umap(x, *args, **kwargs):
        return UMAP(*args, **kwargs).fit_transform(x)
    def do_pca(x, *args, **kwargs):
        return PCA(*args, **kwargs).fit_transform(x)
    def do_tsne(x, *args, **kwargs):
        return TSNE(*args, **kwargs).fit_transform(x)

    if dred == Dred.UMAP:
        return partial(do_umap, n_components=opts.dim, n_neighbors=opts.neighbors)
    if dred == Dred.PCA:
        return partial(do_pca, n_components=opts.dim)
    if dred == Dred.TSNE:
        return partial(do_tsne, n_components=opts.dim, perplexity=opts.perplexity)
    raise ValueError(f"Method {dred:s} not recognized.")


class Cluster(Enum):
    AGGLOMERATIVE = 'Agglomerative'
    KMEANS = 'k-means'
    def __str__(self):
        return self.value

def get_cluster_fnc(cluster, nclusters, opts):
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.cluster import KMeans

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

    def do_kmeans(x, nclusters, *args, **kwargs):
        return np.stack([
            KMeans(*args, n_clusters=n, **kwargs).fit_predict(x)
            for n in nclusters
        ], axis=1) + 1

    if cluster == Cluster.AGGLOMERATIVE:
        return partial(do_agglomerative, nclusters=nclusters)
    if cluster == Cluster.KMEANS:
        return partial(do_kmeans, nclusters=nclusters)
    raise ValueError(f"Method {dred:s} not recognized.")


def get_run_fnc(opts):
    if opts.cluster:
        return get_cluster_fnc(
            opts.cluster_method,
            [int(x) for x in opts.nclusters.split(",")],
            opts)
    return get_dred_fnc(opts.dim_red, opts)


import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--neighbors", metavar="N", type=int, default=256, help="number of neighbors")
parser.add_argument("-d", "--dim", metavar="N", type=int, default=3, help="number of dimensions to project to")
parser.add_argument("paths", metavar="path", nargs="+", type=str, help="input matrix, rows are samples")
parser.add_argument("-o", "--out", metavar="path", type=str, default="", help="prefix for generated output files")
parser.add_argument("-t", "--transpose", action="store_true", help="input matrix, rows are samples")
parser.add_argument("-i", "--individually", action="store_true", help="also run separately for each file")
parser.add_argument("-r", "--relfreq", action="store_true", help="use relative frequency")
parser.add_argument("-p", "--pseudocnt", metavar="FP", type=float, default=0.0, help="pseudo count to add to counts")
parser.add_argument("-f", "--filter", metavar="FP", type=float, default=0.0, help="filter samples with sums less than this value")
parser.add_argument("--dim-red", type=Dred, choices=list(Dred), default=Dred.UMAP, help="method for dimensionality reduction")
parser.add_argument("--perplexity", type=float, default=10, help="preplexity parameter of t-SNE dimensionality reduction")
parser.add_argument("--cluster", action="store_true", help="cluster data (instead of doing dimensionality reduction on it)")
parser.add_argument("--cluster-method", type=Cluster, choices=list(Cluster), default=Cluster.KMEANS, help="clustering method")
parser.add_argument("--nclusters", default="3", help="comma-separated list of number of clusters to use when --cluster is set")

args = parser.parse_args()

run = get_run_fnc(args)

import pandas as pd

matrices = []
colnames = set()
idx = 0
for path in args.paths:
    print("reading path %s" % path)

    matrix = pd.DataFrame(pd.read_csv(path, index_col=0, delimiter='\t')) + args.pseudocnt

    if args.transpose:
        matrix = matrix.transpose()

    print(matrix.shape)

    rs = matrix.sum(axis=1)

    if args.relfreq:
        matrix = (matrix.T / rs).T

    if args.filter > 0.0:
        matrix = matrix.loc[rs >= args.filter, :]
        print(matrix.shape)

    matrices = matrices + [matrix]

    if args.individually:
        print(f"Performing {args.dim_red} individually")
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
print("Union")
print(matrix.shape)
print(f"Performing {args.dim_red} jointly")
embedding = run(matrix)
print(f"Performed {args.dim_red}")
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
