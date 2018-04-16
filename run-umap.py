#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--neighbors", metavar="N", type=int, default=256, help="number of neighbors")
parser.add_argument("-d", "--dim", metavar="N", type=int, default=3, help="number of dimensions to project to")
parser.add_argument("paths", metavar="path", nargs="+", type=str, help="input matrix, rows are samples")
parser.add_argument("-t", "--transpose", action="store_true", help="input matrix, rows are samples")
parser.add_argument("-i", "--individually", action="store_true", help="also perform UMAP separately for each file")
parser.add_argument("-r", "--relfreq", action="store_true", help="use relative frequency")
parser.add_argument("-p", "--pseudocnt", metavar="FP", type=float, default=0.0, help="pseudo count to add to counts")
parser.add_argument("-f", "--filter", metavar="FP", type=float, default=0.0, help="filter samples with sums less than this value")

args = parser.parse_args()

import umap
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
        print("Performing UMAP individually")
        embedding = umap.UMAP(n_components=args.dim, n_neighbors=args.neighbors).fit_transform(matrix)
        df = pd.DataFrame(embedding, index=matrix.index)
        df.to_csv("umap_individual_%03d.txt" % idx, sep="\t")
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
print("Performing UMAP jointly")
embedding = umap.UMAP(n_components=args.dim, n_neighbors=args.neighbors).fit_transform(matrix)
print("Performed UMAP")
pd.DataFrame(embedding).to_csv("umap_joint.txt", sep="\t")

idx = 0
begin = 0
for matrix in matrices:
    nrow, ncol = matrix.shape
    end = begin + nrow

    df = pd.DataFrame(embedding[range(begin, end), :], index=matrix.index)
    df.to_csv("umap_joint_%03d.txt" % idx, sep="\t")

    idx = idx + 1
    begin = end
