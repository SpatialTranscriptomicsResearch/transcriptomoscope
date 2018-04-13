#!/usr/bin/env Rscript
library("deldir")
library("optparse")

option_list = list(
  make_option(c("-c", "--columns"), type="numeric", default=NULL,
              help="number of columns in grid plots [default = ceiling(sqrt(n)), where n is the number of samples]",
              metavar="N"),
  make_option(c("-b", "--blackbg"), type="logical", default=FALSE,
              help="create plots with a black background", action="store_true"),
  make_option(c("-t", "--transpose"), type="logical", default=FALSE,
              help="ensure that genes are rows and columns are spots",
              action="store_true"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="specify output directory",
              action="store"),
  make_option(c("-inv", "--invertgrayscale"), type="logical", default=FALSE,
              help="invert grayscale in output plots (if number of columns is more than 3)",
              action="store_true"),
  make_option(c("-sd", "--selectdims"), type="character", default=NULL,
              help="subselect three columns to compare",
              action="store"),
  make_option(c("-B", "--border"), type="logical", default=FALSE,
              help="draw borders around tesselated objects",
              action="store_true")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser, positional_arguments = c(1, Inf));

#-------------------------------------------------------------------------------

ncols = opt$options$columns
black.bg = opt$options$blackbg
transpose = opt$options$transpose
ncols = opt$options$columns
outdir = opt$options$outdir
bwinv = opt$options$invertgrayscale
sdims = opt$options$selectdims
draw.border = opt$options$border

if (!is.null(outdir)) {
  outdir <- paste0(outdir, "vtess.pdf")
} else {
  outdir <- getwd()
}
if (!is.null(sdims)) {
  ind <- as.integer(strsplit(sdims, split = "x")[[1]])
  stopifnot(length(ind) == 3)
}

paths = opt$args

st.load.matrix <- function(path)
  as.matrix(read.table(file = path, sep="\t", header = T, row.names = 1))

# Load files
d <- list()
for(path in paths) {
  d[[path]] = st.load.matrix(path)
  if (transpose) {
    d[[path]] <- t(d[[path]])
  }
  if (!is.null(sdims)) {
    d[[path]] <- d[[path]][, ind]
  }
}

# Set dimensions of plot
if (ncol(d[[1]]) > 3) {
  n = ncol(d[[1]])
} else {
  n = length(d)
}
if (!is.null(ncols)) {
  nc = ncols
  nr = ceiling(n/ncols)
} else {
  nc = ceiling(sqrt(n))
  nr = ceiling(n/nc)
}

mins = apply(sapply(d, function(x) apply(x, 2, min)), 1, min)
maxs = apply(sapply(d, function(x) apply(x, 2, max)), 1, max)

parse.coords <- function(n, delimiter=coord.delimiter)
  apply(do.call(rbind, strsplit(n, split = "x")), 2, as.numeric)

pdf(file = outdir, width = 6*nc, height = 6*nr)
par(mfrow = c(nr, nc), mar = c(0,0,0,0), bg = ifelse(black.bg, "black", "white"))
for(path in paths) {
  print(path)
  coords <- parse.coords(rownames(d[[path]]))
  x <- coords[, 1]
  y <- coords[, 2]
  y <- 36 - y
  vtess <- deldir(x, y)
  tiles <- tile.list(vtess)
  z = d[[path]]
  ranges = maxs - mins
  z = t((t(z) - mins) / ranges)
  cpi <- chull(x, y)
  cp <- list(x = x[cpi], y = y[cpi])
  
  # Plot each column in Grayscale if there are more than 3 columns
  if (ncol(z) > 3) {
    par(mfrow = c(nr, nc), mar = c(0, 0, 0, 0), bg = ifelse(black.bg, "black", "white"))
    for (i in 1:ncol(z)) {
      v <- z[, i]
      plot(x, y, type = 'n', bty = 'none', axes = FALSE, xlim = c(0, 34), ylim = c(0, 36))
      if (bwinv) {
        v <- 1 - v
      }
      plot(tiles, fillcol = grey(v), showpoints = FALSE, border = ifelse(draw.border, TRUE, NA), clipp = cp, add = TRUE)
    }
  } else {
    plot(x, y, type = 'n', bty = 'none', axes = FALSE, xlim = c(0, 34), ylim = c(0, 36))
    plot(tiles, fillcol = rgb(z), showpoints = FALSE, border = ifelse(draw.border, TRUE, NA), clipp = cp, add = TRUE)
  }
}
dev.off()
