#!/usr/bin/env Rscript
library("deldir")
library("optparse")
library("alphahull")
library("RColorBrewer")

# TODO write out a PDF with a color gradient scale bar

ARRAY.COLS = 31
ARRAY.ROWS = 33

add.grid = FALSE

use.convex.hull = FALSE

option_list = list(
  make_option(c("-1", "--one"), type="logical", default=FALSE, action="store_true",
              help="create one plot by interpreting three columns as red, green, and blue values. use this to plot dimensionality reduction results"),
  make_option(c("-O", "--order"), type="logical", default=TRUE, action="store_false",
              help="change plotting order: pages correspond to samples, and panels to columns; by default pages correspond to columns and panels to samples."),
  make_option(c("-c", "--columns"), type="numeric", default=NULL,
              help="number of columns in grid plots [default = ceiling(sqrt(n)), where n is the number of samples]",
              metavar="N"),
  make_option(c("-m", "--margin"), type="numeric", default=0,
              help="margin size. use 0 to include the frame, and use -1 to exclude the frame [default = 0]",
              metavar="N"),
  make_option(c("-s", "--split"), type="logical", default=FALSE,
              help="put each plot on an individual page", action="store_true"),
  make_option(c("", "--points"), type="logical", default=FALSE,
              help="plot points for the spots", action="store_true"),
  make_option(c("", "--outline"), type="logical", default=FALSE,
              help="draw the outline", action="store_true"),
  make_option(c("-H", "--hull"), type="numeric", default=0.25,
              help="distance to additional points to enlarge hull [default = 0.25]",
              metavar="N"),
  make_option(c("-A", "--ahull"), type="numeric", default=5,
              help="alpha value to use for the alpha-convex hull calculation. 0 corresponds to a standard convex hull [default = 5]",
              metavar="N"),
  make_option(c("-b", "--blackbg"), type="logical", default=FALSE,
              help="create plots with a black background", action="store_true"),
  make_option(c("", "--abs"), type="logical", default=FALSE,
              help="plot absolute values when not using --one", action="store_true"),
  make_option(c("-t", "--transpose"), type="logical", default=FALSE,
              help="ensure that genes are rows and columns are spots",
              action="store_true"),
  make_option(c("", "--pal"), type="character", default="offwhite.to.black",
              help="palette to use. available: green.to.blue, spectral, offwhite.to.black [default = offwhite.to.black]",
              action="store"),
  make_option(c("-o", "--out"), type="character", default="vtess.pdf",
              help="specify output path [default = vtess.pdf]",
              action="store"),
  make_option(c("-i", "--invert"), type="logical", default=FALSE,
              help="invert palette in output plots (if number of columns is more than 3)",
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

one.pic.mode = opt$options$one
plot.order.col = opt$options$order
ncols = opt$options$columns
split.plots = opt$options$split
black.bg = opt$options$blackbg
transpose = opt$options$transpose
pal.choice = opt$options$pal
outpath = opt$options$out
bwinv = opt$options$invert
sdims = opt$options$selectdims
draw.border = opt$options$border
convhull.distance = opt$options$convhull
relative.frequency = !opt$options$abs & !one.pic.mode
margin.offset = opt$options$margin

palettes <- list(
  discrete = brewer.pal(12, "Paired"),
  green.to.blue = colorRamp(brewer.pal(9,"GnBu")),
  the.cols = colorRamp(c(rgb(255,255,217, maxColorValue=255),
    rgb(65,182,196, maxColorValue=255),
    rgb(8, 29, 88, maxColorValue=255)),
    space="Lab"),
  spectral = colorRamp(brewer.pal(9,"Spectral")),
  offwhite.to.black = colorRamp(c(rgb(220,220,220, maxColorValue=255),
      rgb(0, 0, 0, maxColorValue=255)),
    space="Lab")
)

palette <- palettes[[pal.choice]]

if (!is.null(sdims)) {
  ind <- as.integer(strsplit(sdims, split = "x")[[1]])
  stopifnot(length(ind) == 3)
}

paths = opt$args

st.load.matrix <- function(path)
  as.matrix(read.table(file = path, sep="\t", header = T, row.names = 1))

parse.coords <- function(n, delimiter=coord.delimiter)
  apply(do.call(rbind, strsplit(n, split = "x")), 2, as.numeric)

# Load files
d <- list()
coords <- list()
for(path in paths) {
  d[[path]] = st.load.matrix(path)
  if (transpose) {
    d[[path]] <- t(d[[path]])
  }
  if (!is.null(sdims)) {
    d[[path]] <- d[[path]][, ind]
  }
  coords[[path]] <- parse.coords(rownames(d[[path]]))
  if (relative.frequency && pal.choice != "discrete")
    d[[path]] = prop.table(d[[path]], 1)
}

ncolumns = ncol(d[[1]])
# Set dimensions of plot
if(one.pic.mode | plot.order.col) {
  n = length(d)
} else {
  n = ncolumns
}

if (!is.null(ncols)) {
  nc = ncols
  nr = ceiling(n/ncols)
} else {
  nc = ceiling(sqrt(n))
  nr = ceiling(n/nc)
}

if (split.plots) {
  nc = 1
  nr = 1
}

mins = apply(sapply(d, function(x) apply(x, 2, min)), 1, min)
maxs = apply(sapply(d, function(x) apply(x, 2, max)), 1, max)
ranges = maxs - mins

enlarged.convex.hull <- function(x, y, a) {
  X <- c(x, x + a, x, x - a)
  Y <- c(y + a, y, y - a, y)
  pts = cbind(X,Y)
  pts = pts + rnorm(length(pts), 0, 0.05)
  pts = pts[!duplicated(pts),]
  if (use.convex.hull) {
    cpi <- chull(X, Y)
    list(x = X[cpi], y = Y[cpi])
  } else {
    hull <- ahull(pts, alpha=opt$options$ahull)
    indx=hull$arcs[,"end1"]
    list(x = X[indx], y=Y[indx])
  }
}

vtess <- list()
tiles <- list()
z <- list()
cp <- list()

for(path in paths) {
  print(path)
  x <- coords[[path]][, 1]
  y <- coords[[path]][, 2]
  vtess <- deldir(x, y)
  tiles[[path]] <- tile.list(vtess)
  z[[path]] <- d[[path]]
  if (pal.choice != "discrete")
    z[[path]] <- t((t(z[[path]]) - mins) / ranges)

  cp[[path]] <- enlarged.convex.hull(x, y, convhull.distance)
}

pl <- NULL
if(one.pic.mode) {
  pl <- function(fn) {
    par(mfrow = c(nr, nc))
    for(path in paths)
      fn(path, NULL)
  }
} else {
  if(plot.order.col) {
    pl <- function(fn) {
      for(col in 1:ncolumns) {
        par(mfrow = c(nr, nc))
        for(path in paths)
          fn(path, col)
      }
    }
  } else {
    pl <- function(fn) {
      for(path in paths) {
        par(mfrow = c(nr, nc))
        for(col in 1:ncolumns)
          fn(path, col)
      }
    }
  }
}

plot.asp = (ARRAY.ROWS + 2 * margin.offset) / (ARRAY.COLS + 2 * margin.offset)

pdf(file = outpath, width = 6*nc, height = 6*nr*plot.asp)
par(mar = c(0, 0, 0, 0))
if(black.bg)
  par(bg="black")

make.plot = function(path, col) {
  print(c(path, col))
  x <- coords[[path]][, 1]
  y <- coords[[path]][, 2]
  plot(x, y, type = 'n', bty = 'none', axes = FALSE,
       xlim = c(1 + 1 - (margin.offset + 0.5), 1 + ARRAY.COLS + (margin.offset + 0.5)),
       ylim = c(1 + ARRAY.ROWS + (margin.offset + 0.5), 1 + 1 - (margin.offset + 0.5)),
       asp=1)
  if (add.grid) {
    abline(h=2:(ARRAY.ROWS+1))
    abline(h=c(1, ARRAY.ROWS+2), col='red', lwd=2)
    abline(v=2:(ARRAY.COLS+1))
    abline(v=c(1, ARRAY.COLS+2), col='blue', lwd=2)
  }
  if(one.pic.mode) {
    while (ncol(z[[path]]) < 3)
      z[[path]] <- cbind(z[[path]], 0)
    plot(tiles[[path]], fillcol = rgb(z[[path]][,1:3]), showpoints = FALSE,
         border = ifelse(draw.border, TRUE, NA), clipp = cp[[path]], add = TRUE)
  } else {
    v <- z[[path]][, col]
    if (bwinv) {
      v <- 1 - v
    }
    if (pal.choice == "discrete") {
      stopifnot(all(v == as.integer(v)), min(v) >= 1, max(v) <= length(palette))
      cols <- palette[v]
    } else {
      cols <- rgb(palette(v), maxColorValue=255)
    }
    plot(tiles[[path]], fillcol = cols, showpoints = FALSE, border = ifelse(draw.border, TRUE, NA), clipp = cp[[path]], add = TRUE)
  }
  if (opt$options$outline)
    lines(rbind(cp[[path]],cp[[path]][1,]), col='black')
  if (opt$options$points)
    points(coords[[path]], col='black', cex=2)
}

pl(make.plot)
dev.off()
