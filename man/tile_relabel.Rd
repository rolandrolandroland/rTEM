% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rTEM_tiler.R
\name{tile_relabel}
\alias{tile_relabel}
\title{Tile and then relabel}
\usage{
tile_relabel(
  seed,
  pp3_full,
  funcs = c("K", "G", "F", "GXGH"),
  ...,
  host_formula,
  dopant_formula,
  x_tiles = 1,
  y_tiles = 1,
  z_tiles = 1,
  x_trim = 0,
  y_trim = 0,
  z_trim = 0,
  maxKr = 10,
  nKr = 200,
  maxGr = 5,
  nGr = 1000,
  maxGXGHr = 3,
  maxGXHGr = 8,
  nGXr = 1000,
  vside = 0.3,
  K_cor = "border",
  G_cor = "rs",
  F_cor = "rs",
  GXGH_cor = "rs",
  GXHG_cor = "rs"
)
}
\arguments{
\item{seed}{input for `set.seed` function (for reproducibility)}

\item{pp3_full}{input pp3 pattern for relabeling and tiling}

\item{funcs}{summary functions to calculate on relabeling}

\item{host_formula, dopant_formula}{formula for host and dopant marks.  Summary functions
will be calculated from dopant type points}

\item{x_tiles, y_tiles, z_tiles}{number of tiles in each dimension (set equal to 1 for no change)}

\item{z_trim, y_trim, z_trim}{amount to trim off of each dimension.  Trimming will be applied equall to both sides}

\item{K_cor}{edge correction(s) to be used for K function}

\item{G_cor}{edge correction(s) to be used for G function}

\item{F_cor}{edge correction(s) to be used for F function}

\item{GXGH_cor}{edge correction(s) to be used for GXGH function}

\item{GXHG_cor}{edge correction(s) to be used for GXHG function}
}
\description{
This function takes an input pp3 point pattern, relabels it, then
tiles the relabeled point pattern.  Then it calculates summary functions on it
}
