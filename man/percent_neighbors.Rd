% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nndist_analysis_functions.R
\name{percent_neighbors}
\alias{percent_neighbors}
\title{Fraction of all neighbors within a distance that are same type or fraction of all nth nearest neighbors that are same type}
\usage{
percent_neighbors(
  pattern,
  i,
  window = NULL,
  drop_isolated = FALSE,
  value = "NN",
  rads_vec = 1:25,
  nn_vec = 1:5,
  output = "plot",
  location = getwd(),
  image_name = "image",
  unit = "nm",
  round_to = 3,
  ...
)
}
\arguments{
\item{pattern}{Image of type ppp or pp3}

\item{i}{mark for point type to be used}

\item{value}{Determines output.  If = "NN" then returns fraction of all nearest neighbors in nn_vec that are same type.  If = "distance" then returns fraction of
all points within distances in rads_vec that are same type}

\item{rads_vec}{vector of distances to calculate}

\item{nn_vec}{vector of integers values of nearest neighbors}

\item{output}{set to "plot" (default) if only plotting.  Set to "save" if plotting and saving}

\item{location}{location to save image to if `output = save`.  default is working directory.}

\item{image_name}{Name to save image as.  Will have `_NN_dist_perc.png` added to it if `value` = "percent" or `_NN_perc_dist.png` if `value` = "distance"}

\item{unit}{distance unit name to be added to output name}

\item{round_to}{number of decimal places to round values to}
}
\value{
vector of fractions
}
\description{
Calculate either the average fraction of points in object `image` at distances `rads_vec` from each point of type `i` that are type `i` (`value` = "distance") or the average fraction
of each nearest neighbor in `nn_vec` of each point of type `i` that is type `i` `value` = "NN").
}
