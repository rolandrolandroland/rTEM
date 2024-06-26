% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rTEM_multimersim_inhibit.R
\name{create_groups_dev}
\alias{create_groups_dev}
\title{Create Multimers around guest type points}
\usage{
create_groups_dev(
  num_neighbors = 6,
  upp_guest,
  upp_host,
  probs = c(1, 0, 0, 0, 0, 0),
  weights = c(1, 1, 1),
  trans_plane = c("x", "y"),
  trans_frac = 0.5,
  sample_method = "rank",
  group_size = 2,
  exponent = 1
)
}
\arguments{
\item{num_neighbors}{the number of nearest neighbors around each guest type point to
to be considered consider.}

\item{upp_guest}{The guest type points that groups will be created around
Currently, the points in `upp_guest` should not also be part of `upp_host`}

\item{upp_host}{The host type points that the multimers will be selected from}

\item{probs}{vector of probabilities.
For \eqn{probs = c(p_1, p_2, p_3, p_4)}, the probability of the first NN being
selected in \eqn{p_1}, the probability of the second is \eqn{p_2}, and so on}

\item{weights}{vector of length equal to number of dimensions in \emph{upp}.
the weighted distance to each of \emph{num_neighbors} nearest neighbors
is calculated using \eqn{\sqrt{w_1 x^2 + w_2 y^2 + w_3 z^2}},
 where \emph{weights} = (\eqn{w_1, w_2, w_3}). Set to \emph{c(1, 1, 0)} for vertical dimers.
 @param trans_plane plane for translating
 @param trans_frac fraction of distance in plane `trans_plane` to reduce distance by}

\item{sample_method}{if equal to \emph{"rank"}, the probability of a point of rank \emph{x}
being chosen as a guest is \emph{probs[x]}.  If equal to  \emph{"exp"},
the probability of a point of rank \emph{x} being chosen
 as a guest is \emph{probs[x] * exp(-exponent * distances[ranks]))}}

\item{group_size}{a numeric. How large will the groups be}

\item{exponent}{a numeric. If \emph{sample_method = "exp"}, then
this is the value of \emph{exponent} in \emph{exp(-exponent * distances[ranks])}}
}
\description{
This is the function that is called within \code{\link{multimersim}}
to create groups from single points.
}
\details{
Algorithm Steps:
\itemize{
\item{Step 1:} {Calculate `num_neighbors` worth of nearest points in `upp_host` to each
point in `upp_guest`}
\item{Step 2:} Use \code{\link{get_ranks}} to rank each of these neighbors by weighted distance to their
respective guest point.  Note: If weights = c(1, 1, 1) for 3D patterns or c(1, 1) for 2D patterns,
then they will be ranked by distance.
\item{Step 3:} {Use \code{\link{pick_neighbor}} function to pick which of the neighbors will become
guest type.}
\item{Step 4:} {Deal with duplicates.  If `upp_guest` contains points that are close enough to each
other, it is possible they will share nearest neighbors.  If any point is selected twice, then
another shall be selected to maintain the proper number of guests}
\item{Step 5:} {Translate.  All selected neighbors are translated in the `trans_plane` plane so that the
distance to their respective guest point in the `trans_plane` plane is reduced by a fraction of `trans_frac` }}
}
