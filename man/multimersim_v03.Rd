% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rTEM_deprecated.R
\name{multimersim_v03}
\alias{multimersim_v03}
\title{Multimer Simulation Version 03}
\usage{
multimersim_v03(
  guest_pattern = NULL,
  upp,
  output = "guest pattern type",
  n_guests = NA,
  min_thick = NA,
  max_thick = NA,
  ztrim = 0,
  group_size = 2,
  num_neighbors = 6,
  sample_method = "rank",
  weights = c(1, 1, 1),
  probs = c(1, 0, 0, 0),
  intensity_upp = NA
)
}
\arguments{
\item{guest_pattern}{point pattern of class \emph{ppp} or \emph{pp3}.  The final multtimer
pattern will match this pattern in class, intensity, and domain.  If this is left as NULL, then
the domain will match that of \emph{upp}, will be of the class specified in \emph{output},
and have number of guests specified in \emph{n_guest}}

\item{upp}{point pattern of class \emph{ppp} or \emph{pp3} to use as underlying point pattern.
Multimers will be selected from these points.}

\item{output}{leave as \emph{"guest pattern type"} if specifying \emph{guest_pattern}.  Otherwise, set to \emph{"ppp"} for
2D output or \emph{pp3} for 3D output}

\item{n_guests}{leave as \emph{NA} if specifying \emph{UPP}. Otherwise, set to
integer for number of guests in multimer pattern}

\item{min_thick}{if \emph{guest_pattern} is 2d (ppp) and \emph{upp} is 3d (pp3) this
determines the smallest z value to keep before collapsing into 2d.}

\item{max_thick}{if \emph{guest_pattern} is 2d (ppp) and \emph{upp} is 3d (pp3) this
determines the largest z value to keep before collapsing into 2d.}

\item{ztrim}{a numeric.  This will trim this amount from both top and bottom (positive and negative z)
after multimers are generated and before pp3 pattern is collapsed into ppp.
Only applies if \emph{upp} is 3D (\emph{pp3})}

\item{group_size}{size of clusters.  If equal to 1, then all points will be independent
of each other}

\item{num_neighbors}{number of nearest neighbors to select from when
forming dimers, trimers, etc..}

\item{sample_method}{if equal to \emph{"rank"}, the probability of a point of rank \emph{x}
being chosen as a guest is \emph{probs[x]}.  If equal to  \emph{"exp"},
the probability of a point of rank \emph{x} being chosen
 as a guest is \emph{probs[x] * exp(-distances[ranks]))}}

\item{weights}{vector of length equal to number of dimensions in \emph{upp}. the weighted distance to each of \emph{num_neighbors} nearest neighbors
is calculated using \eqn{\sqrt{w_1 x^2 + w_2 y^2 + w_3 z^2}},
 where \emph{weights} = (\eqn{w_1, w_2, w_3}). Set to \emph{c(1, 1, 0)} for vertical dimers.}

\item{probs}{vector of probabilities.  Should sum to \emph{group_size}-1.
For \eqn{probs = c(p_1, p_2, p_3, p_4)}, the probability of the first NN being
selected in \eqn{p_1}, the probability of the second is \eqn{p_2}, and so on}

\item{intensity_upp}{the \emph{upp} will be rescaled to have this intensity before the
marks are assigned. Leave as \emph{NA} to use \emph{upp} as it is}
}
\description{
Under construction. See \code{\link{multimersim}} for stable version. Simulates multimers (groups of two/dimers, three/trimers, etc.) in
underlying point pattern (upp) to match domain and number of points in \emph{guest_pattern}.
}
\details{
algorithm steps:
 \itemize{
 \item{Step 1:} {rescale \emph{upp} to match \emph{intensity_upp}}
 \item{Step 2:} {Select points in the scaled \emph{upp} that are
 inside the domain of \emph{guest_pattern}. }
 \item{Step 3:} {Determine number of guest groups or clusters (for dimers, this is number of guests / 2)
 and assign this many points in the scaled subsetted UPP to be guests.
  These are the "centroids"}
 \item{Step 4:} {Take the \emph{num_neighbors} closest point to each guest}
 \item{Step 5:} {Rank each neighbor by weighted distance to centroid}
 \item{Step 6:} {Using the probabilities in \emph{probs}, select which neighbors
  are to be guests (so that the cluster size is now equal to \emph{group_size})}
 \item{Step 7:} {For any duplicates, redo process so that correct number of guests are present}
 \item{Step 8:} {If \emph{guest_pattern} is 2D and \emph{UPP} is 3D, remove Z coordinate to
 make new pattern 2D}
}
}
