#' Advanced Nearest Neighbour Distance Function G
#' @description Estimates the nearest neighbour distance distribution
#'  Gest from a point pattern in a window of arbitrary shape. Modified from the \emph{spatstat} package
#' function \link[spatstat.explore]{Gest}  to allow for usage of nearest neighbors beyond the first
#' @param k Integer, or integer vector. The algorithm will compute the distance to the kth nearest neighbour.
#' See \link[spatstat.geom]{nndist}
#' @inheritParams spatstat.explore::Gest
#' @md
#' @export
"Gest_nn" <-
  "nearest.neighbour" <-
  function(X, r=NULL, breaks=NULL, k = 1, ..., correction=c("rs", "km", "han"),
           domain=NULL) {
    verifyclass(X, "ppp")
    if(!is.null(domain))
      stopifnot(is.subset.owin(domain, Window(X)))

    ##
    W <- X$window
    npts <- spatstat.geom::npoints(X)
    lambda <- npts/spatstat.geom::area(W)

    ## determine r values
    rmaxdefault <- spatstat.explore::rmax.rule("G", W, lambda)
    breaks <- spatstat.geom::handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
    rvals <- breaks$r
    rmax  <- breaks$max
    zeroes <- numeric(length(rvals))

    ## choose correction(s)
    #  correction.given <- !missing(correction) && !is.null(correction)
    if(is.null(correction)) {
      correction <- c("rs", "km", "han")
    } else correction <- spatstat.geom::pickoption("correction", correction,
                                    c(none="none",
                                      border="rs",
                                      rs="rs",
                                      KM="km",
                                      km="km",
                                      Kaplan="km",
                                      han="han",
                                      Hanisch="han",
                                      cs="han",
                                      ChiuStoyan="han",
                                      best="km"),
                                    multi=TRUE)

    ##  compute nearest neighbour distances
    nnd <- nndist(X$x, X$y, k = k)
    ##  distance to boundary
    bdry <- bdist.points(X)
    ## restrict to subset ?
    if(!is.null(domain)) {
      ok <- inside.owin(X, w=domain)
      nnd <- nnd[ok]
      bdry <- bdry[ok]
    }
    ##  observations
    o <- pmin.int(nnd,bdry)
    ##  censoring indicators
    d <- (nnd <= bdry)

    ## initialise fv object
    df <- data.frame(r=rvals, theo=1-exp(-lambda * pi * rvals^2))
    Z <- fv(df, "r", substitute(G(r), NULL), "theo", . ~ r,
            c(0,rmax),
            c("r", "%s[pois](r)"),
            c("distance argument r", "theoretical Poisson %s"),
            fname="G")

    if("none" %in% correction) {
      ##  UNCORRECTED e.d.f. of nearest neighbour distances: use with care
      if(npts <= 1)
        edf <- zeroes
      else {
        hh <- hist(nnd[nnd <= rmax],breaks=breaks$val,plot=FALSE)$counts
        edf <- cumsum(hh)/length(nnd)
      }
      Z <- bind.fv(Z, data.frame(raw=edf), "hat(%s)[raw](r)",
                   "uncorrected estimate of %s", "raw")
    }
    if("han" %in% correction) {
      if(npts <= 1)
        G <- zeroes
      else {
        ##  uncensored distances
        x <- nnd[d]
        ##  weights
        a <- eroded.areas(W, rvals, subset=domain)
        ## calculate Hanisch estimator
        h <- hist(x[x <= rmax], breaks=breaks$val, plot=FALSE)$counts
        G <- cumsum(h/a)
        G <- G/max(G[is.finite(G)])
      }
      ## add to fv object
      Z <- bind.fv(Z, data.frame(han=G),
                   "hat(%s)[han](r)",
                   "Hanisch estimate of %s",
                   "han")
      ## modify recommended plot range
      attr(Z, "alim") <- range(rvals[G <= 0.9])
    }

    if(any(correction %in% c("rs", "km"))) {
      ## calculate Kaplan-Meier and border correction (Reduced Sample) estimates
      want.rs <- "rs" %in% correction
      want.km <- "km" %in% correction
      if(npts == 0) {
        result <- list(rs=zeroes, km=zeroes, hazard=zeroes, theohaz=zeroes)
      } else {
        result <- km.rs.opt(o, bdry, d, breaks, KM=want.km, RS=want.rs)
        if(want.km)
          result$theohaz <- 2 * pi * lambda * rvals
      }
      wanted <- c(want.rs, rep(want.km, 3L))
      wantednames <- c("rs", "km", "hazard", "theohaz")[wanted]
      result <- as.data.frame(result[wantednames])
      ## add to fv object
      Z <- bind.fv(Z, result,
                   c("hat(%s)[bord](r)", "hat(%s)[km](r)",
                     "hat(h)[km](r)", "h[pois](r)")[wanted],
                   c("border corrected estimate of %s",
                     "Kaplan-Meier estimate of %s",
                     "Kaplan-Meier estimate of hazard function h(r)",
                     "theoretical Poisson hazard function h(r)")[wanted],
                   if(want.km) "km" else "rs")
      ## modify recommended plot range
      attr(Z, "alim") <- with(Z, range(Z$r[Z$km <= 0.9]))
    }
    nama <- names(Z)
    fvnames(Z, ".") <- rev(setdiff(nama, c("r", "hazard", "theohaz")))
    unitname(Z) <- unitname(X)
    return(Z)
  }

