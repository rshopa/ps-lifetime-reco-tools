############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

setBestFitEnvironment <- function( LMAParams )
{
  # initialise/copy routine vars (horten names)
  .tau.rng      <- LMAParams[["tau_rng_ns"]]
  .tau.min.frac <- LMAParams[["min_frac_taus"]]
  .sd.rng       <- LMAParams[["sigma_rng_ns"]]
  .int.min.pPs  <- LMAParams[["min_int_pps"]]
  .int.frac.rng <- LMAParams[["ints_frac_range"]]
  .kde.pars.lst <- LMAParams[["kde_params"]]
  .x.tol        <- LMAParams[["fit_control"]][["x_tol"]]
  .rnd.prec     <- as.integer( LMAParams[["fit_control"]][["round_precision"]] )
  
  ..sd.rng.outliers <- 3.5 # for filtering out outliers in KDE ("deeply" private)
  
  # -----------------------------------------------------------------------
  # Default variables:
  # Sd(1), Mu(2), Idir(3), IpPs(4), IoPs(5), tau_dir(6), tau_oPs(7), Bgr(8)
  # Full table:
  # Sd(1), Mu(2), Idir(3), IpPs(4), IoPs(5), tau_dir(6), tau_pPs(7), tau_oPs(8), Bgr(9)
  # -----------------------------------------------------------------------
  
  
  # -------------------------- reordering --------------------------------------
  # variables:
  # Sd(1), Mu(2), Ishrt(3), IoPs(4), tau_shrt(5), tau_oPs(6), Bgr(7)
  ..reorderTaus2Comp <- function(FitPars)
  {
    ord          <- order(FitPars[5:6])
    Is_ordered   <- FitPars[3:4][ord]
    taus_ordered <- FitPars[5:6][ord]
    # intensity parameters, excluding IpPs
    return( c(FitPars[1:2], Is_ordered, taus_ordered, FitPars[7]) )
  }
  # Sd(1), Mu(2), Idir(3), IpPs(4), IoPs(5), tau_dir(6), tau_oPs(7), Bgr(8)
  ..reorderTaus3Comp <- function(FitPars)
  {
    ord          <- order(FitPars[6:7])
    Is_ordered   <- FitPars[c(3,5)][ord]
    taus_ordered <- FitPars[6:7][ord]
    # intensity parameters, excluding IpPs
    return( c(FitPars[1:2], Is_ordered[1], FitPars[4], Is_ordered[2], taus_ordered, FitPars[8]) )
  }
  # returns 0 (NA), 1 (1 solution), 2 (many solutions), 3 (many but < MinNSolutions)
  ..getSolutionCode <- function( Tab, MinNSolutions )
  {
    n.rows <- nrow(Tab)
    if( is.null(n.rows) ) return( 1 )
    if( n.rows == 0 ) return( 0 )
    else return( if(n.rows < MinNSolutions) 3 else 2 )
  }
  # ------------------------- KDE helpers --------------------------------------
  ..normFitTab <- function(TAB)
  {
    maxs  <- apply(TAB, 2, max)
    mins  <- apply(TAB, 2, min)
    norms <- sapply(1:ncol(TAB), 
                    function(i) (TAB[,i] - mins[i]) / (maxs[i] - mins[i]) )
    return( list( TabNorm = norms, Mins = mins, Maxs = maxs ) )
  }
  ..renormFitVec <- function(Vec, Mins, Maxs) Vec * ( Maxs - Mins ) + Mins
  # Z-score, see sqlpad.io/tutorial/remove-outliers/
  ..noOutliersFactVector <- function(V)
  {
    # if median (from LMAParams[["kde_params"]][["drop_outliers_by_median"]])
    if( .kde.pars.lst[["drop_outliers_by_median"]] ) z <- (V - median(V)) / sd(V)
    else z <- (V - mean(V)) / sd(V)
    return( abs(z) <= ..sd.rng.outliers )
  }
  ..removeOutliersFromTab <- function(Tab)
  {
    fac.no.outliers <- apply( apply(Tab, 2, ..noOutliersFactVector), 1, all )
    return( Tab[fac.no.outliers,] )
  }
  ..removeOutliersFromTabByBCols <- function(Tab, col.num = NULL)
  {
    if(is.null(col.num)) col.num <- ncol( Tab )
    fac.no.outliers <- apply( apply(Tab[,1:col.num], 2, 
                                    ..noOutliersFactVector), 1, all )
    return( Tab[fac.no.outliers,] )
  }
  ..getVarRanges <- function(Tab, MEDIAN = FALSE){
    sds <- apply(Tab, 2, sd)
    mid <- ( if(MEDIAN) apply(Tab, 2, median) else apply(Tab, 2, mean) )
    rbind( mid - ..sd.rng.outliers * sds, mid + ..sd.rng.outliers * sds )
  }
  ..setGridInRanges <- function(VarRng, NptsMax)
  {
    max.diff <- max(apply(VarRng,2,diff))
    dx <- max.diff / NptsMax
    lapply(1:ncol(VarRng), function(k) seq(VarRng[1,k], VarRng[2,k], by = dx))
  }
  ..melt2Vars <- function(v1, v2)
  {
    cbind( rep(v1, length(v2)),
           rep(v2, each= length(v1)) )
  }
  ..melt3Vars <- function(v1, v2, v3)
  {
    len1 <- length(v1)
    len2 <- length(v2)
    len3 <- length(v3)
    cbind(rep(v1, len2*len3),
          rep(rep(v2, each = len1), len3),
          rep(v3, each = len1*len2) )
  }
  ..melt4Vars <- function(v1, v2, v3, v4)
  {
    len1 <- length(v1)
    len2 <- length(v2)
    len3 <- length(v3)
    len4 <- length(v4)
    cbind(rep(v1, len2*len3*len4),
          rep(rep(v2, each = len1),len3 * len4),
          rep(rep(v3, each = len1*len2), len4),
          rep(v4, each = len1*len2*len3) )
  }
  
  # ---------------------------- sorting ---------------------------------------
  # FitTab has 8 columns:
  # RSE(1), Sd(2), Mu(3), Ishrt(4), IoPs(5), tau_shrt(6), tau_oPs(7), Bgr(8)
  sortAndFilt2CompFits <- function( FitTab8Cols, 
                                    AddWeight   = TRUE,
                                    ReorderTaus = TRUE,
                                    Verbose     = TRUE ){
    constraint <- apply(FitTab8Cols, 1, function(l)
    {
      if( !is.na(l[1]) )
      {
        # drop RSE column (first one)
        l.cp <- if(ReorderTaus) ..reorderTaus2Comp(l[2:8]) else l[2:8]
        # check only Sd, not Mu
        F.SIGMA <- l.cp[1] > .sd.rng[1] & l.cp[1] < .sd.rng[2]
        # intensities
        F.INTS <- l.cp[3] > 0 & l.cp[4] > 0
        if( F.INTS )
        {
          i.frac <- l.cp[3] / l.cp[4]   # about 70% / 30%
          F.INTS <- F.INTS & i.frac > .int.frac.rng[1] & i.frac < .int.frac.rng[2]
        }
        F.BGR <- l.cp[7] >= -.x.tol
        # tau's
        # simplified check - tau_short > min, tau_oPs < max
        # no need for more: should be reordered
        F.TAU.RNG  <- l.cp[5] > .tau.rng[1] & l.cp[6] < .tau.rng[2] 
        F.TAU.FRAC <- l.cp[6]/l.cp[5] > .tau.min.frac               # fraction between taus
        return( F.SIGMA & F.INTS * F.BGR * F.TAU.RNG & F.TAU.FRAC )
      } else return(FALSE)
    })
    n.solutions <- sum(constraint)
    if(Verbose) cat("N solutions = ", n.solutions, "\n")
    # out template
    out.lst <- list( NSolutions = n.solutions,
                     FitFraction = n.solutions / length(constraint), 
                     TabFiltered = NA )
    if(n.solutions > 0)
    {
      # returns 7 vars and weights:
      # Sd(1), Mu(2), Ishrt(3), IoPs(4), tau_shrt(5), tau_oPs(6), Bgr(7), RSE_wgt(8)
      RSE.weight <- rep(1, n.solutions)
      if(AddWeight) RSE.weight <- 1 / FitTab8Cols[constraint, 1]^2
      # check if 1 solution or many (vector or tab) and reorder
      if( n.solutions > 1 )
        out.lst[["TabFiltered"]] <- 
          cbind( t(apply(FitTab8Cols[constraint, 2:8], 1, ..reorderTaus2Comp)), 
                 RSE.weight )
      else out.lst[["TabFiltered"]] <- 
          c( ..reorderTaus2Comp(FitTab8Cols[constraint, 2:8]), RSE.weight )
    }
    return( out.lst )
  }
  # FitTab has 9 columns:
  # RSE(1), Sd(2), Mu(3), Idir(4), IpPs(5), IoPs(6), tau_dir(7), tau_oPs(8), Bgr(9)
  sortAndFilt3CompFits <- function( FitTab9Cols, 
                                    FixedTauDir = FALSE, # for 3-component
                                    AddWeight   = TRUE,
                                    ReorderTaus = TRUE,  # forces the order that tau_dir < tau_oPs
                                    Verbose     = TRUE ){
    constraint <- apply(FitTab9Cols, 1, function(l)
    {
      if( !is.na(l[1]) )
      {
        # drop RSE column
        l.cp <- if(FixedTauDir | !ReorderTaus) l[2:9] else ..reorderTaus3Comp(l[2:9])
        # check only Sd, not Mu
        F.SIGMA <- l.cp[1] > .sd.rng[1] & l.cp[1] < .sd.rng[2]
        # intensities
        F.INTS <- l.cp[3] > 0 & #l.cp[4] > 0 & 
          l.cp[5] > 0
        if( F.INTS ){
          i.pPs.norm <- l.cp[4] / sum(l.cp[3:5])
          i.frac     <- l.cp[3] / l.cp[5]
          F.INTS <- F.INTS & i.pPs.norm > .int.min.pPs & 
            i.frac > .int.frac.rng[1] & i.frac < .int.frac.rng[2]
        }
        F.BGR <- l.cp[8] >= -.x.tol
        # tau's
        F.TAU.RNG  <- l.cp[7] > .tau.rng[1] & l.cp[7] < .tau.rng[2]
        F.TAU.FRAC <- l.cp[7]/l.cp[6] > .tau.min.frac # fraction between taus
        if(!FixedTauDir) # add check for tau_dir
          F.TAU.RNG <- F.TAU.RNG & l.cp[6] > .tau.rng[1] & l.cp[6] < .tau.rng[2]
        return( F.SIGMA & F.INTS * F.BGR * F.TAU.RNG & F.TAU.FRAC )
      } else return(FALSE)
    })
    n.solutions <- sum(constraint)
    if(Verbose) cat("N solutions = ", n.solutions, "\n")
    # out template
    out.lst <- list( NSolutions = n.solutions,
                     FitFraction = n.solutions / length(constraint), 
                     TabFiltered = NA )
    if(n.solutions > 0)
    {
      # returns 8 vars and weights:
      # Sd(1), Mu(2), Idir(3), IpPs(4), IoPs(5), tau_dir(6), tau_oPs(7), Bgr(8), RSE_wgt(9)
      RSE.weight <- rep(1, n.solutions)
      if(AddWeight) RSE.weight <- 1 / FitTab9Cols[constraint, 1]^2
      # check if 1 solution or many (vector or tab)
      if( n.solutions > 1)
      {
        if(ReorderTaus) out.lst[["TabFiltered"]] <- 
            cbind( t(apply(FitTab9Cols[constraint, 2:9], 1, ..reorderTaus3Comp)),
                   RSE.weight )
        else out.lst[["TabFiltered"]] <- cbind( FitTab9Cols[constraint, 2:9], RSE.weight )
      } else {
        if(ReorderTaus) out.lst[["TabFiltered"]] <- 
            c( ..reorderTaus3Comp(FitTab9Cols[constraint, 2:9]), RSE.weight )
        else out.lst[["TabFiltered"]] <- c( FitTab9Cols[constraint, 2:9], RSE.weight )
      }
    }
    return( out.lst )
  }
  # TODO: probably redundant
  # sortFitTabByNComponents <- function( FitTabFullParams,
  #                                      ExcludePPs  = TRUE,
  #                                      FixedTauDir = TRUE,
  #                                      Verbose     = TRUE )
  # {
  #   n.cols <- ncol(FitTabFullParams)
  #   if( n.cols == 8 ) # for 2-component model
  #     return( sortAndFilt2CompFits( FitTabFullParams,
  #                                   AddWeight   = TRUE,
  #                                   ReorderTaus = TRUE,
  #                                   Verbose     = Verbose ) )
  #   else {
  #     out.cols <- if(ExcludePPs) c(1:7, 9:10) else (1:n.cols)
  #     return( sortAndFilt3CompFits( FitTabFullParams[,out.cols],          # exclude tau_pPs?
  #                                   FixedTauDir = FixedTauDir,                   # for 3-component
  #                                   AddWeight   = TRUE,
  #                                   ReorderTaus = TRUE,
  #                                   Verbose     = Verbose ) )
  #   }
  # }
  # for Sd, Mu [and Tau_short]
  sortFilt2Or3Cols <- function( TabParams,
                                AddWeight    = TRUE,
                                Verbose      = FALSE )
  {
    n.cols    <- ncol(TabParams)
    n.entries <- nrow(TabParams)
    # vectors of c(RSE, Sd, Mu, [Tau_short])
    # factor for tau_short
    F.FINITE <- apply(apply(TabParams, 2, is.finite), 1, all)
    F.TAU    <- rep(TRUE, n.entries)
    if( n.cols == 4 ) F.TAU  <- TabParams[,4] > .tau.rng[1] & TabParams[,4] < .tau.rng[2] 
    F.SIGMA  <- TabParams[,2] > .sd.rng[1] & TabParams[,2] < .sd.rng[2]
    F.CONSTR <- F.FINITE & F.TAU & F.SIGMA
    n.solutions <- sum(F.CONSTR)
    if(Verbose) cat("N solutions =", n.solutions, "\n")
    # out template
    out.lst <- list( NSolutions  = n.solutions,
                     FitFraction = n.solutions / n.entries, 
                     TabFiltered = NA )
    if(n.solutions > 0)
    {
      # returns 7 vars and weights:
      # Sd(1), Mu(2), Ishrt(3), IoPs(4), tau_shrt(5), tau_oPs(6), Bgr(7), RSE_wgt(8)
      RSE.weight <- rep(1, n.solutions)
      if(AddWeight) RSE.weight <- 1 / TabParams[F.CONSTR, 1]^2
      # check if 1 solution or many (vector or tab)
      if( n.solutions > 1)
        out.lst[["TabFiltered"]] <- cbind( TabParams[F.CONSTR, 2:n.cols], RSE.weight )
      else out.lst[["TabFiltered"]] <- c( TabParams[F.CONSTR, 2:n.cols], RSE.weight )
    }
    return( out.lst )
  }
  
  # ------------------- MAIN FUNCTIONS: fing best fit --------------------------
  # for Sd, Mu [and Tau_short]
  findQuickBestByRSEFor2Or3Cols <- function( FitTabParams,
                                             Verbose = FALSE )
  {
    if(Verbose) cat("Quick find by RSE:\t")
    n.cols       <- ncol(FitTabParams)
    filtered.lst <- sortFilt2Or3Cols( FitTabParams, AddWeight = TRUE, Verbose = Verbose )
    if( filtered.lst[["NSolutions"]] > 0 )
    {
      # check if 1 solution or many (vector or tab)
      if( filtered.lst[["NSolutions"]] > 1)
        return( as.numeric( 
          filtered.lst[["TabFiltered"]][which.max(filtered.lst[["TabFiltered"]][, n.cols]), 1:(n.cols-1)] 
        ) )
      else return( as.numeric( filtered.lst[["TabFiltered"]][1:(n.cols-1)] ) )
    }
    else return( NA )
  }
  
  findBestByRSEFromFilteredTab <- function( ListFiltTabWWgt )                 # optional: precision to round
  {
    n.cols <- NA
    if( ListFiltTabWWgt[["NSolutions"]] == 1 ) n.cols <- length(ListFiltTabWWgt[["TabFiltered"]]) # flag for 2- or 3-component
    if( ListFiltTabWWgt[["NSolutions"]] > 1 )  n.cols <- ncol(ListFiltTabWWgt[["TabFiltered"]]) # flag for 2- or 3-component
    # dummy template
    out.list <- list( FitFraction = ListFiltTabWWgt[["FitFraction"]], 
                      BestFit     = NA, 
                      MinRSE      = NA )
    if( ListFiltTabWWgt[["NSolutions"]] > 0 )
    {
      if( ListFiltTabWWgt[["NSolutions"]] == 1 )
      {
        best.WGT <- ListFiltTabWWgt[["TabFiltered"]][n.cols]
        best.out <- as.numeric(ListFiltTabWWgt[["TabFiltered"]][1:(n.cols-1)])
      }
      else {
        best.id  <- which.max(ListFiltTabWWgt[["TabFiltered"]][,n.cols])
        best.WGT <- ListFiltTabWWgt[["TabFiltered"]][,n.cols][best.id] # reverse from weight
        best.out <- as.numeric(ListFiltTabWWgt[["TabFiltered"]][best.id, 1:(n.cols-1)])
      }
      if( n.cols == 9 ) out.list[["BestFit"]] <- round( ..reorderTaus3Comp(best.out), .rnd.prec )
      else out.list[["BestFit"]] <- round( ..reorderTaus2Comp(best.out), .rnd.prec )
      out.list[["MinRSE"]] <- sqrt( 1 / best.WGT )
    }
    return( out.list )  
  }
  
  
  # --------------------------- KDE section ------------------------------------
  # info helpers
  
  ..infoKDE <- function( KDEobj, Mode = NULL )
  {
    if( Mode == "Tab" ){
      cat("Quick find by KDE:\t")
      if( is.na(KDEobj[1]) )
      {
        cat("no solutions, returns NA.\n")
        return(NA)
      } else cat( "N solutions =", nrow(KDEobj), "\n" )
    }
    if( Mode == "Hpi" )
    {
      cat("Hpi Matrix (dscalar):\n")
      print( KDEobj )
    }
  }
  # eg. c(Sd, Mu, Tau_shrt, RSE.wgt)
  findBestFitKDEFor3ParamTab <- function( TabFiltWWgts,
                                          Verbose = TRUE )
  {
    if( Verbose ) ..infoKDE( TabFiltWWgts, Mode = "Tab" )
    kde.weights   <- TabFiltWWgts[,4] # RSE.wgt
    lst.fits.norm <- ..normFitTab( ..removeOutliersFromTab( TabFiltWWgts[,1:3] ) )
    
    if( nrow(lst.fits.norm[[1]]) < .kde.pars.lst[["min_n_solutions"]] ){
      cat("too few points, returns NA.\n")
      return(NA)
    }
    # the rest (Ps decay params), SLOW!!!
    kde.tab  <- lst.fits.norm[[1]] # PI vars, depending on nr of components
    H.dsc    <- Hpi( kde.tab, pilot = "dscalar" )
    if( Verbose ) ..infoKDE( H.dsc, Mode = "Hpi" )
    kde.vars <- ..setGridInRanges(..getVarRanges(kde.tab), 
                                  NptsMax = .kde.pars.lst[["max_axis_length_3d"]])
    kde.vars <- ..melt3Vars( kde.vars[[1]], kde.vars[[2]], kde.vars[[3]] )        # convert to melted tab (see library(reshape2) )
    est.kde  <- estimate_kde_for_pts_wgt( kde.vars, kde.tab, H.dsc, kde.weights ) # from kdeCPP.cpp
    best.pi  <- kde.vars[which.max(est.kde),] # 4 vars
    # renormalise and merge
    return( ..renormFitVec( best.pi, lst.fits.norm[[2]], lst.fits.norm[[3]] ) )
  }
  
  # eg. c(Sd, Mu, RSE.wgt)
  findBestFitKDEFor2ParamTab <- function( TabFiltWWgts,
                                          Verbose = TRUE )
  {
    if(Verbose) ..infoKDE( TabFiltWWgts, Mode = "Tab" )
    kde.weights   <- TabFiltWWgts[,3]
    lst.fits.norm <- ..normFitTab( ..removeOutliersFromTab( TabFiltWWgts[,1:2] ) )
    if( nrow(lst.fits.norm[[1]]) < .kde.pars.lst[["min_n_solutions"]] ){
      cat("too few points, returns NA.\n")
      return(NA)
    }
    # the rest (Ps decay params), SLOW!!!
    kde.tab  <- lst.fits.norm[[1]] # PI vars, depending on nr of components
    H.dsc    <- Hpi( kde.tab, pilot = "dscalar" )
    if( Verbose ) ..infoKDE( H.dsc, Mode = "Hpi" )
    kde.vars <- ..setGridInRanges(..getVarRanges(kde.tab), 
                                  NptsMax = .kde.pars.lst[["max_axis_length_2d"]])
    kde.vars <- ..melt2Vars( kde.vars[[1]], kde.vars[[2]] ) # convert to melted tab (see library(reshape2) )
    est.kde  <- estimate_kde_for_pts_wgt( kde.vars, kde.tab, H.dsc, kde.weights )
    best.pi  <- kde.vars[which.max(est.kde),] # 4 vars
    # renormalise and merge
    return( ..renormFitVec( best.pi, lst.fits.norm[[2]], lst.fits.norm[[3]] ) )
  }

  rm( list = c("LMAParams") )
  return( environment() )
}
