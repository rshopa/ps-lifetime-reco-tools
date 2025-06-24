############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

# WARNING: need to load those libraries and compile OpenMP KDE function
# library(minpack.lm)
# library(parallel)
# library(ks)
# library(jsonlite)
# source("setParamsEnvr.R")
# source("setLMAEnvr.R")
# Rcpp::sourceCpp("kdeCPP.cpp")
# source("setBestFitEnvr.R")

setPsLifetimeLMAFitEnvironment <- function()
{
  .pars.env     <- setParamsEnvr() 
  .lma.env      <- NA
  .best.fit.env <- NA
  
  configureSetupFromJSON <- function( ... )
  {
    args <- list(...)
    if( length(args) == 3 )
    {
      .pars.env[["setParams"]]( args[["HST.PARAMS"]], 
                                args[["INI.GUESS.PARAMS"]],
                                args[["LMA.PARAMS"]] )
    } else .pars.env[["setParamsCombinedJSON"]]( args[[1]] )
    .lma.env      <<- setLMMultiCompFitEnvironment( .pars.env )
    .best.fit.env <<- setBestFitEnvironment( .pars.env[["lma.params"]] )
  }
  
  # defines spectrum and updates tau_dir from matter density (in g/mL, e.g. 1.0 for water)
  .definePsLifetimeSpectrum <- function( DtAxisNs, HistAU, MatterDensitygmL = NULL )
  {
    .lma.env[["defineSpectrumParams"]]( DtAxisNs, HistAU, .pars.env[["dt.hst.params"]] )
    if( !is.null(MatterDensitygmL) ) 
      .lma.env[["p_tau.dir"]] <<- round( 1e-3 * 819.151 * exp(-0.517436 * MatterDensitygmL), 3 )
  }
  # helper to find best fit (makes code more DRY)
  ..getBestFromUnfiltLMATab <- function( LMATab,
                                         AddWeight = TRUE,
                                         Verbose   = TRUE, 
                                         ... ) # optional: ReorderTaus (2-comp) or FixedTauDir (3-comp)
  {
    args <- list(...)
    filt.lst <- if( ncol(LMATab) == 8 )
      .best.fit.env[["sortAndFilt2CompFits"]]( LMATab, AddWeight, args[["ReorderTaus"]], Verbose )
    else 
      # NOTE: for 3-component, tau_pPs is inserted as col no 8, need to be removed
      .best.fit.env[["sortAndFilt3CompFits"]]( LMATab[,c(1:7,9:10)], 
                                               args[["FixedTauDir"]], 
                                               AddWeight, 
                                               args[["ReorderTaus"]], 
                                               Verbose )
    return( .best.fit.env[["findBestByRSEFromFilteredTab"]](filt.lst) )
  }
  # helpers for LMA2 or LMA3 output
  ..calcRSEsForFittedParams <- function( Params, NComponents )
  {
    rse.linear <- NA
    rse.log    <- NA
    if( is.finite(Params[["BestFit"]][1]) )
    {
      fit.pars <- NULL
      if( NComponents == 2 ) fit.pars <- Params[["BestFit"]]
      else fit.pars <- c( Params[["BestFit"]][1:6], 
                          .lma.env[["p_tau.pPs"]], Params[["BestFit"]][7:8] )
      rse.linear <- .lma.env[[".estimateRSELinear"]](fit.pars )
      rse.log    <- .lma.env[[".estimateRSELog"]]( fit.pars )
    }
    return( c(rse.linear, rse.log) )
  }
  ..makeOutListForLMA <- function( DirBest, RSEBest, KDEBest, NComponents )
  {
    rses.dir <- ..calcRSEsForFittedParams( DirBest, NComponents )
    rses.rse <- ..calcRSEsForFittedParams( RSEBest, NComponents )
    rses.kde <- ..calcRSEsForFittedParams( KDEBest, NComponents )
    return( 
      list( byDirectFit = list( par = DirBest[["BestFit"]], 
                                RSElin = rses.dir[1], RSElog = rses.dir[2] ),
            byRSE       = list( par = RSEBest[["BestFit"]],        
                                RSElin = rses.rse[1], RSElog = rses.rse[2] ),
            byKDE       = list( par = KDEBest[["BestFit"]],        
                                RSElin = rses.kde[1], RSElog = rses.kde[2] ) )
    )
  }

  # ---------------------------- MAIN FUNCTIONS --------------------------------
  fitHistByTwoStageLMA2OMP <- function( DtAxisNs, 
                                        HistAU, 
                                        Verbose = TRUE )
  {
    .definePsLifetimeSpectrum( DtAxisNs, HistAU ) # tau_dir not used here
    # first fit, only for Sd, Mu and tau_short/fast: no full span, no errors
    lmfit <- 
      .lma.env[["runSeededLMA2AllFreeLinOMP"]]( .pars.env[["ini.guess.df"]], 
                                                LowThrldDt = .pars.env[["dt.hst.params"]][["low_threshold_for_first_fit_ns"]] ) 
    # filtered list
    fltlst <- .best.fit.env[["sortFilt2Or3Cols"]]( lmfit ) # add weight TRUE by default, no verbose needed here
    # find best outcomes for RSE-based and KDE-besed mathod
    RSEbest <- .best.fit.env[["findQuickBestByRSEFor2Or3Cols"]]( lmfit, Verbose = Verbose )
    KDEbest <- .best.fit.env[["findBestFitKDEFor3ParamTab"]]( fltlst[["TabFiltered"]], Verbose = Verbose )
    
    # second fit (in log scale), using two prevoius outcomes
    lmfit   <- .lma.env[["runSeededLMA2LogOMP"]]( .pars.env[["ini.guess.df"]], RSEbest )
    RSEbest <- ..getBestFromUnfiltLMATab( lmfit, TRUE, Verbose, ReorderTaus = FALSE ) # add RSE weight, no reorder taus
    lmfit   <- .lma.env[["runSeededLMA2LogOMP"]]( .pars.env[["ini.guess.df"]], KDEbest )
    KDEbest <- ..getBestFromUnfiltLMATab( lmfit, TRUE, Verbose, ReorderTaus = FALSE) # add RSE weight, no reorder taus

    # recalculate first stage for direct fit: no lower cut, full span, no errors
    lmfitDir  <- .lma.env[["runSeededLMA2AllFreeLinOMP"]]( .pars.env[["ini.guess.df"]], 
                                                           FullSpan = TRUE, ReturnFullTable = TRUE )
    DIRbest   <- ..getBestFromUnfiltLMATab( lmfitDir, TRUE, Verbose, ReorderTaus = TRUE ) # add RSE weight, reorder taus
    return( ..makeOutListForLMA( DIRbest, RSEbest, KDEbest, 2 )
    )
  }
  
  # 3-component with competing fits:
  # 1 - TO CHECK: (tau_dir is fixed), linear scale over all span
  # 2 - tau_dir is free, 2 OMPs with refinement
  # 3 - (optional) tau_oPs is fixed, the rest refined
  # IRatioMode - Intensity fractions: 0 (else): all free, 1: fix I_oPs / I_pPs = 3, 2: fix all (60 : 10 : 10)
  fitHistByMultiStageLMA3IratioSwitchableOMP <- function(  DtAxisNs, 
                                                           HistAU, 
                                                           IRatioMode         = 0L,     # 0 (else): all free, 1: fix I_oPs / I_pPs = 3, 2: fix all (60 : 10 : 10)
                                                           FixedTauDir        = FALSE,  # TODO: only for modes 0 and 2
                                                           RefineWFixedTauoPs = TRUE,   # TODO: only for modes 0 and 2
                                                           MatterDensitygmL   = NULL,   # if tau_dir is to be modified according to density
                                                           RefinedDtRange     = NULL,   # if arbitrary Dt range is needed for the third fit
                                                           Verbose            = TRUE )
  {
    .definePsLifetimeSpectrum( DtAxisNs, HistAU, MatterDensitygmL ) # tau_dir not used here
    
    # ---------- direct fit: linear scale, full span ----------
    # select proper function and option for reorder taus:
    reorder.taus <- if( IRatioMode == 1L | IRatioMode == 2L ) FALSE else TRUE   # true only for all free
    fit.name     <- "runSeededLMA3AllFreeLinOMP"
    if( IRatioMode == 1L ) fit.name <- "runSeededLMA3AllFreeLinITausRatioFixedOMP"
    if( IRatioMode == 2L ) fit.name <- "runSeededLMA3AllFreeLinIRatioFixedOMP"
    if( FixedTauDir & IRatioMode != 1L )
      fit.name <- if( IRatioMode == 2L ) "runSeededLMA3TauDirLinIRatioFixedOMP"
    else "runSeededLMA3TauDirLinOMP"
    
    lmfit.dir <- .lma.env[[fit.name]]( .pars.env[["ini.guess.df"]],
                                       FullSpan        = TRUE,
                                       LowThrldDt      = NULL,
                                       ReturnFullTable = TRUE ) # no errors (AddErrors = FALSE)
    DIR.best <- .best.fit.env[["findBestByRSEFromFilteredTab"]](
      .best.fit.env[["sortAndFilt3CompFits"]]( lmfit.dir[,c(1:7,9:10)], FALSE, TRUE, reorder.taus, Verbose ) ) # non-fixed Tau_dir, Add RSE weight ...
    
    # ----------------- multi stage fit, tau_dir is free ----------------------
    # I - linear stage
    # re-check and update fit function name 
    if( FixedTauDir & IRatioMode == 0L ) fit.name <- "runSeededLMA3AllFreeLinOMP"
    if( FixedTauDir & IRatioMode == 2L ) fit.name <- "runSeededLMA3AllFreeLinIRatioFixedOMP"
    lmfit <- .lma.env[[fit.name]]( .pars.env[["ini.guess.df"]],
                                   FullSpan        = FALSE, 
                                   LowThrldDt      = NULL, 
                                   ReturnFullTable = FALSE ) # no errors (AddErrors = FALSE)
    # now find the best for RSE-based and KDE-based approach
    RSE.best <- .best.fit.env[["findQuickBestByRSEFor2Or3Cols"]]( lmfit, Verbose = Verbose )   # from fit tab
    filt.lst <- .best.fit.env[["sortFilt2Or3Cols"]]( lmfit, AddWeight = TRUE )                 # Verbose = FALSE
    KDE.best <- .best.fit.env[["findBestFitKDEFor2ParamTab"]]( filt.lst[["TabFiltered"]],
                                                               Verbose = Verbose  ) # from list of sorted
    # II - log stage
    # !!! slow: executed twice !!!
    if( IRatioMode == 1L ) fit.name <- "runSeededLMA3FreeTauDirLogITausRatioFixedOMP"
    else if( IRatioMode == 2L ) fit.name <- "runSeededLMA3FreeTauDirLogIRatioFixedOMP"
    else fit.name <- "runSeededLMA3FreeTauDirLogOMP"
    lmfit.rse <- .lma.env[[fit.name]]( .pars.env[["ini.guess.df"]], BestGauss = RSE.best ) 
    lmfit.kde <- .lma.env[[fit.name]]( .pars.env[["ini.guess.df"]], BestGauss = KDE.best ) 
    # NOTE: RSE and KDE here denote how the first two parameters were acquired (Sd, Mu)
    # Add RSE weight, non-fixed Tau_dir, reorder only for all-free intensity ratios
    RSE.best <- ..getBestFromUnfiltLMATab( lmfit.rse, TRUE, Verbose, FixedTauDir = FALSE, ReorderTaus = reorder.taus ) 
    KDE.best <- ..getBestFromUnfiltLMATab( lmfit.kde, TRUE, Verbose, FixedTauDir = FALSE, ReorderTaus = reorder.taus )
    
    # TODO: add refined fit for mode 1: only I_oPs / I_pPs is fixed
    # III (optional) - refine all parameters with fixed tau_oPs
    if( RefineWFixedTauoPs & IRatioMode != 1L &
        is.finite(RSE.best[["BestFit"]][1]) &
        is.finite(KDE.best[["BestFit"]][1]) ) # ignore if NA's in the previous
    {
      fit.name <- if( IRatioMode == 2L ) "runSeededLMA3TauoPsLinIRatioFixedOMP"
      else "runSeededLMA3TauoPsLinOMP"
      lmfit.rse <- .lma.env[[fit.name]]( .pars.env[["ini.guess.df"]],
                                         TauoPs  = RSE.best[["BestFit"]][7],
                                         DtRange = RefinedDtRange )
      lmfit.kde <- .lma.env[[fit.name]]( .pars.env[["ini.guess.df"]],
                                         TauoPs  = KDE.best[["BestFit"]][7],
                                         DtRange = RefinedDtRange )
      # Add RSE weight, non-fixed Tau_dir, reorder only for all-free intensity ratios
      RSE.bestRef <- ..getBestFromUnfiltLMATab( lmfit.rse, TRUE, Verbose, FixedTauDir = FALSE, ReorderTaus = reorder.taus ) 
      KDE.bestRef <- ..getBestFromUnfiltLMATab( lmfit.kde, TRUE, Verbose, FixedTauDir = FALSE, ReorderTaus = reorder.taus ) 
      # replace by refined only of non-NA!
      if( is.finite(RSE.bestRef[["BestFit"]][1]) ) RSE.best <- RSE.bestRef
      if( is.finite(KDE.bestRef[["BestFit"]][1]) ) KDE.best <- KDE.bestRef
    }
    return( ..makeOutListForLMA( DIR.best, RSE.best, KDE.best, 3 ) )
  }
  return( environment() )
}
  
