############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

setLMMultiCompFitEnvironment <- function( ParamsEnvr )
{
  # initialise routine vars (just to shorten names)
  .n.cores             <- as.integer( ParamsEnvr[["n.cores"]] )
  .x.tol               <- ParamsEnvr[["lma.params"]][["fit_control"]][["x_tol"]]
  .dt.rng              <- ParamsEnvr[["dt.hst.params"]][["dt_range_cut_ns"]]
  .bgr.ign.rng.ns      <- ParamsEnvr[["dt.hst.params"]][["bgr_range_to_ignore_ns"]]
  .log.cut.range       <- ParamsEnvr[["dt.hst.params"]][["log_cut_range_ns"]]
  .drop.thrshd.for.log <- ParamsEnvr[["dt.hst.params"]][["int_threshold_for_log_to_drop_au"]]
  .n.seeds             <- as.integer( ParamsEnvr[["n.guess"]] )
  .nls.control         <-  nls.lm.control( maxiter = ParamsEnvr[["lma.params"]][["fit_control"]][["n_iter_max"]], 
                                           ftol    = ParamsEnvr[["lma.params"]][["fit_control"]][["x_tol"]], 
                                           ptol    = ParamsEnvr[["lma.params"]][["fit_control"]][["x_tol"]],
                                           gtol    = ParamsEnvr[["lma.params"]][["fit_control"]][["x_tol"]] )
  # "public" vars to fit (names are important)
  p_tau.dir    <- ParamsEnvr[["dt.hst.params"]][["tau_dir_ns"]]     # for 3 components
  p_tau.pPs    <- ParamsEnvr[["dt.hst.params"]][["tau_pps_ns"]]     # for 3 components
  p_bgr.avg    <- NA     # average background
  p_Dt         <- NA
  p_hst.norm   <- NA
  p_DF.lin.fit <- NA
  p_DF.log.fit <- NA     # will be data frame
  
  # ----------------------------------------------------------------------------
  # helper functions
  .assessPIBgrFromTails <- function( DtHistVec,
                                     DtAxis,
                                     Median = TRUE ) # median or mean
  {
    if( Median ) 
      return( median(DtHistVec[DtAxis < .bgr.ign.rng.ns[1] | DtAxis > .bgr.ign.rng.ns[2]], na.rm = T ) )
    else 
      return( mean(DtHistVec[DtAxis < .bgr.ign.rng.ns[1] | DtAxis > .bgr.ign.rng.ns[2]], na.rm = T ) )
  }
  # confidence interval
  # TODO: consider making external
  ..confidenceInterval <- function(Vector, Interval, Median = TRUE)
  {
    # Standard deviation of sample
    vec_sd <- sd(Vector)
    # Sample size
    n <- length(Vector)
    # Mean of sample
    vec_mid <- ( if(Median) median(Vector) else mean(Vector) )
    # Error according to t distribution
    error <- qt((Interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
    # Confidence interval as a vector
    result <- c("lower" = vec_mid - error, "upper" = vec_mid + error)
    return(result)
  }
  ..makeVectorOfLowerBounds <- function( NVars, IniVecDf )
  {
    lower <- -Inf
    if( IniVecDf[["Bounded"]] ) lower <- -IniVecDf[["PosConstrErr"]]
    return( c(lower, -Inf, rep(lower, NVars - 2L)) )
  }
  ..shrinkHistByDtRange <- function( HstLength, LowThrldDt = NULL, HighThrldDt = NULL )
  {
    FACT.DT <- rep( TRUE, HstLength )  
    if( !is.null(LowThrldDt) )  FACT.DT <- FACT.DT & p_DF.lin.fit[["x"]] >= LowThrldDt
    if( !is.null(HighThrldDt) ) FACT.DT <- FACT.DT & p_DF.lin.fit[["x"]] <= HighThrldDt
    return( p_DF.lin.fit[FACT.DT,] )
  }

  # ----------------------------------------------------------------------------
  
  # spectrum pre-definition, both internal and external function
  defineSpectrumParams <- function( TAxNs, HstVec, HistParams )          # HistParams is external (not from ParamsEnvr[["dt.hst.params"]])
  {
    # probably low importance: no major improvement in fitting all-coincidence data
    if( HistParams[["reduce_bgr_for_randoms"]] )
    {
      hst.negative <-  HstVec[TAxNs < 0]
      bottom.CI    <- ..confidenceInterval( hst.negative, 0.9544, FALSE )[["lower"]] # use mean
      median.hst   <- median(hst.negative)
      min.hst      <- min(hst.negative)
      HstVec       <- HstVec - max( min.hst, min(bottom.CI, median.hst) )
      HstVec[HstVec < 0] <- 0
    }
    # TODO: make ImageTools2024.R with filters, resize, normalise etc
    # normalise to sum before cut
    HstVec <- HstVec / sum(HstVec)
    
    FACT.DT <- if( is.null(.dt.rng) ) rep(TRUE, length(TAxNs))
    else TAxNs >= .dt.rng[1] & TAxNs <= .dt.rng[2] 
    FAC.POS <- HstVec > .drop.thrshd.for.log
    FAC.LOG <- TAxNs >= .log.cut.range[1] & TAxNs <= .log.cut.range[2]
    
    p_Dt         <<- TAxNs[FACT.DT]
    p_hst.norm   <<- HstVec[FACT.DT]
    # data frames - for LMA fits
    p_DF.lin.fit <<- data.frame( x = p_Dt, y = p_hst.norm )
    p_DF.log.fit <<- data.frame( x = p_Dt[FAC.LOG & FAC.POS], y = log(p_hst.norm[FAC.LOG & FAC.POS]) )
    
    # IMPORTANT: assess background W I T H O U T filtering for log!
    p_bgr.avg <<- if(HistParams[["calc_bgr_from_full_span"]]) .assessPIBgrFromTails( HstVec, TAxNs, TRUE )
    else .assessPIBgrFromTails( p_hst.norm, p_Dt, TRUE )
  }
  # ------------------ FITTING MODELS ------------------------------------------
  # shrt == fast (various component notation for merging pPs + dir)
  ..makeLMFormula2CmpBgr <- function(Bg = p_bgr.avg)
  {
    as.formula(paste(
      "y ~ Ishrt * Lshrt * exp( Lshrt * (Mu + 0.5 * Lshrt * Sd^2 - x) ) * pnorm( (Mu - x) / Sd + Lshrt * Sd, lower.tail = F) +
           IoPs * LoPs * exp( LoPs * (Mu + 0.5 * LoPs * Sd^2 - x) ) * pnorm( (Mu - x) / Sd + LoPs * Sd, lower.tail = F) + ",
      Bg
    ))
  }
  ..makeLMFormula2CmpLogTauShrtBlurBgr <- function(Sd, Mu, TauShrt, Bg = p_bgr.avg)
  {
    Lshrt    <- 1 / TauShrt
    SdInv    <- 1 / Sd
    MuBySd   <- Mu / Sd
    as.formula(paste(
      "y ~ log( Ishrt * ", Lshrt, " * exp( ", Lshrt * Mu + 0.5 * (Lshrt * Sd)^2, " - ", Lshrt, " * x ) * pnorm(", MuBySd + Lshrt * Sd, " - x * ", SdInv, ", lower.tail = F) + 
                IoPs * LoPs * exp( LoPs * (", Mu, " + 0.5 * LoPs * ", Sd^2, " - x) ) * pnorm(", MuBySd, " + LoPs * ", Sd, " - x * ", SdInv, ", lower.tail = F) + ",
      Bg, ")"
    ))
  }
  ..makeLMFormula3CmpTausBgr <- function(Bg = p_bgr.avg)
  {
    Ldir <- 1/p_tau.dir
    LpPs <- 1/p_tau.pPs
    as.formula(paste(
      "y ~ Idir * ", Ldir, " * exp( ", Ldir, " * (Mu - x) + ", 0.5 * Ldir^2, " * Sd^2 ) * pnorm( (Mu - x) / Sd + ", Ldir, " * Sd, lower.tail = F) +
           IpPs * ", LpPs, " * exp( ", LpPs, " * (Mu - x) + ", 0.5 * LpPs^2, " * Sd^2 ) * pnorm( (Mu - x) / Sd + ", LpPs, " * Sd, lower.tail = F) +
           IoPs * LoPs * exp( LoPs * (Mu + 0.5 * LoPs * Sd^2 - x) ) * pnorm( (Mu - x) / Sd + LoPs * Sd, lower.tail = F) + ",
      Bg
    ))
  }
  ..makeLMFormula3CmpTausBgrFixedIratio <- function(Bg = p_bgr.avg)
  {
    Ldir <- 1/p_tau.dir
    LpPs <- 1/p_tau.pPs
    as.formula(paste(
      "y ~ Nrm * ( ", 0.6 * Ldir, " * exp( ", Ldir, " * (Mu - x) + ", 0.5 * Ldir^2, " * Sd^2 ) * pnorm( (Mu - x) / Sd + ", Ldir, " * Sd, lower.tail = F) +
                   ", 0.1 * LpPs, " * exp( ", LpPs, " * (Mu - x) + ", 0.5 * LpPs^2, " * Sd^2 ) * pnorm( (Mu - x) / Sd + ", LpPs, " * Sd, lower.tail = F) +
                   0.3 * LoPs * exp( LoPs * (Mu + 0.5 * LoPs * Sd^2 - x) ) * pnorm( (Mu - x) / Sd + LoPs * Sd, lower.tail = F) + ",
      Bg, " )"
    ))
  }
  ..makeLMFormula3CmpTaupPsBgr <- function(Bg = p_bgr.avg)
  {
    LpPs <- 1/p_tau.pPs
    as.formula(paste(
      "y ~ Idir * Ldir * exp( Ldir * (Mu + 0.5 * Ldir * Sd^2 - x) ) * pnorm( (Mu - x) / Sd + Ldir * Sd, lower.tail = F) +
           IpPs * ", LpPs, " * exp( ", LpPs, " * (Mu - x) + ", 0.5 * LpPs^2, " * Sd^2 ) * pnorm( (Mu - x) / Sd + ", LpPs, " * Sd, lower.tail = F) +
           IoPs * LoPs * exp( LoPs * (Mu + 0.5 * LoPs * Sd^2 - x) ) * pnorm( (Mu - x) / Sd + LoPs * Sd, lower.tail = F) + ",
      Bg
    ))
  }  
  ..makeLMFormula3CmpTaupPsBgrFixedIratio <- function(Bg = p_bgr.avg)
  {
    LpPs <- 1/p_tau.pPs
    as.formula(paste(
      "y ~  Nrm * ( 0.6 * Ldir * exp( Ldir * (Mu + 0.5 * Ldir * Sd^2 - x) ) * pnorm( (Mu - x) / Sd + Ldir * Sd, lower.tail = F) +
                    ", 0.1 * LpPs, " * exp( ", LpPs, " * (Mu - x) + ", 0.5 * LpPs^2, " * Sd^2 ) * pnorm( (Mu - x) / Sd + ", LpPs, " * Sd, lower.tail = F) +
                    0.3 * LoPs * exp( LoPs * (Mu + 0.5 * LoPs * Sd^2 - x) ) * pnorm( (Mu - x) / Sd + LoPs * Sd, lower.tail = F) + ",
      Bg, " )"
    ))
  }
  ..makeLMFormula3CmpTaupPsBgrFixedITausRatio <- function(Bg = p_bgr.avg)
  {
    Ldir <- 1/p_tau.dir
    LpPs <- 1/p_tau.pPs
    as.formula(paste(
      "y ~ Idir * Ldir * exp( Ldir * (Mu + 0.5 * Ldir * Sd^2 - x) ) * pnorm( (Mu - x) / Sd + Ldir * Sd, lower.tail = F) +
           IoPs / 3 * ", LpPs, " * exp( ", LpPs, " * (Mu - x) + ", 0.5 * LpPs^2, " * Sd^2 ) * pnorm( (Mu - x) / Sd + ", LpPs, " * Sd, lower.tail = F )  +
           IoPs * LoPs * exp( LoPs * (Mu + 0.5 * LoPs * Sd^2 - x) ) * pnorm( (Mu - x) / Sd + LoPs * Sd, lower.tail = F) + ",
      Bg
    ))
  }
  # --- log scale ---
  ..makeLMFormula3CmpLogTaupPsBlurBgr <- function(Sd, Mu, Bg = p_bgr.avg)
  {
    LpPs   <- 1 / p_tau.pPs
    SdInv  <- 1 / Sd
    MuBySd <- Mu / Sd
    as.formula(paste(
      "y ~ log( Idir * Ldir * exp( Ldir * (", Mu, " + 0.5 * Ldir * ", Sd^2, " - x) ) * pnorm(", MuBySd, " + Ldir * ", Sd, " - x * ", SdInv, ", lower.tail = F) + 
                IpPs * ", LpPs, " * exp( ", LpPs * Mu + 0.5 * (LpPs * Sd)^2, " - ", LpPs, " * x ) * pnorm(", MuBySd + LpPs * Sd, " - x * ", SdInv, ", lower.tail = F) + 
                IoPs * LoPs * exp( LoPs * (", Mu, " + 0.5 * LoPs * ", Sd^2, " - x) ) * pnorm(", MuBySd, " + LoPs * ", Sd, " - x * ", SdInv, ", lower.tail = F) + ",
      Bg, ")"
    ))
  }  
  ..makeLMFormula3CmpLogTaupPsBlurBgrFixedIratio <- function(Sd, Mu, Bg = p_bgr.avg)
  {
    LpPs   <- 1 / p_tau.pPs
    SdInv  <- 1 / Sd
    MuBySd <- Mu / Sd
    as.formula(paste(
      "y ~ log( Nrm * ( 0.6 * Ldir * exp( Ldir * (", Mu, " + 0.5 * Ldir * ", Sd^2, " - x) ) * pnorm(", MuBySd, " + Ldir * ", Sd, " - x * ", SdInv, ", lower.tail = F) + 
                       ", 0.1 * LpPs, " * exp( ", LpPs * Mu + 0.5 * (LpPs * Sd)^2, " - ", LpPs, " * x ) * pnorm(", MuBySd + LpPs * Sd, " - x * ", SdInv, ", lower.tail = F) + 
                       0.3 * LoPs * exp( LoPs * (", Mu, " + 0.5 * LoPs * ", Sd^2, " - x) ) * pnorm(", MuBySd, " + LoPs * ", Sd, " - x * ", SdInv, ", lower.tail = F) + ",
      Bg, ") )"
    ))
  }
  ..makeLMFormula3CmpLogTaupPsBlurBgrFixedITausRatio <- function(Sd, Mu, Bg = p_bgr.avg)
  {
    LpPs   <- 1 / p_tau.pPs
    SdInv  <- 1 / Sd
    MuBySd <- Mu / Sd
    as.formula(paste(
      "y ~ log( Idir * Ldir * exp( Ldir * (", Mu, " + 0.5 * Ldir * ", Sd^2, " - x) ) * pnorm(", MuBySd, " + Ldir * ", Sd, " - x * ", SdInv, ", lower.tail = F) + 
                IoPs *", LpPs / 3, " * exp( ", LpPs * Mu + 0.5 * (LpPs * Sd)^2, " - ", LpPs, " * x ) * pnorm(", MuBySd + LpPs * Sd, " - x * ", SdInv, ", lower.tail = F) + 
                IoPs * LoPs * exp( LoPs * (", Mu, " + 0.5 * LoPs * ", Sd^2, " - x) ) * pnorm(", MuBySd, " + LoPs * ", Sd, " - x * ", SdInv, ", lower.tail = F) + ",
      Bg, ")"
    ))
  }  
  # warning: no log!
  ..makeLMFormula3CmpTauoPsBgr <- function(TauoPs, Bg = p_bgr.avg)
  {
    LoPs   <- 1/TauoPs
    LpPs   <- 1/p_tau.pPs
    as.formula(paste(
      "y ~ Idir * Ldir * exp( Ldir * (Mu + 0.5 * Ldir * Sd^2 - x) ) * pnorm( (Mu - x) / Sd + Ldir * Sd, lower.tail = F) +
           IpPs * ", LpPs, " * exp( ", LpPs, " * (Mu - x) + ", 0.5 * LpPs^2, " * Sd^2 ) * pnorm( (Mu - x) / Sd + ", LpPs, " * Sd, lower.tail = F) +
           IoPs * ", LoPs, " * exp( ", LoPs, " * (Mu - x) + ", 0.5 * LoPs^2, " * Sd^2 ) * pnorm( (Mu - x) / Sd + ", LoPs, " * Sd, lower.tail = F) + ", 
      Bg
    ))
  }
  ..makeLMFormula3CmpTauoPsBgrFixedIratio <- function(TauoPs, Bg = p_bgr.avg)
  {
    LoPs   <- 1/TauoPs
    LpPs   <- 1/p_tau.pPs
    as.formula(paste(
      "y ~ Nrm * ( 0.6 * Ldir * exp( Ldir * (Mu + 0.5 * Ldir * Sd^2 - x) ) * pnorm( (Mu - x) / Sd + Ldir * Sd, lower.tail = F) +
                   ", 0.1 * LpPs, " * exp( ", LpPs, " * (Mu - x) + ", 0.5 * LpPs^2, " * Sd^2 ) * pnorm( (Mu - x) / Sd + ", LpPs, " * Sd, lower.tail = F) +
                   ", 0.3 * LoPs, " * exp( ", LoPs, " * (Mu - x) + ", 0.5 * LoPs^2, " * Sd^2 ) * pnorm( (Mu - x) / Sd + ", LoPs, " * Sd, lower.tail = F) + ", 
      Bg, " )"
    ))
  }
  # ----------------------------------------------------------------------------
  .oneCompEMG <- function( I, Tau, Sd, Mu, X )
  {
    Lambda <- 1 / Tau
    LBySd2 <- Lambda * Sd^2
    I * Lambda * exp( ( Mu - X + 0.5 * LBySd2 ) * Lambda ) * 
        pnorm( ( Mu - X + LBySd2 ) / Sd, lower.tail = F )
  }
  # all free
  # P = c(Sd, Mu, Ishrt, IoPs, tauShrt, taupPs, tauoPs, Bgr)
  .predict2CompGeneral <- function( P, X )
  {
    .oneCompEMG( P[3], P[5], P[1], P[2], X ) +   # direct + pPs
    .oneCompEMG( P[4], P[6], P[1], P[2], X ) + P[7]
  }
  # P = c(Sd, Mu, Idir, IpPs, IoPs, tauDir, taupPs, tauoPs, Bgr)
  .predict3CompGeneral <- function( P, X )
  {
    .oneCompEMG( P[3], P[6], P[1], P[2], X ) +   # direct
    .oneCompEMG( P[4], P[7], P[1], P[2], X ) +  # pPs
    .oneCompEMG( P[5], P[8], P[1], P[2], X ) + P[9]
  }
  # ----------------------------------------------------------------------------
  # WARNING: Params must include all taus: tau_dir, tau_pPs, tau_oPs (7 or 9 vars)
  .estimateRSELinear <- function( Params, 
                                  TAxNs   = NULL, 
                                  HstNorm = NULL )
  {
    n.vars  <- length(Params)
    dt      <- if( is.null(TAxNs) ) p_Dt else TAxNs
    hst.lin <- if( is.null(HstNorm) ) p_hst.norm else HstNorm
    predicted <- if(n.vars == 7) .predict2CompGeneral( Params, dt ) 
    else .predict3CompGeneral( Params, dt )
    Errs <- predicted - hst.lin
    IS.FINITE <- is.finite(Errs)
    return( sqrt( 
      sum( Errs[IS.FINITE]^2 / (sum(IS.FINITE) - n.vars - 1L) ) 
    ) )
  }
  .estimateRSELog <- function( Params, 
                               TAxNs      = NULL, 
                               HstNorm    = NULL )
  {
    n.vars  <- length(Params)
    dt      <- NA
    hst.log <- NA
    if( !is.null(HstNorm) & !is.null(TAxNs)  )
    {
      FAC.LOG <- HstNorm > .drop.thrshd.for.log     # above the threshold
      dt      <- TAxNs[ FAC.LOG ]
      hst.log <- log( HstNorm[ FAC.LOG ] )
    } else {
      dt      <- p_DF.log.fit[["x"]]
      hst.log <- p_DF.log.fit[["y"]]
    }
    predicted <- if(n.vars == 7) .predict2CompGeneral( Params, dt ) 
    else .predict3CompGeneral( Params, dt )
    Errs <- log( predicted ) - hst.log
    IS.FINITE <- is.finite(Errs)
    return( sqrt( 
      sum( Errs[IS.FINITE]^2 / (sum(IS.FINITE) - n.vars - 1L) ) 
    ) )
  }
  
  # --------- helpers for making multi-thread LMA tabs -------------------------
  
  ..makeTabLMA2AllFreeLin <- function( LMAResult, ReturnFullTable, AddErrors )
  {
    ncol.out <- if( AddErrors ) 7L else 4L # no full table
    if( ReturnFullTable ) 
      ncol.out <- if( AddErrors ) 14L else 8L # RSE + 7 + [6 errs]
    out.vec  <- rep( as.numeric(NA), ncol.out )
    
    if( class(LMAResult) == "nls" )
    {
      lm.out     <- summary(LMAResult)[["parameters"]][,1:2]
      lm.params  <- lm.out[,1]
      # convert lambdas -> taus and add background (p_bgr.avg must be pre-calculated)
      fit.params <- c( lm.params[1:4], 1/lm.params[5:6], p_bgr.avg ) 
      out.vec[1] <- .estimateRSELinear( fit.params, TAxNs   = p_DF.lin.fit[["x"]], 
                                        HstNorm = p_DF.lin.fit[["y"]] ) 
      
      if( ReturnFullTable ) out.vec[2:8] <- fit.params
      else                  out.vec[2:4] <- c( fit.params[1:2], min( fit.params[5:6] ) )
      if( AddErrors )
      {
        lm.errs <- summary(LMAResult)[["parameters"]][,2]
        if( ReturnFullTable )
          out.vec[9:14] <- c( lm.errs[1:4], 
                              abs( 1 / lm.params[5] - 1 / (lm.params[5] - lm.errs[5]) ),
                              abs( 1 / lm.params[6] - 1 / (lm.params[6] - lm.errs[6]) ) )
        else {
          tauS.id <- which.max( lm.params[5:6] ) + 4L # max lambda means min tau
          out.vec[5:7] <- c( lm.errs[1:2], abs( 1 / lm.params[tauS.id] - 1 / (lm.params[tauS.id] - lm.errs[tauS.id]) ) )
        }
      }
    } 
    return( out.vec ) 
  }
  
  ..makeTabLMA2Log <- function( LMAResult, VecGaussAndTauShort, AddErrors )
  {
    if( class(LMAResult) == "nls" ){
      lm.out       <- summary(LMAResult)[["parameters"]][,1:2]
      params.final <- c( VecGaussAndTauShort[1:2],
                         lm.out[1:2,1],
                         VecGaussAndTauShort[3],
                         1 / lm.out[3,1],
                         p_bgr.avg )
      if( AddErrors )
      {
        lm.errs <- lm.out[,2]
        params.final <- c( params.final, 
                           lm.errs[1:2], 
                           abs( 1 / lm.out[3,1] - 1 / ( lm.out[3,1] - lm.errs[3]) ) )
      }
      return( as.numeric(c( .estimateRSELog( params.final ), params.final )) )
    } else {
      if(AddErrors) return( rep( as.numeric(NA), 11)) else return( rep( as.numeric(NA), 8))
    }
  }
  
  # TODO: make DRY (too similar functions)
  ..makeTabLMA3TauDirLin <- function( LMAResult, ReturnFullTable, AddErrors )
  {
    ncol.out <- if( AddErrors ) 5L else 3L # no full table
    if( ReturnFullTable ) 
      ncol.out <- if( AddErrors ) 16L else 10L # no full table
    out.vec  <- rep( as.numeric(NA), ncol.out )
    if( class(LMAResult) == "nls" )
    {
      lm.out     <- summary(LMAResult)[["parameters"]][,1:2]
      lm.params  <- lm.out[,1]
      # convert lambda_oPs -> tau_oPs, add other taus and background (must be pre-calculated)
      fit.params <- c( lm.params[1:5], p_tau.dir, p_tau.pPs, 1 / lm.params[6], p_bgr.avg )
      out.vec[1] <- .estimateRSELinear( fit.params, TAxNs   = p_DF.lin.fit[["x"]], 
                                                    HstNorm = p_DF.lin.fit[["y"]] ) 
      if( ReturnFullTable ) 
        out.vec[2:10]  <- fit.params
      else out.vec[2:3] <- lm.params[1:2]
      if( AddErrors )
      {
        lm.errs <- lm.out[,2]
        if( ReturnFullTable ) 
          out.vec[11:16] <- c( lm.errs[1:5], abs( 1 / lm.params[6] - 1 / (lm.params[6] - lm.errs[6]) ) )
        else out.vec[4:5] <- lm.errs[1:2]
      }
    } 
    return( out.vec )
  }
  # with fixed intensity ratios
  ..makeTabLMA3TauDirLinIRatioFixed <- function( LMAResult, ReturnFullTable, AddErrors )
  {
    ncol.out <- if( AddErrors ) 5L else 3L # no full table
    if( ReturnFullTable ) 
      ncol.out <- if( AddErrors ) 14L else 10L # no full table
    out.vec  <- rep( as.numeric(NA), ncol.out )
    if( class(LMAResult) == "nls" )
    {
      lm.out     <- summary(LMAResult)[["parameters"]][,1:2]
      lm.params  <- lm.out[,1]
      # convert lambda_oPs -> tau_oPs, add other taus and background (must be pre-calculated)
      fit.params <- c( lm.params[1:2],  lm.params[3] * c(0.6, 0.1, 0.3),
                       p_tau.dir, p_tau.pPs, 1 / lm.params[4], p_bgr.avg )
      out.vec[1] <- .estimateRSELinear( fit.params, TAxNs   = p_DF.lin.fit[["x"]], 
                                                    HstNorm = p_DF.lin.fit[["y"]] ) 
      if( ReturnFullTable ) 
        out.vec[2:10]  <- fit.params
      else out.vec[2:3] <- lm.params[1:2]
      if( AddErrors ){
        lm.errs <- lm.out[,2]
        if( ReturnFullTable ) 
          out.vec[11:14] <- c( lm.errs[1:3], abs( 1 / lm.params[4] - 1 / (lm.params[4] - lm.errs[4]) ) )
        else out.vec[4:5] <- lm.errs[1:2]
      }
    } 
    return( out.vec )
  }

  ..makeTabLMA3AllFreeLin <- function( LMAResult, ReturnFullTable, AddErrors )
  {
    ncol.out <- if( AddErrors ) 5L else 3L # no full table
    if( ReturnFullTable ) 
      ncol.out <- if( AddErrors ) 17L else 10L
    out.vec  <- rep( as.numeric(NA), ncol.out )
    if( class(LMAResult) == "nls" )
    {
      lm.out     <- summary(LMAResult)[["parameters"]][,1:2]
      lm.params  <- lm.out[,1]
      # convert lambdas -> taus, add tau_pPs and background (must be pre-calculated)
      fit.params <- c( lm.params[1:5], 1 / lm.params[6], p_tau.pPs, 1 / lm.params[7], p_bgr.avg )
      out.vec[1] <- .estimateRSELinear( fit.params, TAxNs   = p_DF.lin.fit[["x"]], 
                                                    HstNorm = p_DF.lin.fit[["y"]] ) 
      if( ReturnFullTable ) 
        out.vec[2:10]  <- fit.params
      else out.vec[2:3] <- lm.params[1:2]
      if( AddErrors ){
        lm.errs   <- lm.out[,2]
        if( ReturnFullTable ) 
          out.vec[11:17] <- c( lm.errs[1:5], 
                               abs( 1 / lm.params[6] - 1 / (lm.params[6] - lm.errs[6]) ),
                               abs( 1 / lm.params[7] - 1 / (lm.params[7] - lm.errs[7]) ))
        else out.vec[4:5] <- lm.errs[1:2]
      }
    } 
    return( out.vec )
  }
  # with fixed intensity ratios
  ..makeTabLMA3AllFreeLinIRatioFixed <- function( LMAResult, ReturnFullTable, AddErrors )
  {
    ncol.out <- if( AddErrors ) 5L else 3L # no full table
    if( ReturnFullTable ) 
      ncol.out <- if( AddErrors ) 15L else 10L
    out.vec  <- rep( as.numeric(NA), ncol.out )
    if( class(LMAResult) == "nls" )
    {
      lm.out     <- summary(LMAResult)[["parameters"]][,1:2]
      lm.params  <- lm.out[,1]
      # convert lambdas -> taus, add tau_pPs and background (must be pre-calculated)
      fit.params <- c( lm.params[1:2],  lm.params[3] * c(0.6, 0.1, 0.3),
                       1 / lm.params[4], p_tau.pPs, 1 / lm.params[5], p_bgr.avg )
      out.vec[1] <- .estimateRSELinear( fit.params, TAxNs   = p_DF.lin.fit[["x"]], 
                                                    HstNorm = p_DF.lin.fit[["y"]] ) 
      if( ReturnFullTable ) 
        out.vec[2:10]  <- fit.params
      else out.vec[2:3] <- lm.params[1:2]
      if( AddErrors ){
        lm.errs   <- lm.out[,2]
        if( ReturnFullTable ) 
          out.vec[11:15] <- c( lm.errs[1:3], 
                               abs( 1 / lm.params[4] - 1 / (lm.params[4] - lm.errs[4]) ),
                               abs( 1 / lm.params[5] - 1 / (lm.params[5] - lm.errs[5]) ))
        else out.vec[4:5] <- lm.errs[1:2]
      }
    } 
    return( out.vec )
  }
  # with (partially) fixed intensity ratios for taus
  ..makeTabLMA3AllFreeLinITausRatioFixed <- function( LMAResult, ReturnFullTable, AddErrors )
  {
    ncol.out <- if( AddErrors ) 5L else 3L # no full table
    if( ReturnFullTable ) 
      ncol.out <- if( AddErrors ) 16L else 10L
    out.vec  <- rep( as.numeric(NA), ncol.out )
    if( class(LMAResult) == "nls" )
    {
      lm.out     <- summary(LMAResult)[["parameters"]][,1:2]
      lm.params  <- lm.out[,1]
      # convert lambdas -> taus, add tau_pPs and background (must be pre-calculated)
      fit.params <- c( lm.params[1:3], lm.params[4] / c(3, 1),
                       1 / lm.params[5], p_tau.pPs, 1 / lm.params[6], p_bgr.avg )
      out.vec[1] <- .estimateRSELinear( fit.params, TAxNs   = p_DF.lin.fit[["x"]], 
                                                    HstNorm = p_DF.lin.fit[["y"]] ) 
      if( ReturnFullTable ) 
        out.vec[2:10]   <- fit.params
      else out.vec[2:3] <- lm.params[1:2]
      if( AddErrors ){
        lm.errs   <- lm.out[,2]
        if( ReturnFullTable ) 
          out.vec[11:16] <- c( lm.errs[1:4], 
                               abs( 1 / lm.params[5] - 1 / (lm.params[5] - lm.errs[5]) ),
                               abs( 1 / lm.params[6] - 1 / (lm.params[6] - lm.errs[6]) ))
        else out.vec[4:5] <- lm.errs[1:2]
      }
    } 
    return( out.vec )
  }

  # ------------------ MAIN FUNCTIONS ------------------------------------------
  # WARNING: the spectrum must be predefined

  # --- 2-component ---
  # all free and return either all or only Gauss and tau_short
  runSeededLMA2AllFreeLinOMP <- function( INIGuessDF,                # ini.guess.df
                                          FullSpan        = FALSE,   # if true, ignores ini.guess.df[["UpThrldDt"]]
                                          LowThrldDt      = NULL,    # can be taken from dt.hst.params
                                          ReturnFullTable = FALSE,   
                                          AddErrors       = FALSE )  # if true, adds columns with errors
  {
    formula.low.fit <- ..makeLMFormula2CmpBgr()        # default Bg = p_bgr.avg (must be pre-calculated)
    hst.size        <- length( p_Dt )
    # main loop (OpenMP)
    Reduce("rbind", mclapply(seq_len(.n.seeds), function(i)
    {
      ini.vec.df <- INIGuessDF[i,] # data frame
      DF.local <- if( FullSpan ) p_DF.lin.fit 
                  else ..shrinkHistByDtRange( hst.size, LowThrldDt, ini.vec.df[["UpThrldDt"]] )
      lower.params <- ..makeVectorOfLowerBounds( 6L, ini.vec.df )
      result <- try({
        lm.fit.blur <- nlsLM( formula.low.fit,
                              data  = DF.local,
                              start = list( Sd    = ini.vec.df[["Sd"]],
                                            Mu    = ini.vec.df[["Mu"]],
                                            Ishrt = ini.vec.df[["Idir"]],
                                            IoPs  = ini.vec.df[["IoPs1"]],
                                            Lshrt = 1 / ini.vec.df[["Tau_dir"]],
                                            LoPs  = 1 / ini.vec.df[["Tau_oPs1"]] ),
                              trace   = FALSE,
                              lower   = lower.params,
                              control = .nls.control )
      }, silent = TRUE)
      return( ..makeTabLMA2AllFreeLin( result, ReturnFullTable, AddErrors ) )
    }, mc.cores = .n.cores))
  }
  # WARNING: the spectrum must be predefined
  runSeededLMA2LogOMP <- function( INIGuessDF,
                                   VecGaussAndTauShort,   # c(Sd, Mu, tau_short)
                                   AddErrors = FALSE )    # if true, adds columns with errors
  {
    formula.log.fit <- ..makeLMFormula2CmpLogTauShrtBlurBgr( Sd = VecGaussAndTauShort[1],
                                                             Mu = VecGaussAndTauShort[2],
                                                             TauShrt = VecGaussAndTauShort[3] ) # default Bg = bgr.avg (must be pre-calculated)
    Reduce("rbind", mclapply(seq_len(.n.seeds), function(i)
    {
      ini.vec.df <- INIGuessDF[i,] # data frame
      lower <- if( ini.vec.df[["Bounded"]] ) -ini.vec.df[["PosConstrErr"]] else -Inf
      result     <- try({
        lm.fit.log <- nlsLM( formula.log.fit,
                             data  = p_DF.log.fit,
                             start = list( Ishrt = ini.vec.df[["Idir"]],
                                           IoPs  = ini.vec.df[["IoPs1"]],
                                           LoPs  = 1 / ini.vec.df[["Tau_oPs1"]] ),
                             trace   = FALSE,
                             lower   = rep( lower, 3 ),
                             control = .nls.control )
      }, silent = TRUE)
      return( ..makeTabLMA2Log( result, VecGaussAndTauShort, AddErrors ) )
    }, mc.cores = .n.cores))
  }
  
  # --- 3-component ---
  
  # all free except tau_dir is fixed, return either all or only Gauss
  runSeededLMA3TauDirLinOMP <- function( INIGuessDF,
                                         FullSpan        = FALSE,     # if true, ignores INIS[["UpThrldDt"]]
                                         LowThrldDt      = NULL,
                                         ReturnFullTable = FALSE,
                                         AddErrors       = FALSE )    # if true, adds columns with errors
  {
    formula.low.fit <- ..makeLMFormula3CmpTausBgr()                             # default Bg = p_bgr.avg (must be pre-calculated)
    hst.size        <- length( p_Dt )
    Reduce("rbind", mclapply(seq_len(.n.seeds), function(i)
    {
      ini.vec.df <- INIGuessDF[i,] # data frame
      DF.local <- if( FullSpan ) p_DF.lin.fit 
                  else ..shrinkHistByDtRange( hst.size, LowThrldDt, ini.vec.df[["UpThrldDt"]] )
      lower.params <- ..makeVectorOfLowerBounds( 6L, ini.vec.df )
      result <- try({
        lm.fit.blur <- nlsLM( formula.low.fit,
                              data  = DF.local,
                              start = list( Sd    = ini.vec.df[["Sd"]],
                                            Mu    = ini.vec.df[["Mu"]],
                                            Idir  = ini.vec.df[["Idir"]],
                                            IpPs  = ini.vec.df[["IpPs"]],
                                            IoPs  = ini.vec.df[["IoPs1"]],
                                            LoPs  = 1 / ini.vec.df[["Tau_oPs1"]] ),
                              trace   = FALSE,
                              lower   = lower.params,
                              control = .nls.control )
      }, silent = TRUE)
      return( ..makeTabLMA3TauDirLin( result, ReturnFullTable, AddErrors ) )
    }, mc.cores = .n.cores))
  }
  # all free
  runSeededLMA3AllFreeLinOMP <- function( INIGuessDF,
                                          FullSpan        = FALSE,     # if true, ignores INIS[["UpThrldDt"]]
                                          LowThrldDt      = NULL,
                                          ReturnFullTable = FALSE,    # if false, returns only Sd and Mu
                                          AddErrors       = FALSE )
  {
    formula.low.fit <- ..makeLMFormula3CmpTaupPsBgr()                    # default Bg = p_bgr.avg (must be pre-calculated)
    hst.size        <- length( p_Dt )
    Reduce("rbind", mclapply(seq_len(.n.seeds), function(i)
    {
      ini.vec.df <- INIGuessDF[i,] # data frame
      DF.local <- if( FullSpan ) p_DF.lin.fit 
                  else ..shrinkHistByDtRange( hst.size, LowThrldDt, ini.vec.df[["UpThrldDt"]] )
      lower.params <- ..makeVectorOfLowerBounds( 7L, ini.vec.df )
      result <- try({
        lm.fit.blur <- nlsLM( formula.low.fit,
                              data  = DF.local,
                              start = list( Sd    = ini.vec.df[["Sd"]],
                                            Mu    = ini.vec.df[["Mu"]],
                                            Idir  = ini.vec.df[["Idir"]],
                                            IpPs  = ini.vec.df[["IpPs"]],
                                            IoPs  = ini.vec.df[["IoPs1"]],
                                            Ldir  = 1 / ini.vec.df[["Tau_dir"]],
                                            LoPs  = 1 / ini.vec.df[["Tau_oPs1"]] ),
                              trace   = FALSE,
                              lower   = lower.params,
                              control = .nls.control )
      }, silent = TRUE)
      return( ..makeTabLMA3AllFreeLin( result, ReturnFullTable, AddErrors ) )
    }, mc.cores = .n.cores))
  }
  
  runSeededLMA3FreeTauDirLogOMP <- function( INIGuessDF, BestGauss ) # BestGauss = c(Sd, Mu)
  {
    formula.log.fit <- ..makeLMFormula3CmpLogTaupPsBlurBgr( Sd = BestGauss[1],
                                                            Mu = BestGauss[2] ) # default Bg = bgr.avg (must be pre-calculated)
    Reduce("rbind", mclapply(seq_len(.n.seeds), function(i)
    {
      ini.vec.df <- INIGuessDF[i,] # data frame
      lower <- if( ini.vec.df[["Bounded"]] ) -ini.vec.df[["PosConstrErr"]] else -Inf
      result     <- try({
        lm.fit.log <- nlsLM( formula.log.fit,
                             data  = p_DF.log.fit,
                             start = list( Idir = ini.vec.df[["Idir"]],
                                           IpPs = ini.vec.df[["IpPs"]],
                                           IoPs = ini.vec.df[["IoPs1"]],
                                           Ldir = 1 / ini.vec.df[["Tau_dir"]],
                                           LoPs = 1 / ini.vec.df[["Tau_oPs1"]] ),
                             trace   = FALSE,
                             lower   = rep( lower, 5 ),
                             control = .nls.control )
      }, silent = TRUE)
      # no need for separate function
      if( class(result) == "nls" )
      {
        lm.params    <- summary(lm.fit.log)[["parameters"]]
        params.final <- c( BestGauss[1], BestGauss[2],
                           lm.params[1:3,1],
                           1 / lm.params[4,1], p_tau.pPs, 1 / lm.params[5,1],
                           p_bgr.avg )
        return( as.numeric(c( .estimateRSELog( params.final ), params.final )) )
      } else return( rep( as.numeric(NA), 10))
    }, mc.cores = .n.cores))
  }
  
  # (optional) refined third-stage fit with fixed tau_oPs
  # former function runCalcSeededLM3CompFixedTauoPsOMP()
  runSeededLMA3TauoPsLinOMP <- function( INIGuessDF, 
                                         TauoPs,
                                         DtRange = NULL )
  {
    formula.lin.fit <- ..makeLMFormula3CmpTauoPsBgr( TauoPs = TauoPs ) # default Bg = p_bgr.avg (must be pre-calculated)
    hst.size        <- length( p_Dt )
    # warning: different than for other stages, sets before OMP run
    DF.local   <- if( is.null(DtRange) ) p_DF.lin.fit 
                  else ..shrinkHistByDtRange( hst.size, DtRange[1], DtRange[2] )
    Reduce("rbind", mclapply(seq_len(.n.seeds), function(i)
    {
      ini.vec.df   <- INIGuessDF[i,] # data frame
      lower.params <- ..makeVectorOfLowerBounds( 6L, ini.vec.df )
      result <- try({
        lm.fit.oPs <- nlsLM( formula.lin.fit,
                             data  = DF.local,
                             start = list( Sd   = ini.vec.df[["Sd"]],
                                           Mu   = ini.vec.df[["Mu"]],
                                           Idir = ini.vec.df[["Idir"]],
                                           IpPs = ini.vec.df[["IpPs"]],
                                           IoPs = ini.vec.df[["IoPs1"]],
                                           Ldir = 1 / ini.vec.df[["Tau_dir"]] ),
                             trace   = FALSE,
                             lower   = lower.params,
                             control = .nls.control )
      }, silent = TRUE)
      if( class(result) == "nls" ){
        lm.params    <- summary(lm.fit.oPs)[["parameters"]][,1]
        params.final <- c( lm.params[1:5],
                           1 / lm.params[6],
                           p_tau.pPs, TauoPs,
                           p_bgr.avg )
        return( as.numeric(c( .estimateRSELinear( params.final, TAxNs   = p_DF.lin.fit[["x"]], 
                                                                HstNorm = p_DF.lin.fit[["y"]] ),
                              params.final )) )
      } else return( rep( as.numeric(NA), 10 ) )
    }, mc.cores = .n.cores))
  }
  
  # ----------------------------------------------------------------------------
  #      LMA3 with fixed intensity ratio I_dir : I_pPs : I_oPs = 6 : 1 : 3
  # ----------------------------------------------------------------------------
  # all free except tau_dir is fixed, return either all or only Gauss
  runSeededLMA3TauDirLinIRatioFixedOMP <- function( INIGuessDF,
                                                    FullSpan        = FALSE,     # if true, ignores INIS[["UpThrldDt"]]
                                                    LowThrldDt      = NULL,
                                                    ReturnFullTable = FALSE,
                                                    AddErrors       = FALSE )    # if true, adds columns with errors
  {
    formula.low.fit <- ..makeLMFormula3CmpTausBgrFixedIratio()                             # default Bg = p_bgr.avg (must be pre-calculated)
    hst.size        <- length( p_Dt )
    Reduce("rbind", mclapply(seq_len(.n.seeds), function(i)
    {
      ini.vec.df <- INIGuessDF[i,] # data frame
      DF.local <- if( FullSpan ) p_DF.lin.fit 
                  else ..shrinkHistByDtRange( hst.size, LowThrldDt, ini.vec.df[["UpThrldDt"]] )
      lower.params <- ..makeVectorOfLowerBounds( 4L, ini.vec.df )
      result <- try({
        lm.fit.blur <- nlsLM( formula.low.fit,
                              data  = DF.local,
                              start = list( Sd    = ini.vec.df[["Sd"]],
                                            Mu    = ini.vec.df[["Mu"]],
                                            Nrm   = ini.vec.df[["Idir"]],
                                            LoPs  = 1 / ini.vec.df[["Tau_oPs1"]] ),
                              trace   = FALSE,
                              lower   = lower.params,
                              control = .nls.control )
      }, silent = TRUE)
      return( ..makeTabLMA3TauDirLinIRatioFixed( result, ReturnFullTable, AddErrors ) )
    }, mc.cores = .n.cores))
  }
  
  # all free
  runSeededLMA3AllFreeLinIRatioFixedOMP <- function( INIGuessDF,
                                                     FullSpan        = FALSE,     # if true, ignores INIS[["UpThrldDt"]]
                                                     LowThrldDt      = NULL,
                                                     ReturnFullTable = FALSE,    # if false, returns only Sd and Mu
                                                     AddErrors       = FALSE )
  {
    formula.low.fit <- ..makeLMFormula3CmpTaupPsBgrFixedIratio()                    # default Bg = p_bgr.avg (must be pre-calculated)
    hst.size        <- length( p_Dt )
    Reduce("rbind", mclapply(seq_len(.n.seeds), function(i)
    {
      ini.vec.df <- INIGuessDF[i,] # data frame
      DF.local <- if( FullSpan ) p_DF.lin.fit 
                  else ..shrinkHistByDtRange( hst.size, LowThrldDt, ini.vec.df[["UpThrldDt"]] )
      lower.params <- ..makeVectorOfLowerBounds( 5L, ini.vec.df )
      result <- try({
        lm.fit.blur <- nlsLM( formula.low.fit,
                              data  = DF.local,
                              start = list( Sd    = ini.vec.df[["Sd"]],
                                            Mu    = ini.vec.df[["Mu"]],
                                            Nrm   = ini.vec.df[["Idir"]],
                                            Ldir  = 1 / ini.vec.df[["Tau_dir"]],
                                            LoPs  = 1 / ini.vec.df[["Tau_oPs1"]] ),
                              trace   = FALSE,
                              lower   = lower.params,
                              control = .nls.control )
      }, silent = TRUE)
      return( ..makeTabLMA3AllFreeLinIRatioFixed( result, ReturnFullTable, AddErrors ) )
    }, mc.cores = .n.cores))
  }
  
  runSeededLMA3FreeTauDirLogIRatioFixedOMP <- function( INIGuessDF, BestGauss ) # BestGauss = c(Sd, Mu)
  {
    formula.log.fit <- 
      ..makeLMFormula3CmpLogTaupPsBlurBgrFixedIratio( Sd = BestGauss[1],
                                                      Mu = BestGauss[2] ) # default Bg = bgr.avg (must be pre-calculated)
    Reduce("rbind", mclapply(seq_len(.n.seeds), function(i)
    {
      ini.vec.df <- INIGuessDF[i,] # data frame
      lower <- if( ini.vec.df[["Bounded"]] ) -ini.vec.df[["PosConstrErr"]] else -Inf
      result     <- try({
        lm.fit.log <- nlsLM( formula.log.fit,
                             data  = p_DF.log.fit,
                             start = list( Nrm = ini.vec.df[["Idir"]],
                                           Ldir = 1 / ini.vec.df[["Tau_dir"]],
                                           LoPs = 1 / ini.vec.df[["Tau_oPs1"]] ),
                             trace   = FALSE,
                             lower   = rep( lower, 3 ),
                             control = .nls.control )
      }, silent = TRUE)
      # no need for separate function
      if( class(result) == "nls" )
      {
        lm.params    <- summary(lm.fit.log)[["parameters"]]
        params.final <- c( BestGauss[1], BestGauss[2],
                           lm.params[1,1] * c(0.6, 0.1, 0.3),
                           1 / lm.params[2,1], p_tau.pPs, 1 / lm.params[3,1],
                           p_bgr.avg )
        return( as.numeric(c( .estimateRSELog( params.final ), params.final )) )
      } else return( rep( as.numeric(NA), 10))
    }, mc.cores = .n.cores))
  }
  
  # (optional) refined third-stage fit with fixed tau_oPs
  # former function runCalcSeededLM3CompFixedTauoPsOMP()
  runSeededLMA3TauoPsLinIRatioFixedOMP <- function( INIGuessDF, 
                                                    TauoPs,
                                                    DtRange = NULL )
  {
    formula.lin.fit <- ..makeLMFormula3CmpTauoPsBgrFixedIratio( TauoPs = TauoPs ) # default Bg = p_bgr.avg (must be pre-calculated)
    hst.size        <- length( p_Dt )
    # warning: different than for other stages, sets before OMP run
    DF.local   <- if( is.null(DtRange) ) p_DF.lin.fit 
                  else ..shrinkHistByDtRange( hst.size, DtRange[1], DtRange[2] )
    Reduce("rbind", mclapply(seq_len(.n.seeds), function(i)
    {
      ini.vec.df   <- INIGuessDF[i,] # data frame
      lower.params <- ..makeVectorOfLowerBounds( 4L, ini.vec.df )
      result <- try({
        lm.fit.oPs <- nlsLM( formula.lin.fit,
                             data  = DF.local,
                             start = list( Sd   = ini.vec.df[["Sd"]],
                                           Mu   = ini.vec.df[["Mu"]],
                                           Nrm  = ini.vec.df[["Idir"]],
                                           Ldir = 1 / ini.vec.df[["Tau_dir"]]  ),
                             trace   = FALSE,
                             lower   = lower.params,
                             control = .nls.control )
      }, silent = TRUE)
      if( class(result) == "nls" ){
        lm.params    <- summary(lm.fit.oPs)[["parameters"]][,1]
        params.final <- c( lm.params[1:2], lm.params[3] * c(0.6, 0.1, 0.3),
                           1 / lm.params[4],
                           p_tau.pPs, TauoPs,
                           p_bgr.avg )
        return( as.numeric(c( .estimateRSELinear( params.final, TAxNs   = p_DF.lin.fit[["x"]], 
                                                                HstNorm = p_DF.lin.fit[["y"]] ),
                              params.final )) )
      } else return( rep( as.numeric(NA), 10 ) )
    }, mc.cores = .n.cores))
  }
  
  # ----------------------------------------------------------------------------
  #      LMA3 with (partially) fixed intensity ratio I_pPs : I_oPs = 1 : 3
  # ----------------------------------------------------------------------------
  runSeededLMA3AllFreeLinITausRatioFixedOMP <- function( INIGuessDF,
                                                         FullSpan        = FALSE,     # if true, ignores INIS[["UpThrldDt"]]
                                                         LowThrldDt      = NULL,
                                                         ReturnFullTable = FALSE,    # if false, returns only Sd and Mu
                                                         AddErrors       = FALSE )
  {
    formula.low.fit <- ..makeLMFormula3CmpTaupPsBgrFixedITausRatio()            # default Bg = p_bgr.avg (must be pre-calculated)
    hst.size        <- length( p_Dt )
    Reduce("rbind", mclapply(seq_len(.n.seeds), function(i)
    {
      ini.vec.df <- INIGuessDF[i,] # data frame
      DF.local <- if( FullSpan ) p_DF.lin.fit 
      else ..shrinkHistByDtRange( hst.size, LowThrldDt, ini.vec.df[["UpThrldDt"]] )
      lower.params <- ..makeVectorOfLowerBounds( 6L, ini.vec.df )
      result <- try({
        lm.fit.blur <- nlsLM( formula.low.fit,
                              data  = DF.local,
                              start = list( Sd    = ini.vec.df[["Sd"]],
                                            Mu    = ini.vec.df[["Mu"]],
                                            Idir  = ini.vec.df[["Idir"]],
                                            IoPs  = ini.vec.df[["IoPs1"]],
                                            Ldir  = 1 / ini.vec.df[["Tau_dir"]],
                                            LoPs  = 1 / ini.vec.df[["Tau_oPs1"]] ),
                              trace   = FALSE,
                              lower   = lower.params,
                              control = .nls.control )
      }, silent = TRUE)
      return( ..makeTabLMA3AllFreeLinITausRatioFixed( result, ReturnFullTable, AddErrors ) )
    }, mc.cores = .n.cores))
  }
  runSeededLMA3FreeTauDirLogITausRatioFixedOMP <- function( INIGuessDF, BestGauss ) # BestGauss = c(Sd, Mu)
  {
    formula.log.fit <- 
      ..makeLMFormula3CmpLogTaupPsBlurBgrFixedITausRatio( Sd = BestGauss[1],
                                                          Mu = BestGauss[2] ) # default Bg = bgr.avg (must be pre-calculated)
    Reduce("rbind", mclapply(seq_len(.n.seeds), function(i)
    {
      ini.vec.df <- INIGuessDF[i,] # data frame
      lower <- if( ini.vec.df[["Bounded"]] ) -ini.vec.df[["PosConstrErr"]] else -Inf
      result     <- try({
        lm.fit.log <- nlsLM( formula.log.fit,
                             data  = p_DF.log.fit,
                             start = list( Idir  = ini.vec.df[["Idir"]],
                                           IoPs  = ini.vec.df[["IoPs1"]],
                                           Ldir  = 1 / ini.vec.df[["Tau_dir"]],
                                           LoPs  = 1 / ini.vec.df[["Tau_oPs1"]] ),
                             trace   = FALSE,
                             lower   = rep( lower, 4 ),
                             control = .nls.control )
      }, silent = TRUE)
      # no need for separate function
      if( class(result) == "nls" )
      {
        lm.params    <- summary(lm.fit.log)[["parameters"]]
        params.final <- c( BestGauss[1], BestGauss[2],
                           lm.params[1], lm.params[2] / c(3, 1),
                           1 / lm.params[3,1], p_tau.pPs, 1 / lm.params[4,1],
                           p_bgr.avg )
        return( as.numeric(c( .estimateRSELog( params.final ), params.final )) )
      } else return( rep( as.numeric(NA), 10))
    }, mc.cores = .n.cores))
  }  

  rm( list = c("ParamsEnvr") ) 
  return( environment() )
}