############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################
# WARNING: requires imageTools.R

require("jsonlite", quietly = TRUE)
require("minpack.lm", quietly = TRUE)
require("parallel", quietly = TRUE)
require("ks", quietly = TRUE)

current.dir <- getwd()
# extract directory of a script
args <- commandArgs(trailingOnly = FALSE)
# solution from here: stackoverflow.com/a/32016824/538403
script.dir  <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
# again, but no system info
args <- commandArgs(trailingOnly = TRUE)
# HELP & usage
USAGE.LINE <- paste("\nUsage:\n",
                    "Rscript [--vanilla] runMultiVoxelFitLMA2.R", "-p <json_params>\n", 
                    "-dt <min,max,n_bins>", "-s <lifetime_spectra>\n",
                    "-vox-ids <voxel_coordinates>", "-o <output_prefix>\n",
                    "[-vox-size <x_mm,y_mm,z_mm>]", "[-preserve-na]", "[-na-to-zero]\n",
                    "[-v] [> <output_log>]\n\n")
if( length(args) == 0 | args[1] == "-h" | args[1] == "--help" | args[1] == "-?")
{
  cat("\n--- HELP ---\n")
  cat( USAGE.LINE )
  cat("-p/--params : path to JSON file with combined parameters\n")
  cat("-s/--hists : path to RDS or a 4-byte binary with histograms of time delay spectra (as rows in a tab)\n")
  cat("-dt : time delay scale as [min,max,no of bins] (in nanoseconds)\n")
  cat("-vox-ids : path to RDS or a 2-byte binary with voxel IDs (a tab of n_histograms x 3)\n")
  cat("-o/--out : output prefix (for JSON and RDS)\n")
  cat("\nOptional:\n")
  cat("-vox-size : (only for info) size of a voxel [x,y,z] (in milimetres) \n")
  cat("-preserve-na : voxels with non-converged fits are coerced to NAs (default - to median)\n")
  cat("-na-to-zero : voxels with non-converged fits are coerced to 0 (default - to median)\n")
  cat("-v/--verbose : shows info during each LMA fitting\n\n")
  stop("Stopped.")
}
params.json   <- NULL
output.prefix <- NULL
hists.path    <- NULL
voxels.path   <- NULL
dt.params     <- NULL
vox.size      <- NULL     
coerce.na.to  <- "median" # default: NAs are replaced by a median from other voxels
verbose       <- FALSE

for( i in seq_along(args) )
{
  if(args[i] == "-p" | args[i] == "--params") params.json <- args[i + 1]
  else if(args[i] == "-s" | args[i] == "--hists") hists.path <- args[i + 1]
  else if(args[i] == "-vox-ids" ) voxels.path <- args[i + 1]
  else if(args[i] == "-dt") dt.params <- as.numeric(strsplit(args[i + 1], ",")[[1]])
  else if(args[i] == "-o" | args[i] == "--out") output.prefix <- args[i + 1]
  else if(args[i] == "-vox-size") vox.size <- as.numeric(strsplit(args[i + 1], ",")[[1]])
  else if(args[i] == "-preserve-na") coerce.na.to <- NULL
  else if(args[i] == "-na-to-zero") coerce.na.to <- "zero"
  else if(args[i] == "-v" | args[i] == "--verbose") verbose <- TRUE
}

if( is.null(params.json) | is.null(output.prefix) | is.null(voxels.path) |
    is.null(hists.path) | length(dt.params) != 3 )
{
  cat( USAGE.LINE )
  stop("Improper input!")
}

##############################################
cat("\nCurrent directory:", current.dir, "\n")
cat("\nPath to JSON with parameters:", params.json, "\n")
cat("Path to the histograms:", hists.path, "\n")
cat("Path to the voxel IDs:", voxels.path, "\n")
cat("Output prefix:", output.prefix, "\n\n")

cat("Processing source files... ")
source( file.path(script.dir, "miscTools.R") ) # load additional functions
source( file.path(script.dir, "setParamsEnvr.R") )
Rcpp::sourceCpp( file.path(script.dir, "../cpp/kdeCPP.cpp") )
source( file.path(script.dir, "setLMAEnvr.R") )
source( file.path(script.dir, "setBestFitEnvr.R") )
source( file.path(script.dir, "setPsLifetimeLMAFitEnvr.R") )

dt.ns  <- seq(dt.params[1], dt.params[2], length.out = dt.params[3])
n.bins <- length(dt.ns)

# --- read histograms and voxels ---
read.sp.as.RDS <- fileExtentionLowercase(hists.path) == ".rds"
read.vx.as.RDS <- fileExtentionLowercase(voxels.path) == ".rds"
hists.tab <- 
  if(read.sp.as.RDS) readRDS(hists.path) else importTableFromBinaryFile(hists.path, 
                                                                        n.bins, "double", 4L)
n.hsts <- nrow( hists.tab )
vox.ids.tab <- 
  if(read.vx.as.RDS) readRDS(voxels.path) else importTableFromBinaryFile(voxels.path, 
                                                                         3, "integer", 2L)
if( n.hsts != nrow( vox.ids.tab ) ) stop("Incopatible data for histograms and voxel IDs!")

cat("Done!\nLaunching multi-voxel minimisation...\n")
# set fit environment
fit.env <- setPsLifetimeLMAFitEnvironment()
fit.env[["configureSetupFromJSON"]]( params.json )

# TODO: (WARNING: No RSE, R2_adj, Xi2_red yet)
vec.NA.template <- rep( as.numeric(NA), n.hsts )
df.template     <- data.frame( Sd       = vec.NA.template,
                               Mu       = vec.NA.template,
                               I_shrt   = vec.NA.template,
                               I_oPs    = vec.NA.template,
                               tau_shrt = vec.NA.template,
                               tau_oPs  = vec.NA.template,
                               Bg       = vec.NA.template,
                               RSElin   = vec.NA.template,
                               RSElog   = vec.NA.template )
out <- list( byKDE       = df.template,
             byRSE       = df.template,
             byDirectFit = df.template ) 
rm(df.template) # clean memory

START.TIME <- Sys.time()
# main multi-voxel loop
for( i in seq_len(n.hsts) )
{
  LMA2fit <- fit.env[["fitHistByTwoStageLMA2OMP"]]( dt.ns, hists.tab[i,], Verbose = verbose )
  # save results to list
  for( nm in names(out) ) if( !is.na(LMA2fit[[nm]][["par"]][1]) )
  {
    out[[nm]][["Sd"]][i]       <- as.numeric( LMA2fit[[nm]][["par"]][1] )
    out[[nm]][["Mu"]][i]       <- as.numeric( LMA2fit[[nm]][["par"]][2] )
    out[[nm]][["I_shrt"]][i]   <- as.numeric( LMA2fit[[nm]][["par"]][3] )
    out[[nm]][["I_oPs"]][i]    <- as.numeric( LMA2fit[[nm]][["par"]][4] )
    out[[nm]][["tau_shrt"]][i] <- as.numeric( LMA2fit[[nm]][["par"]][5] )
    out[[nm]][["tau_oPs"]][i]  <- as.numeric( LMA2fit[[nm]][["par"]][6] )
    out[[nm]][["Bg"]][i]       <- as.numeric( LMA2fit[[nm]][["par"]][7] )
    out[[nm]][["RSElin"]][i]   <- as.numeric( LMA2fit[[nm]][["RSElin"]] )
    out[[nm]][["RSElog"]][i]   <- as.numeric( LMA2fit[[nm]][["RSElog"]] )
  }
  # show on screen
  cat("\n\t====== Voxel no.", i, "=====\n\tTau_oPs, ns (rse, kde), I_shrt / I_oPs (rse, kde), Bgr: | ", 
      sprintf( "%3.3f", out[["byRSE"]][["tau_oPs"]][i] ), " ",
      sprintf( "%3.3f", out[["byKDE"]][["tau_oPs"]][i] ), " | ",
      sprintf( "%4.3f", out[["byRSE"]][["I_shrt"]][i] / out[["byRSE"]][["I_oPs"]][i] ), " ",
      sprintf( "%3.3f", out[["byKDE"]][["I_shrt"]][i] / out[["byKDE"]][["I_oPs"]][i] ), " | ",
      sprintf( "%3.4g", out[["byRSE"]][["Bg"]][i] ), "\n\n" )
  if(i %% 10 == 0) 
  {
    ELAPSED.TIME <- difftime(Sys.time(), START.TIME, units = "secs")
    cat(sprintf("( Elapsed time for %d voxels: %02d:%02d:%02d )\n\n", i, 
                as.integer(ELAPSED.TIME) %/% 3600, 
                (as.integer(ELAPSED.TIME) %% 3600) %/% 60, 
                as.integer(ELAPSED.TIME) %% 60))
  }
}

# --- export results --
cat("\nResults written to files:\n")
rds.path  <- paste0(output.prefix,".rds")
saveRDS( out, file = rds.path )
cat(rds.path, "\n" )

for( nm in names(out) )
{
  out.path <- paste0(output.prefix,"_",nm,".dat")
  best.img <- makeArrayFromValsForVoxels( out[[nm]][["tau_oPs"]], vox.ids.tab, 
                                          ReplaceNAsBy = coerce.na.to )
  writeBin(c(best.img), out.path, size = 4L, endian = "little")
  cat(out.path, "\n" )
}

# best image
# TODO: consider making separate function if too similar with LMA3 
best.val.unr <- vec.NA.template
for( i in seq_len(n.hsts) )
{
  min.id <- which.min( c(out[["byRSE"]][["RSElog"]][i], 
                         out[["byKDE"]][["RSElog"]][i], 
                         out[["byDirectFit"]][["RSElog"]][i]) ) 
  # which.min() ignores NAs or returns integer(0) if all are NA
  if(length(min.id) == 0) best.val.unr[i] <- NA
  else best.val.unr[i] <- c(out[["byRSE"]][["tau_oPs"]][i], 
                            out[["byKDE"]][["tau_oPs"]][i], 
                            out[["byDirectFit"]][["tau_oPs"]][i])[min.id]
}
best.img <- makeArrayFromValsForVoxels( best.val.unr, vox.ids.tab, 
                                        ReplaceNAsBy = coerce.na.to )
out.path <- paste0(output.prefix,"_bestRSElog.dat")
writeBin(c(best.img), out.path, size = 4L, endian = "little")
cat(out.path, "\n" )

# save "interfile" with info 
interfile.name <- paste0(output.prefix,".hdr")
saveMapInfoToInterfile(interfile.name, vox.size, dim(best.img))
cat(interfile.name, "\n" )
