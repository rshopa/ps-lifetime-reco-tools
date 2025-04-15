############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

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
                    "Rscript [--vanilla] fitOneHistogramLMA3.R", "-p <json_params>", 
                    "-s <lifetime_spectrum_ascii>", "-o <output_prefix>\n",
                    "[-fixed-ints-ratio] [-three-stage] [-d <material_density>]",
                    "[-v] [> <output_log>]\n\n")
if( length(args) == 0 | args[1] == "-h" | args[1] == "--help" | args[1] == "-?")
{
  cat("\n--- HELP ---\n")
  cat( USAGE.LINE )
  cat("-p/--params : path to JSON file with combined parameters\n")
  cat("-s/--hist : path to ASCII with a histogram of time delay spectrum (2-column [dt, ns | hst, au])\n")
  cat("-o/--out : output prefix (for JSON and RDS)\n")
  cat("\nOptional:\n")
  cat("-fixed-ints-ratio : fixes intensity ratio - direct : p-Ps : o-Ps = 60% : 10% : 30%\n")
  cat("-three-stage : adds a third-stage fit, linear scale with fixed tau_oPs\n")
  cat("-d/--density : in [g/mL], adjusts tau_dir by material density\n")
  cat("-v/--verbose : shows info during LMA fitting\n\n")
  stop("Stopped.")
}
params.json   <- NULL
output.prefix <- NULL
hist.path     <- NULL
# optional
int.ratio.fixed <- FALSE
three.stage     <- FALSE
matter.density  <- NULL
verbose         <- FALSE
# TODO: make as an external function
# Loop through arguments
for( i in seq_along(args) )
{
  if(args[i] == "-p" | args[i] == "--params") params.json <- args[i + 1]
  else if(args[i] == "-s" | args[i] == "--hist") hist.path <- args[i + 1]
  else if(args[i] == "-o" | args[i] == "--out") output.prefix <- args[i + 1]
  else if(args[i] == "-fixed-ints-ratio") int.ratio.fixed <- TRUE
  else if(args[i] == "-three-stage") three.stage <- TRUE
  else if(args[i] == "-d" | args[i] == "--density") 
    matter.density <- as.numeric(args[i + 1])
  else if(args[i] == "-v" | args[i] == "--verbose") verbose <- TRUE
}
if( is.null(params.json) | is.null(output.prefix) | is.null(hist.path) ) 
{
  cat( USAGE.LINE )
  stop("Improper input!")
}

cat("\nCurrent directory:", current.dir, "\n")
cat("\nPath to JSON with parameters:", params.json, "\n")
cat("Path to the ASCII histogram:", hist.path, "\n")
cat("Output prefix:", output.prefix, "\n\n")

cat("Processing source files... ")
source( file.path(script.dir,"setParamsEnvr.R") )
Rcpp::sourceCpp( file.path(script.dir,"../cpp/kdeCPP.cpp") )
source( file.path(script.dir,"setLMAEnvr.R") )
source( file.path(script.dir,"setBestFitEnvr.R") )
source( file.path(script.dir,"setPsLifetimeLMAFitEnvr.R") )

cat("Done!\nLaunching minimisation...\n")
# set fit environment
fit.env <- setPsLifetimeLMAFitEnvironment()
fit.env[["configureSetupFromJSON"]]( params.json )
# update range for the third stage (if needed)
dt.range.3rd.stage <- if(three.stage) 
  fit.env[[".pars.env"]][["dt.hst.params"]][["third_stage_dt_range_cut_ns"]] else NULL

# load spectrum data (data frame)
DF.hist <- read.table( hist.path, header = FALSE, col.names = c("Dt", "Hst") )
# launch minimisation
LMA3fit <- 
  fit.env[["fitHistByMultiStageLMA3IratioSwitchableOMP"]]( DF.hist[["Dt"]], 
                                                           DF.hist[["Hst"]], 
                                                           IRatioFixed = int.ratio.fixed,
                                                           RefineWFixedTauoPs = three.stage,
                                                           MatterDensitygmL = matter.density,
                                                           RefinedDtRange = dt.range.3rd.stage,
                                                           Verbose = verbose )
cat("Done!\n")
if(verbose)
{
  cat("\nOutput (LMA3):\n")
  print(str(LMA3fit, digits.d = 4, vec.len = 10))
}
# save output as JSON and RDS
json.path <- paste0(output.prefix,".json")
rds.path  <- paste0(output.prefix,".rds")
write_json( LMA3fit, pretty = TRUE, auto_unbox = T, digits = 6, path = json.path )
saveRDS( LMA3fit, file = rds.path )
cat("\nResults written to files:\n", json.path, "\n", rds.path, "\n" )
