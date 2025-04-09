############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

# --- helpers for processing images or binaries ----
# extracts file extension as a tail + dot, eg ".txt"
fileExtentionLowercase <- function(Path, 
                                   ExtLen = 3)
{
  len <- nchar(Path)
  tolower(substr(Path, len-ExtLen, len))
}
# normalises vector or array
normalize <- function( X, 
                       Type = "peak" )
{
  if(Type == "peak") return( (X - min(X) ) / (max(X) - min(X)) ) 
  else if(Type == "zero") return( X / max(abs(X)) ) else return( X / sum(X) )
}
# imports vector from binary file
importBinaryFromFile <- function( F.Path, 
                                  Var.Type = "double", 
                                  Val.Size = 8L, 
                                  Signed   = TRUE )
{
  n.elements <- file.info(F.Path)[["size"]] / Val.Size
  # import from file
  to.read <- file(F.Path, "rb")
  out <- readBin( to.read, 
                  Var.Type, n.elements, 
                  Val.Size, endian = "little", signed = Signed )
  close(to.read)
  return(out)
}
# imports table with given no of columns and unknown no of rows
importTableFromBinaryFile <- function( F.Path, 
                                       N.Cols,
                                       Var.Type = "double", 
                                       Val.Size = 8L, 
                                       Signed   = TRUE )
{
  unrolled <- importBinaryFromFile( F.Path, Var.Type, Val.Size, Signed )
  return( array(unrolled, c(length(unrolled) / N.Cols, N.Cols)) ) 
}

# (for axes): extracts axis step, e.g. voxel size
getAxisStep <- function(Axis)
{
  median(diff(sort(unique(Axis))))
}
# sorts unique axis values and extends by adding Edge to head and tail
flattenAndExtendAxis <- function(Axis, Edge = 0)
{
  ax.span <- max(abs(Axis))*c(-1,1)
  d.ax    <- getAxisStep(Axis)
  return( seq(ax.span[1] - Edge*d.ax, ax.span[2] + Edge*d.ax, by = d.ax) )
}
# makes image of vectors with values for specific voxels
makeArrayFromValsForVoxels <- function(Vals, 
                                       VoxIDs, 
                                       AddEdge = 2, 
                                       ReplaceNAsBy = NULL,
                                       Verbose = FALSE)
{
  axes    <- apply(VoxIDs, 2, flattenAndExtendAxis, Edge = AddEdge)
  map.out <- array(0, sapply(axes, length))
  for( k in 1:nrow(VoxIDs) )
  {
    if(Verbose) cat(paste("Voxel no:", k, "\r"))
    map.out[which(VoxIDs[k,1] == axes[[1]]),
            which(VoxIDs[k,2] == axes[[2]]),
            which(VoxIDs[k,3] == axes[[3]])] <- Vals[k]
  }
  if( !is.null(ReplaceNAsBy) )
  {
    if(ReplaceNAsBy == "median")
      map.out[is.na(map.out)] <- median( map.out[!is.na(map.out) & map.out >0] )
    else map.out[is.na(map.out)] <- 0
  }
  return( map.out )
}
# saves info in INTERFILE format 
saveMapInfoToInterfile <- function(FilePath, VoxSize, MapSize)
{
  # Open the file for writing
  FileConn <- file(FilePath, open = "wt")
  writeLines("!INTERFILE :=\n", FileConn)
  writeLines("!GENERAL DATA :=", FileConn)
  writeLines("imagedata byte order := LITTLEENDIAN\n", FileConn)
  writeLines("!total number of images := 1", FileConn)
  writeLines("number of dimensions := 3", FileConn)
  writeLines(paste("!matrix size [1] :=", MapSize[1]), FileConn)
  writeLines(paste("!matrix size [2] :=", MapSize[2]), FileConn)
  writeLines(paste("!matrix size [3] :=", MapSize[3]), FileConn)
  writeLines("!number format := float", FileConn)
  writeLines("!number of bytes per pixel := 4", FileConn)
  vox.scale <- if(is.null(VoxSize)) rep("unknown", 3) else VoxSize
  writeLines(paste("scaling factor (mm/pixel) [1] :=", vox.scale[1]), FileConn)
  writeLines(paste("scaling factor (mm/pixel) [2] :=", vox.scale[2]), FileConn)
  writeLines(paste("scaling factor (mm/pixel) [3] :=", vox.scale[3]), FileConn)
  close(FileConn)
}
