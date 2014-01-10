##' Read multimala output files
##'
##' multimala generates simple streamed output files. This loads them
##' into raster objects
##' @title Read Multimala outputs
##' @param filename one of the output files
##' @param grid grid specification
##' @param mask optional mask
##' @return a raster of the values from the file
##' @author Barry S Rowlingson
##' @export
read_multimala_output <- function(filename, grid, mask){
    size = grid$size
    data <- matrix(
        scan(filename),ncol=2*size)
    data = t(data[1:size,size:1])
    r = raster::raster(data,
                    xmn=min(grid$x),
                    xmx=max(grid$x),
                    ymn=min(grid$y),
                    ymx=max(grid$y)
                    )    
    projection(r)="+init=epsg:27700"

    if(!missing(mask)){
        r[!mask] <- NA
    }
    
    return(r)
}
