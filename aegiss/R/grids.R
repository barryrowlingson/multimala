##' Grid for original Aegiss data
##'
##' x, y, size for aegiss grid
##' @title Hampshire Grid Data
##' @return a list of xgrid, ygrid, size, using OSGB coords
##' @author Barry S Rowlingson
##' @export
hampshire <- function(){
    xgrid <- seq(402000+(89000/128/2),491000-(89000/128/2),l=128)
    ygrid <- seq(87000+(89000/128/2),176000-(89000/128/2),l=128)
    size <- 128
    return(list(x=xgrid, y=ygrid, size=size))
}
