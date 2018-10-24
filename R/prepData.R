
#' Preprocessing of the tracking data
#'
#' @param trackData A dataframe of the tracking data, including at least coordinates
#' (either longitude/latitude values or cartesian coordinates), and optionnaly a field \code{ID}
#' (identifiers for the observed individuals). Additionnal fields are considered as covariates.
#' Note that, if the names of the coordinates are not "x" and "y", the \code{coordNames} argument
#' should specified.
#' @param type \code{'LL'} if longitude/latitude provided (default), \code{'UTM'} if easting/northing.
#' @param coordNames Names of the columns of coordinates in the data frame. Default: \code{c("x","y")}.
#' @param LLangle Logical. If TRUE, the turning angle is calculated with \code{geosphere::bearing}
#' (default), else calculated with \code{atan2}.
#'
#' @return An object \code{moveData}, i.e. a dataframe of:
#' \item{ID}{The ID(s) of the observed animal(s)}
#' \item{step}{The step lengths - in kilometers if longitude/latitude provided, and in the metrics of
#' the data otherwise}
#' \item{angle}{The turning angles (if any) - in radians}
#' \item{x}{Either Easting or longitude (or e.g. depth for 1D data)}
#' \item{y}{Either Northing or latitude (all zero if 1D data)}
#' \item{...}{Covariates (if any)}
#'
#' @examples
#' coord1 <- c(1,2,3,4,5,6,7,8,9,10)
#' coord2 <- c(1,1,1,2,2,2,1,1,1,2)
#' trackData <- data.frame(coord1=coord1,coord2=coord2)
#' d <- prepData(trackData,type='UTM',coordNames=c("coord1","coord2"))
#'
#' @export
#' @importFrom geosphere distGeo distRhumb bearing bearingRhumb

prepData <- function(trackData, type=c('LL','UTM'), coordNames=c("x","y"), LLangle=NULL)
{
    # check arguments
    type <- match.arg(type)
    if(length(which(coordNames %in% names(trackData)))<2)
        stop("Check the columns names of your coordinates.")

    if(!is.null(trackData$ID))
        ID <- as.character(trackData$ID) # homogenization of numeric and string IDs
    else
        ID <- rep("Animal1",nrow(trackData)) # default ID if none provided

    if(length(which(is.na(ID)))>0)
        stop("Missing IDs")

    if(!is.null(LLangle) & !is.logical(LLangle))
        stop("'LLangle' must be logical.")
    if(is.null(LLangle))
        LLangle <- type=='LL' # TRUE if type=='LL', FALSE otherwise


    # remove tracks with less than two observations
    for(zoo in unique(ID)) {
        if(length(which(ID==zoo))<2) {
            trackData <- trackData[-which(ID==zoo),]
            ID <- ID[-which(ID==zoo)]
            warning(paste("Track",zoo,"only contains one observation,",
                          "and will be removed from the data."))
        }
    }

    x <- trackData[,coordNames[1]]
    y <- trackData[,coordNames[2]]

    # Pre-create data of the appropriate size to prevent slow rbinding
    nbPoints <- length(ID)
    data <- data.frame(ID=character(nbPoints),
                       step=numeric(nbPoints),
                       angle=numeric(nbPoints))
    levels(data$ID) <- unique(ID)

    nbAnimals <- length(unique(ID))

    # check that each animal's observations are contiguous
    for(i in 1:nbAnimals) {
        ind <- which(ID==unique(ID)[i])
        if(length(ind)!=length(ind[1]:ind[length(ind)]))
            stop("Each animal's obervations must be contiguous.")
    }

    eucDistance <- function(xyMatrix) {
      xyD <- diff(xyMatrix)
      sqrt(xyD[,1]^2 + xyD[,2]^2)
    }

    for(zoo in 1:nbAnimals) {
        nbObs <- length(which(ID==unique(ID)[zoo])) # number of observations for animal zoo
        i1 <- which(ID==unique(ID)[zoo])[1] # index of 1st obs for animal zoo
        i2 <- i1+nbObs-1 # index of last obs for animal zoo
        
        pt_matrix <- matrix(c(x[i1:i2], y[i1:i2]), ncol=2)
        if (type == 'LL') {
          step <- distGeo(pt_matrix)
        } else {
          step <- eucDistance(pt_matrix)
        }
        meters_per_km <- 1000
        step <- c(NA, step / meters_per_km)


        if (LLangle) {
          bear <- bearing(pt_matrix)
          angle <- c(NA, -diff(bear)/180*pi)
        } else {
          diffMatrix <- diff(pt_matrix)
          angle <- diff(atan(diffMatrix[,2]/diffMatrix[,1]))
          nonNAangles <- angle[!is.na(angle)]
          while(any(nonNAangles <= -pi) | any(nonNAangles > pi)) {
            angle <- ifelse(angle <= -pi, angle + 2*pi,
                       ifelse(angle > pi, angle - 2*pi, angle))
          }
          angle <- c(NA, angle, NA)
        }


        # d = data for one individual
        d <- data.frame(ID=rep(unique(ID)[zoo],nbObs),
                        step=step,
                        angle=angle)

        # append individual data to output
        data[i1:i2,] <- d
    }

    # Vectorized NA filler
    fillMissing <- function(covar) {
      N <- length(covar)
      naPos <- which(is.na(covar))
      
      # All NA or all non-NA
      if (length(naPos) %in% c(0, N)) return(covar)

      nonNaPos <- which(!is.na(covar))
      intervals  <- findInterval(naPos, nonNaPos, all.inside = TRUE)
      leftPos <- nonNaPos[pmax(1, intervals)]
      rightPos <- nonNaPos[pmin(N, intervals+1)]
      leftDist  <- naPos - leftPos
      rightDist <- rightPos - naPos

      covar[naPos] <- ifelse(leftDist <= rightDist, covar[leftPos], covar[rightPos])
      
      covar
    }
    # identify covariate columns
    covsCol <- which(names(trackData)!="ID" & names(trackData)!=coordNames[1] & names(trackData)!=coordNames[2])
    if(length(covsCol)>0) {
        covs <- data.frame(trackData[,covsCol]) # to prevent error if nbCovs==1
        colnames(covs) <- names(trackData)[covsCol]
        for (i in 1:length(covsCol)) { 
          covs[,i] <- fillMissing(covs[,i])
        }

    } else covs <- NULL

    data <- cbind(data,x=x,y=y)
    if(!is.null(covs))
        data <- cbind(data,covs)
    return(moveData(data))
}
