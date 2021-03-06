
#' Plot \code{moveData}
#' @method plot moveData
#'
#' @param x An object \code{moveData}
#' @param animals Vector of indices or IDs of animals for which information will be plotted.
#' Default: \code{NULL} ; all animals are plotted.
#' @param compact \code{TRUE} for a compact plot (all individuals at once), \code{FALSE} otherwise
#' (default -- one individual at a time).
#' @param ask If \code{TRUE}, the execution pauses between each plot.
#' @param breaks Histogram parameter. See \code{hist} documentation.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' # data is a moveData object (as returned by prepData), automatically loaded with the package
#' data <- example$data
#'
#' plot(data,compact=TRUE,breaks=20,ask=FALSE)
#'
#' @export
#' @importFrom graphics abline axis hist mtext par plot points
#' @importFrom grDevices hcl

plot.moveData <- function(x,animals=NULL,compact=FALSE,ask=TRUE,breaks="Sturges",...)
{
    data <- x
    # check arguments
    if(length(data)<1) stop("The data input is empty.")
    if(is.null(data$ID) | is.null(data$x) | is.null(data$step))
        stop("Missing field(s) in data.")

    nbAnimals <- length(unique(data$ID))

    if(all(data$y==0) & compact) {
        warning("One-dimensional data cannot plotted with the option 'compact'.")
        compact <- FALSE
    }

    ##################################
    ## Define animals to be plotted ##
    ##################################
    if(is.null(animals)) # all animals are plotted
        animalsInd <- 1:nbAnimals
    else {
        if(is.character(animals)) { # animals' IDs provided
            animalsInd <- NULL
            for(zoo in 1:length(animals)) {
                if(length(which(unique(data$ID)==animals[zoo]))==0) # ID not found
                    stop("Check animals argument.")

                animalsInd <- c(animalsInd,which(unique(data$ID)==animals[zoo]))
            }
        }

        if(is.numeric(animals)) { # animals' indices provided
            if(length(which(animals<1))>0 | length(which(animals>nbAnimals))>0) # index out of bounds
                stop("Check animals argument.")

            animalsInd <- animals
        }
    }

    # graphical parameters
    par(mar=c(5,4,4,2)-c(0,0,2,1)) # bottom, left, top, right
    par(ask=ask)

    if(is.null(data$angle) | all(data$y==0)) # only step length is provided
    {
        for(zoo in animalsInd) {
            ################
            ## Plot track ##
            ################*
            ID <- unique(data$ID)[zoo]
            x <- data$x[which(data$ID==ID)]
            plot(x, type="o",lwd=1.3, xlab="time", ylab="x", pch=20)

            ##########################################
            ## Plot steps time series and histogram ##
            ##########################################
            par(mfrow=c(1,2))
            step <- data$step[which(data$ID==ID)]
            # step length time series
            plot(step,type="l",xlab="t",ylab="step length",
                 ylim=c(0,max(step,na.rm=T)))
            # step length histogram
            hist(step,xlab="step length",main="",col="grey",border="white",breaks=breaks)
            mtext(paste("Animal ID:",ID),side=3,outer=TRUE,padj=2)
        }
    } else {
        # step length and turning angle are provided

        if(compact) {
            ################################
            ## Map of all animals' tracks ##
            ################################
            par(mfrow = c(1,1))

            # determine bounds
            ind <- which(data$ID %in% unique(data$ID)[animalsInd])
            xmin <- min(data$x[ind],na.rm=T)
            xmax <- max(data$x[ind],na.rm=T)
            ymin <- min(data$y[ind],na.rm=T)
            ymax <- max(data$y[ind],na.rm=T)
            # make sure that x and y have same scale
            if(xmax-xmin>ymax-ymin) {
                ymid <- (ymax+ymin)/2
                ymax <- ymid+(xmax-xmin)/2
                ymin <- ymid-(xmax-xmin)/2
            } else {
                xmid <- (xmax+xmin)/2
                xmax <- xmid+(ymax-ymin)/2
                xmin <- xmid-(ymax-ymin)/2
            }

            if(length(animalsInd)>7) {
                # to make sure that all colours are distinct (emulate ggplot default palette)
                hues <- seq(15, 375, length = length(animalsInd) + 1)
                colors <- hcl(h = hues, l = 65, c = 100)[1:length(animalsInd)]
            } else # color-blind friendly palette
                colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

            ID <- unique(data$ID)[animalsInd[1]]
            x <- data$x[which(data$ID==ID)]
            y <- data$y[which(data$ID==ID)]
            # plot the first track
            plot(x,y,type="o",pch=20,lwd=1.3,col=colors[1], cex=0.5,
                 xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="x",ylab="y")

            # add other tracks
            for(zoo in animalsInd[-1]) {
                ID <- unique(data$ID)[zoo]
                x <- data$x[which(data$ID==ID)]
                y <- data$y[which(data$ID==ID)]
                points(x,y,type="o",pch=20,lwd=1.3,col=colors[zoo],cex=0.5)
            }
        }

        for(zoo in animalsInd) {
            ID <- unique(data$ID)[zoo]
            x <- data$x[which(data$ID==ID)]
            y <- data$y[which(data$ID==ID)]
            step <- data$step[which(data$ID==ID)]
            angle <- data$angle[which(data$ID==ID)]

            if(!compact) {
                ################################
                ## Map of each animal's track ##
                ################################
                par(mfrow=c(1,1))

                # map of the animal's track
                plot(x,y,type="o",lwd=1.3,xlab="x",ylab="y",pch=20,cex=0.5,asp=1)

                mtext(paste("Animal ID:",ID),side=3,outer=TRUE,padj=2)
            }

            ##################################
            ## Steps and angles time series ##
            ##################################
            par(mfrow=c(2,2))
            plot(step,type="h",xlab="t",ylab="step length",
                 ylim=c(0,max(step,na.rm=T)))
            plot(angle,type="h",xlab="t",ylab="turning angle (radians)",
                 ylim=c(-pi,pi),yaxt="n")
            axis(2, at = c(-pi, -pi/2, 0, pi/2, pi),
                 labels = expression(-pi, -pi/2, 0, pi/2, pi))
            abline(h=c(-pi,0,pi),lty=2)

            #################################
            ## Steps and angles histograms ##
            #################################
            hist(step,xlab="step length",main="", col="grey",border="white",breaks=breaks)

            h <- hist(angle,breaks=breaks,plot=FALSE) # to define the breaks

            hist(angle,xlab="turning angle (radians)",main="", col="grey",border="white",
                 breaks=seq(-pi,pi,length=length(h$breaks)),xaxt="n")
            axis(1, at = c(-pi, -pi/2, 0, pi/2, pi),
                 labels = expression(-pi, -pi/2, 0, pi/2, pi))

            mtext(paste("Animal ID:",ID),side=3,outer=TRUE,padj=2)
        }
    }

    # set graphical parameters back to default
    par(ask=FALSE)
    par(mfrow=c(1,1))
    par(mar=c(5,4,4,2)+0.1)
}
