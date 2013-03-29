####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function to set range of colors for heatplot:


##Input:
### upper.lim, lower.lim: limits for aberration calling
### op: list with other plot parameters

##Output:
### colors: a vector giving a range of color nuances
### intervals: vector giving intervals between lower.lim and upper.lim corresponding to different color nuances

##Required by:
### plotHeatmap 

##Requires:
### none

colorSetup <- function(upper.lim,lower.lim,op){
	#Colors:
	if(is.na(op$n.col)){
    op$n.col <- 50
	}
	#Range of colours in shades:
	colors <- colorRampPalette(op$colors, space="rgb")(op$n.col)

	#Divide into sequence according to values of limits:
	lim.seq <- seq(lower.lim,upper.lim,length.out=op$n.col-2)
	lim.seq <- c(-Inf,lim.seq,Inf)
	n.seq <- length(lim.seq)

	#Get matrix with intervals:
	low <- lim.seq[1:(n.seq-1)]
	high <- lim.seq[2:n.seq]
	intervals <- matrix(c(low,high),ncol=2,nrow=n.seq-1,byrow=FALSE)

	return(list(colors=colors,intervals=intervals))
}

