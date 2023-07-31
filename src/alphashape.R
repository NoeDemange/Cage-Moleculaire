#Script File for Calculating the Envelope of the Molecule.
#Library "alphashape3d" is used to compute the envelope.
library("alphashape3d", attach.required = FALSE)

#' @title Rashape3d
#' @description Calculates the envelope of a molecule using the "alphashape3d" library.
#'
#' This function calculates the envelope of a 3D molecule using the "alphashape3d" library.
#' The input data should be provided in a matrix format, where each row represents a 3D point.
#' The calculated envelope is returned as a list containing relevant information about the envelope.
#'
#' @param data A matrix representing 3D points of the molecule.
#' @param alpha The alpha value used for the alphashape calculation.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{triang}{A matrix representing the triangles forming the envelope.}
#'   \item{edge}{A matrix representing the edges of the envelope.}
#'   \item{vertex}{A matrix representing the vertices of the envelope.}
#'   \item{x}{A matrix containing the coordinates of the vertices.}
#'   \item{alpha}{The alpha value used in the alphashape calculation.}
#' }
#'
#' @seealso \code{\link{ashape3d}}
#'
#' @importFrom alphashape3d ashape3d
Rashape3d <- function(data, alpha) {

	# Convert input data to a matrix with 3 columns
	data = matrix(data, ncol=3)

	# Calculate the alpha shape using ashape3d function from alphashape3d library
	as3d <- ashape3d(data, alpha = alpha)
	
	# Extract the relevant information from the alpha shape result
	ret <- list(
		triang=as3d$triang[as3d$triang[, 9] == 2 | as3d$triang[, 9] == 3,],
	 	edge=as3d$edge[as3d$edge[,8] == 2 | as3d$edge[,8] == 3, c("ed1", "ed2")],
	  vertex=as3d$vertex[as3d$vertex[,5] == 2 | as3d$vertex[,5] == 3, c("v1")],
	  x=as3d$x[as3d$vertex[,5] == 2 | as3d$vertex[,5] == 3,],
	  alpha=as3d$alpha
	)

	return (ret)
}
