#Fichier de script parmettant de calculer l'enveloppe de la molécule.
#Librairie "alphashape3d" utilisée pour calculer l'enveloppe.
library("alphashape3d")

Rinashape3d <- function(triang, x, alpha, points) {
  #print(triang)
  
	triang = matrix(triang, ncol=9)
	triang[triang[,9] == 3,9] = 2
	colnames(triang) <- c("tr1", "tr2", "tr3", "on.ch", "attached", "rhoT", "muT", "MuT", "fc:2")
	x = matrix(x, ncol=3)
	#print(x)
	#as3d <- ashape3d(x, alpha = alpha)
	#plot(as3d)
	#Sys.sleep(20)
	
	points = matrix(points, ncol=3)
	as3d <- list(triang = triang, x = x, alpha = alpha)
	#print(as3d)
	print(points)
	insh <- inashape3d(as3d, 1, points=points)
	print(insh)
	return (insh) 
}

Rashape3d <- function(data, alpha) {

	data = matrix(data, ncol=3)
	#print(data)
	as3d <- ashape3d(data, alpha = alpha)
	#plot(as3d)
	#Sys.sleep(20)
	
	ma <- list(
		triang=as3d$triang[as3d$triang[, 9] == 2 | as3d$triang[, 9] == 3,],
	 	edge=as3d$edge[as3d$edge[,8] == 2 | as3d$edge[,8] == 3, c("ed1", "ed2")],
	  vertex=as3d$vertex[as3d$vertex[,5] == 2 | as3d$vertex[,5] == 3, c("v1")],
	  x=as3d$x[as3d$vertex[,5] == 2 | as3d$vertex[,5] == 3,],
	  alpha=as3d$alpha
	)
	pr <- list(
		triang=as3d$triang,
	 	edge=as3d$edge,
	  vertex=as3d$vertex,
	  x=as3d$x,
	  alpha=as3d$alpha
	)
	
	ret<- list(
	  ma=ma,
	  pr=pr
	)
	#points <- c(-3.8818, 7, 2.8605, -4, 3.4921, 2)
	#insh = Rinashape3d(as3d$triang, as3d$x, as3d$alpha, points)
	
	#print(ret)
	return (ret)
}
