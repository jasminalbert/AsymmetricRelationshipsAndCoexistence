cmMap <- function(dat){
	dat$map[dat$map>1] <- rep(1, length(dat$map[dat$map>1]))
	dat$map[dat$map<(-1)] <- rep(-1, length(dat$map[dat$map<(-1)]))
	
	a <- unique(dat$a); P <- unique(dat$P); Tbar <- unique(dat$Tbar)
	dims <- list(a=round(a,3), P=round(P,3), Tbar=round(Tbar,3))
	
	len <- unlist(lapply(dims, length))
	vars <- dims[len!=1]
	
	map <- matrix(dat$map, nrow=length(vars[[1]]), ncol=length(vars[[2]]), dimnames=vars, byrow=T)
	
	return(map)
}

cmContour <- function(map, ncolor=51){
	x <- dimnames(map)[2]
	y <- dimnames(map)[1]
	main <- expression(paste(Delta[i]^"[EC]","/IGR"))
	cm <- cm.colors(ncolor)
	
	image2D(z=t(map), y=as.numeric(y[[1]]), x=as.numeric(x[[1]]), contour=FALSE, col=cm, xlab=names(x), ylab=names(y), main=main)
	
	contour(z=t(map), y=as.numeric(y[[1]]), x=as.numeric(x[[1]]),add=TRUE, col='grey50')
}

#test
aP4 <- readRDS("ampP4.RDS")
map4 <- cmMap(aP4)
cmContour(map4)