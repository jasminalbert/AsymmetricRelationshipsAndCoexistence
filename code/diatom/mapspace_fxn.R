source("getEpsilons_fxn.R")
#parmlist	list length 3 of vectors a, P, and Tbar
#			one vector should be constant length 1
mapspace <- function(parmlist, sims){
	
	len <- unlist(lapply(parmlist, length))
	varlen <- len[len!=1]
	vars <- parmlist[len!=1]
	rows <- varlen[[1]]*varlen[[2]]
	res <- matrix(NA, ncol=5, nrow=rows)
	r <- 1
	
	a <- parmlist[[1]]
	P <- parmlist[[2]]
	Tb <- parmlist[[3]]
	
	for (i in 1:length(a)){
		for (j in 1:length(P)){
			for (k in 1:length(Tb)){
				ep1 <- getep2(a[i], P[j], Tb[k], sims, 1)
				ep2 <- getep2(a[i], P[j], Tb[k], sims, 2)
				Delta1 <- ep1-ep2
				
				res[r,] <- c(a[i], P[j], Tb[k], sum(Delta1), Delta1$epsECbrk)
				r <- r+1
			}	
		}
	}
	res <- data.frame(res)
	colnames(res) <- c("a", "P", "Tbar", "IGR", "ATA")
	res$eff <- res$IGR - res$ATA
	res$map <- res$ATA/res$IGR
	
	map <- matrix(res$map, nrow=varlen[[1]], ncol=varlen[[2]], dimnames=lapply(vars, round, 3), byrow=TRUE)
	heatmap(map, Colv=NA, Rowv=NA)
	return(res)
}

parmlist <- list(seq(4,4.5,0.05), seq(170,190,5), 18)
mapspace(parmlist, 200)

a <- seq(3.5, 6, length.out=100)
P <- 60
Tb <- seq(16.5, 19, length.out=100)
parmlist <- list(a,P,Tb)
mapspace(parmlist, 500) #run this overnight





