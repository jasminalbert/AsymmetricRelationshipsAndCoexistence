#figure 5 panels e-f plotting

##libraries used (invoked with ::): plot3D, graphics


#par(mar=c(2,2.2,2,2.2), mgp=c(2,.8,0))
#par(mar=c(3,4,3,2))

#makes matrix from col of data.frame and param values as ij
MAT <- function(dat, vals){
	a <- unique(dat$a); 
	P <- unique(dat$P); 
	Tbar <- unique(dat$Tbar)
	dims <- list(a=round(a,2), P=round(P,2), 	Tbar=round(Tbar,2))
	len <- unlist(lapply(dims, length))
	vars <- dims[len!=1]

	mat <- matrix(vals, nrow=length(vars[[1]]), ncol=length(vars[[2]]), dimnames=vars)
	return(mat)
}

#GWRmat <- MAT(aTb, aTb$GWR)
#ATAmat <- MAT(aTb, aTb$ATA)

#ck	plot colkey?
map <- function(mat, add=F, col=hcl.colors(2, "Light Grays"), ck=list(plot=FALSE), ct=F){
	z <- t(mat); 
	x <- dimnames(mat)[2];y <- dimnames(mat)[1]
	plot3D::image2D(z=z, y=as.numeric(y[[1]]), x=as.numeric(x[[1]]), contour=F, xlab='', ylab='', col=col, add=add,NAcol=rgb(0,0,0,0), colkey=ck, xaxt='n', yaxt='n')
	if (ct==T){
		nlvl <- 8
		zlim <- range(z, finite=T)
		lvls <- pretty(zlim, nlvl)
		plot3D::contour2D(z=z, y=as.numeric(y[[1]]), x=as.numeric(x[[1]]), add=T, col="black", drawlabels=F, nlevels=nlvl, colkey=list(plot=F), lwd=par("lwd")+0.2)
		plot3D::contour2D(z=z, y=as.numeric(y[[1]]), x=as.numeric(x[[1]]), add=T, col=col[round(seq(1,length(col), len=length(lvls)))], drawlabels=F, nlevels=nlvl,colkey=F)#list(cex.axis=0.8, mgp=c(3,.7,0))
		plot3D::contour2D(z=z, y=as.numeric(y[[1]]), x=as.numeric(x[[1]]), add=T, col="black", labels=round(exp(lvls),5), nlevels=nlvl, colkey=list(plot=F), lty=0)
	} 
}; 
#coex	T:draw coexistence, F:draw ATA rescue
draw <- function(ATAmat, GWRmat, coex=T){
	Effmat <- GWRmat
	
	if (coex==T){
		Effmat[Effmat<0] <- NA
	} else {
		Effmat[!(ATAmat>GWRmat & GWRmat>0) ] <- NA}
	
	rows <- as.numeric(dimnames(Effmat)[[1]]);
	cols <- as.numeric(dimnames(Effmat)[[2]])
	if ( !sum(is.na(Effmat)) ){
		verts <- expand.grid(range(cols),range(rows))
		verts[3:4,] <- verts[4:3,]
		pl <- as.matrix(verts)
	} else {
		xu <- NA; xl <- NA; xr <- NA; xlft <- NA 
		yu <- NA; yl <- NA; yr <- NA; ylft <- NA
		cu <- 1; cl <- 1; cr <- 1; clft <- 1
		nr <- nrow(Effmat); nc <- ncol(Effmat)
		for (j in nc:1){
			for(i in nr:1){
				if(!is.na(Effmat[i,j])){
					if(j==1){
						ylft[clft] <- i; 
						xlft[clft] <- j; clft<-clft+1
					}
					if(j==nc){
						yr[cr] <- i; 
						xr[cr] <- j; cr<-cr+1
					}
					if(i==nr){
						yu[cu] <- i; 
						xu[cu] <- j; cu<-cu+1
					} else if(is.na(Effmat[i+1,j])){
						yu[cu] <- i; 
						xu[cu] <- j; cu<-cu+1
					} else if(i==1){
						yu[cu] <- i; 
						xu[cu] <- j; cu<-cu+1
					} else if(is.na(Effmat[i-1,j])){
						yl[cl] <- i; 
						xl[cl] <- j; cl<-cl+1
					}
				}
			}
		}
	#pl <- matrix(c(cols[xu], rev(cols[xlft]), rev(cols[xl]), cols[xr], rows[yu], rev(rows[ylft]), rev(rows[yl]), rows[yr]), ncol=2)
		pl <- matrix(c(cols[xu], NA, rev(cols[xl]), cols[xr], rows[yu], NA, rev(rows[yl]), rows[yr]), ncol=2)
		pl <- pl[!is.na(pl)[,1],]
	}
	return(pl)
}
#maps ATA effect / abs GWR
map1 <- function(ATAmat, GWRmat, col, cob=T){
	mat <- ATAmat/abs(GWRmat)	
	mat1 <- mat; mat1[mat1<0] <- NA
	mat1[mat1>10] <- 10
	lmat1 <- log(mat1)
	r1 <- range(mat1, na.rm=T)
	
	#mat2 <- mat; mat2[mat2>0] <- NA
	#absmat2 <- abs(mat2)
	#labsmat2 <- log(absmat2)
	
	coex <- draw(ATAmat, GWRmat)
	ATA <- draw(ATAmat, GWRmat, coex=F)
	
	#col <- hcl.colors(n,"Geyser")
	#col1 <- col[(n/2 + 1):n]; col2 <- col[(n/2):round(n/4)]
	
	map(log(mat1), col=col, ct=T, ck=F)
	if(cob==T){
		graphics::polygon(coex, border=rgb(0,0,1,0.6), lwd=2)
		graphics::polygon(ATA, density=15, border=NA, col=rgb(0,0,0,0.6))
	}
	#if(sum(is.na(mat2))<length(mat2)){
	#	map(labsmat2, T, col=col2, ct=T, ck=F)
	#	r2 <- range(absmat2, na.rm=T)
	#}
	return(list(col1=col1, col2=rev(col2), lim1=r1))
}; 





#mat1 <- ATAmat/abs(GWRmat); range(mat1)
#col1 <- hcl.colors(100,"Geyser")
#map(log(mat1), col=col1[51:100])
#dim(mat1)

#checks
#map1(ATAmat[10:30,10:30], GWRmat[10:30,10:30], 100)
#map1(ATAmat[31:100,1:70], GWRmat[31:100,1:70], 100)
#map1(ATAmat[1:50,1:50], GWRmat[1:50,1:50], 100)


#map(GWRmat, ck=NULL)
#map(ATAmat)
#map(ATAmat/abs(GWRmat))
#range(ATAmat/abs(GWRmat))

#hist(mat)

#map(ATAmat, col=col,F)
#map(abs(GWRmat), col=col,F)
#map(ATAmat/abs(GWRmat), col=col,F)
#range(ATAmat/abs(GWRmat))
#mat22 <- mat2[!is.na(mat2)]
#mat11 <- mat1[!is.na(mat1)]
#range(log(-mat22))
#range(log(mat11))











