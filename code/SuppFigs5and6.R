#If desired, this script can be used to produce versions of 
#fig 5 and fig 6 where Cyclotella psuedosteligaria 
#is the rare species instead of Fragilaria crotonensis

##libraries used (invoked with ::): none

### sourcing ###
source("./diatom/diatomDecomp_fxns.R")
source("./diatom/fig5_fxns.R")
source("./diatom/fig5def_fxns.R")

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
	} else if (sum(is.na(Effmat))==length(Effmat)){
		pl <- c(NA,NA)
	} else {
		bottom <- cbind(rev(cols),rows[1],1,-1)
		left <- cbind(cols[1], rows,-1,1)[-1,]
		top <- cbind(cols, rows[length(rows)],-1,-1)[-1,]
		right <- cbind(cols[length(cols)], rev(rows),1,-1)[-1,]
	
		border <- rbind(bottom, left, top, right)
		x <- NA; y <- NA; cx <- 1; cy <- 1
		for (p in 1:nrow(border)){
			row <- which(rownames(Effmat)==border[p,2])
			col <- which(colnames(Effmat)==border[p,1])
			if (!is.na(Effmat[row,col])){
				x[cx] <- border[p,1]; cx <- cx+1
				y[cy] <- border[p,2]; cy <- cy+1
			} else if (!is.na(Effmat[row+border[p,3], col+border[p,4]])){
				x[cx] <- cols[col+border[p,4]]; cx <- cx+1
				y[cy] <- rows[row+border[p,3]]; cy <- cy+1 
			}
		}
		pl <- cbind(x,y)
	}
	return(pl)	
}	

### location to save figs ###
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}

### define parameters (invader=2 this time) ###
parms <- c('a'=6, 'P'=60, 'Tbar'=18, 'time'=3000, 'reps'=200, 'invader'=2)


### supplemental version of fig 5 ###

## location of numeric results ## 
numRes_loc <- "../results_numeric/"
dat_loc <- paste0(numRes_loc, "fig5dat2/")

#if data is missing, make them
if (dir.exists(dat_loc)==FALSE){
  dir.create(dat_loc)
  
  #define sets of variables
  a <- seq(1,6,length.out=100)
  Tbar <- seq(16,18,length.out = 100)
  P <- seq(51,199.5,length.out=100)
  
  #makes results as folder of .RDS files
  dat5(dat_loc, a_vec=a, P_vec=P, T_vec=Tbar, parms=parms)
}
## name of figure file ##
fig5_2 <- paste(fig_loc,"fig5_2.pdf",sep="")

## makes figure ##
fig5(fig5_2, dat_loc, invader=2) 

### get and save standard error ###
fig5_2se <- getSE(dat_loc)
fig5_2maxse <- max(fig5_2se)
cat("\nmaximum standard error in figure five (SI) is", fig5_2maxse)
fig5_2maxse_loc <- paste0(numeric_results_loc, "fig5_2maxse.RDS")
saveRDS(fig5_2maxse, file=fig5_2maxse_loc)


### supplemental version of fig 6 ###

## location of numeric results ## 
dat_loc <- paste0(numRes_loc, "fig6dat2/")
dat_loc <- paste0(numRes_loc, "out91923/")

#if data is missing, make them
if (dir.exists(dat_loc)==FALSE){
  dir.create(dat_loc)
  
  #define sets of variables
  a <- seq(3.5,6,length.out=100)
  P <- seq(51,199.5,length.out=100)
  Tbar <- seq(16,19,length.out=100)
  
   #makes results as folder of .RDS files
  dat6(dat_loc, a_vec=a, P_vec=P, T_vec=Tbar, parms)
}

#setting up figure

xat <- list(c(16, 16.5, 17, 17.5, 18), seq(4,6,len=5),c(16, 16.5, 17, 17.5, 18))
xlab <- list(c(16, NA, NA, NA, 18), c(50, 67.5, NA, 102.5, 120),c(50, 67.5, NA, 102.5, 120))
yat <- list(c(4,5,6), seq(50,120,len=3),seq(50, 120, len=5))
ylab <- list(c(4,NA,6), c(4,NA,6),T=c(16, NA, NA, NA, 18))

labs <- c("amplitude (a)", "Period (P)", expression(paste("mean temperature (",theta[0],")" ) ))
ncol <- 100
col1 <- hcl.colors(ncol,"YlOrRd", rev=T, alpha=.8)
## name of figure file ##
fig6_2 <- paste0(fig_loc,"fig6_2.pdf")


files <- list.files(dat_loc)
mats <- list()
rmats <- list()
for (f in seq_along(files)){
	dat <- readRDS(paste0(dat_loc,files[f]))$dat
	GWRmat <- MAT(dat, dat$GWR)
	ATAmat <- MAT(dat, dat$ATA)
	mat <- ATAmat/abs(GWRmat)
	rmats[[f]] <- mat
	mats[[f]] <- list(ATA=ATAmat, GWR=GWRmat)
}
ranges <- log(rbind(sapply(rmats, min), rep(10,3)))
diffs <- apply(ranges, 2, diff)
standiff <- round(diffs/max(diffs),2)
stanRan <- rbind(ranges[1,]/min(ranges[1,]), ranges[2,]/max(ranges[2,]))
stanRan <- round(stanRan,2)
col1 <- col1
col2 <- col1[(ncol*(1-standiff[2])):ncol]
col3 <- col1[(ncol*(1-standiff[3])):ncol]
cols <- list(col1,col2,col3)
mats[[1]] <- lapply(mats[[1]],t)
pdf(fig6_2, width=5)
par(mfrow=c(3,1), mar=c(2,1,1,0), oma=c(.5,3,0,3),mgp=c(2,.5,0))
for (m in 1:length(mats)){
	mat <- mats[[m]]
	map1(t(mat$ATA), t(mat$GWR), col=col1)#cols[[m]]) 
	labs <- names(dimnames(mat[[1]]))
	graphics::title(ylab=labs[2], cex.lab=1.3, line=1.1, xpd=NA)
	graphics::title(xlab=labs[1], cex.lab=1.3, line=1.1, xpd=NA)
	graphics::axis(1, lwd.ticks=2); 
	graphics::axis(2, lwd.ticks=2,mgp=c(2,.5,0))
	#graphics::mtext(paste0("(", letters[5],")"), 3, -1.5, adj=0.985, cex=1.1)
}
dev.off()

## makes figure ##
#fig6(fig6_2, dat_loc, invader=2) 

### get and save standard error ###
fig6_2se <- getSE(dat_loc)
fig6_2maxse <- max(fig6_2se)
cat("\nmaximum standard error in figure six (SI) is", fig6_2maxse)
fig6_2maxse_loc <- paste0(numeric_results_loc, "fig6_2maxse.RDS")
saveRDS(fig6_2maxse, file=fig6_2maxse_loc)








