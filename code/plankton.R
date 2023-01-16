#looking at Ceratium species time series
dat_loc <- "../plankton/planktontimeseries"

#extract C. fusus data and temp from region 25 of link data*
#plot against negative temp
#temp =E. C. fusus= competition
#n=12 in Gosh2020 is area 25 in Sheppard2019

#C. fusus data (sp1)
Cerfus_loc <- paste0(dat_loc, "/Ceratium_fusus.csv")
cerfus <- read.csv(Cerfus_loc)
head(cerfus)
str(cerfus)
nrow(cerfus)

#C. furca (sp2)
Cerfur_loc <- paste0(dat_loc, "/Ceratium_furca.csv")
cerfur <- read.csv(Cerfur_loc)

#C. tripos (sp3)
Certri_loc <- paste0(dat_loc, "/Ceratium_tripos.csv")
certri <- read.csv(Certri_loc)

#C. macroceros (sp4)
Cermac_loc <- paste0(dat_loc, "/Ceratium_macroceros.csv")
cermac <- read.csv(Cermac_loc)

#predator Calanus finmarchicus
Calfin_loc <- paste0(dat_loc,"/Calanus_finmarchicus.csv")
calfin <- read.csv(Calfin_loc)

#temp data
#different files for different seasonal averages
temp_loc <- paste0(dat_loc, "/temp1to12.csv")
tempSp_loc <- paste0(dat_loc, "/temp3to5.csv")
tempSpSu_loc <- paste0(dat_loc, "/temp3to9.csv")
tempSu_loc <- paste0(dat_loc, "/temp6to8.csv")
tempF_loc <- paste0(dat_loc, "/temp9to11.csv")

temp <- read.csv(temp_loc)
tempSp <- read.csv(tempSp_loc)
tempSpSu <- read.csv(tempSpSu_loc)
tempSu <- read.csv(tempSu_loc)
tempF <- read.csv(tempF_loc)
head(temp)
dim(temp)==dim(cerfus)
dim_t <- dim(temp)
dimnames <- list(n=1:26, year=1:57, season=c("all","Sp","SpSu","Su","F"))

tempdat <- unlist(c(temp, tempSp, tempSpSu, tempSu, tempF))
tmps <- array(tempdat, dim = c(dim_t[1], dim_t[2],5), dimnames=dimnames)

n <- 25

pdf("cerfus_temp.pdf")
par(mar=(rep(0.8,4)),mfrow=c(5,2),mgp=c(1.5,0.2,0))
for (S in 1:dim(tmps)[3]){
	x <- unlist(cerfus[,-1][cerfus[,1]==n,])
	y <- -tmps[,-1,S][tmps[,1,]==n]
	T <- length(x)
	plot(x,y)
	#plot(rank(x), rank(y))
	plot(rank(x)/T, rank(y)/T)
	lines(c(min(x), max(x)), c(max(y), min(y)))
}
dev.off()


pdf("cerfus_temp1.pdf")
par(mar=c(2.5,2.5,1,1),mfrow=c(1,1),mgp=c(1.5,0.2,0))
x <- unlist(cerfus[,-1][cerfus[,1]==n,])
y <- -tmps[,-1,S][tmps[,1,]==n]
T <- length(x)
plot(rank(x)/T, rank(y)/T, pch=20, col="darkgrey", cex=2, ylim=c(-.15,1.15),xlim=c(-.15,1.15))
dev.off()

#49 = 26
pdf("cerfus_cerfur_49.pdf")
n <- 49
par(mar=c(2.5,2.5,1,1),mfrow=c(1,1),mgp=c(1.5,0.2,0))
x <- unlist(cerfus[,-1][cerfus[,1]==n,])
y <- unlist(cerfur[,-1][cerfur[,1]==n,])
T <- length(x)
plot(rank(x)/T, rank(y)/T, pch=20, col="darkgrey", cex=2, ylim=c(-.15,1.15),xlim=c(-.15,1.15))
dev.off()

#25=12
pdf("cerfus_temp_25.pdf")
par(mar=c(2.5,2.5,1,1),mfrow=c(1,1),mgp=c(1.5,0.2,0))
n <- 25
x <- unlist(cerfus[,-1][cerfus[,1]==n,])
y <- -tmps[,-1,S][tmps[,1,]==n]
T <- length(x)
plot(rank(x)/T, rank(y)/T, pch=20, col="darkgrey", cex=2, ylim=c(-.15,1.15),xlim=c(-.15,1.15))
dev.off()
#good for use; LT

#49 = 26
pdf("certri_cermac_49.pdf")
par(mar=c(2.5,2.5,1,1),mfrow=c(1,1),mgp=c(1.5,0.2,0))
n <- 49
x <- unlist(certri[,-1][certri[,1]==n,])
y <- unlist(cermac[,-1][cermac[,1]==n,])
T <- length(x)
plot(rank(x)/T, rank(y)/T, pch=20, col="darkgrey", cex=2, ylim=c(-.15,1.15),xlim=c(-.15,1.15))
dev.off()
#good for use; RT

#23=10
pdf("cerfus_cerfur_23.pdf")
par(mar=c(2.5,2.5,1,1),mfrow=c(1,1),mgp=c(1.5,0.2,0))
n <- 23
x <- unlist(cerfus[,-1][cerfus[,1]==n,])
y <- unlist(cerfur[,-1][cerfur[,1]==n,])
T <- length(x)
plot(rank(x)/T, rank(y)/T, pch=20, col="darkgrey", cex=2, ylim=c(-.15,1.15),xlim=c(-.15,1.15))
dev.off()
#good LT

#23=10
pdf("cerfur_temp_23.pdf")
par(mar=c(2.5,2.5,1,1),mfrow=c(1,1),mgp=c(1.5,0.2,0))
n <- 23
x <- unlist(cerfur[,-1][cerfur[,1]==n,])
y <- -tmps[,-1,S][tmps[,1,]==n]
T <- length(x)
plot(rank(x)/T, rank(y)/T, pch=20, col="darkgrey", cex=2, ylim=c(-.15,1.15),xlim=c(-.15,1.15))
dev.off()
#good LT

#40=23
pdf("cerfus_cerfur_40.pdf")
par(mar=c(2.5,2.5,1,1),mfrow=c(1,1),mgp=c(1.5,0.2,0))
n <- 40
x <- unlist(cerfus[,-1][cerfus[,1]==n,])
y <- unlist(cerfur[,-1][cerfur[,1]==n,])
T <- length(x)
plot(rank(x)/T, rank(y)/T, pch=20, col="darkgrey", cex=2, ylim=c(-.15,1.15),xlim=c(-.15,1.15))
dev.off()
#good RT

Rxvec <- list(certri49=unlist(certri[,-1][certri[,1]==49,]),
			cerfus40=unlist(cerfus[,-1][cerfus[,1]==40,]),
			certri40=unlist(certri[,-1][certri[,1]==40,]))
Ryvec <- list(cermac49=unlist(cermac[,-1][cermac[,1]==49,]),
			cerfur40=unlist(cerfur[,-1][cerfur[,1]==40,]),
			calfin40=-unlist(calfin[,-1][calfin[,1]==40,]))
Lxvec <- list(cerfus25=unlist(cerfus[,-1][cerfus[,1]==25,]),
			cerfus23=unlist(cerfus[,-1][cerfus[,1]==23,]),
			cerfur23=unlist(cerfur[,-1][cerfur[,1]==23,]))
Lyvec <- list(temp25=-tmps[,-1,S][tmps[,1,]==25],
			cerfur23=unlist(cerfur[,-1][cerfur[,1]==23,]),
			temp23=-tmps[,-1,S][tmps[,1,]==23])

xlist <- append(Rxvec, Lxvec); ylist <- append(Ryvec, Lyvec)
par(mar=c(2.5,2.5,1,1),mfrow=c(2,3),mgp=c(1.5,0.2,0))
for(p in 1:length(xlist)){
	x <- xlist[[p]]; y <- ylist[[p]]
	T <- length(x)
	plot(rank(x)/T, rank(y)/T, pch=20, col="darkgrey", cex=2, ylim=c(-.15,1.15),xlim=c(-.15,1.15), xlab=names(xlist)[p], ylab=names(ylist)[p])
}

#best RT: cerfus x cerfur 40 (23 in gosh)
loc40 <- data.frame(cerfus=unlist(cerfus[,-1][cerfus[,1]==40,]),cerfur=unlist(cerfur[,-1][cerfur[,1]==40,]))
#best LT: cerfus x cerfur 23 (10 in gosh)
loc23 <- data.frame(cerfus=unlist(cerfus[,-1][cerfus[,1]==23,]),cerfur=unlist(cerfur[,-1][cerfur[,1]==23,]))
ceratium1x2<-list(loc23=loc23,loc40=loc40)
save(ceratium1x2, file="../results_numeric/ceratium1x2.RData")



T <- length(x)
x <- unlist(cerfus[,-1][cerfus[,1]==n,])
y <- -tmps[,-1,S][tmps[,1,]==n]
par(mfrow=c(1,1))
plot(rank(x)/T, rank(y)/T)
abline(1,-1)
xnr <- rank(x)/T; ynr <- rank(y)/T
xl <- xnr[(xnr+ynr)<1]
yl <- ynr[(xnr+ynr)<1]
xu <- xnr[(xnr+ynr)>1]
yu <- ynr[(xnr+ynr)>1]
points(xl,yl, pch=20, col="red", cex=0.5)
points(xu,yu, pch=20, col="blue", cex=0.5)

topl <- sum((xl-mean(xnr))*(yl-mean(ynr)))
btml <- (T-1)*sqrt(var(xnr)*var(ynr))
corl <- topl/btml

topu <- sum((xu-mean(xnr))*(yu-mean(ynr)))
btmu <- (T-1)*sqrt(var(xnr)*var(ynr))
coru <- topu/btmu
corl-coru
cor(xnr,ynr, method="spearman")

hist(ynr); hist(xnr)
hist(qnorm(ynr)); hist(qnorm(xnr))
plot(qnorm(xnr), qnorm(ynr))
xnr <- qnorm(xnr); ynr<-qnorm(ynr)
xnr <- xnr[c(-22,-30)]; ynr <- ynr[c(-22,-30)]
xl <- xnr[(xnr+ynr)<0]
yl <- ynr[(xnr+ynr)<0]
xu <- xnr[(xnr+ynr)>0]
yu <- ynr[(xnr+ynr)>0]
points(xl,yl, pch=20, col="red", cex=0.5)
points(xu,yu, pch=20, col="blue", cex=0.5)


topl <- sum((xl-mean(xnr))*(yl-mean(ynr)))
btml <- (T-1)*sqrt(var(xnr)*var(ynr))
corl <- topl/btml


topu <- sum((xu-mean(xnr, na.rm=T))*(yu-mean(ynr, na.rm=T)))
btmu <- (T-1)*sqrt(var(xnr)*var(ynr))
coru <- topu/btmu
corl-coru;corl+coru
cor(xnr,ynr, method="spearman")













