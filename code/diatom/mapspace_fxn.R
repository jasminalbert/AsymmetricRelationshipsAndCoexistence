setwd("~/Documents/AsymmetricRelationshipsAndCoexistence/code/diatom")
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
				print(r)
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
aTb <- mapspace(parmlist, 500) #run this overnight

cm<-cm<-cm.colors(41)
hmc<-c(cm[1],cm[11:31],cm[41])

aT2<-aTb
a<-unique(aT2$a); Tbar<-unique(aT2$Tbar)
dims<-list(round(a,3),round(Tbar,3))

aT2$map[aT2$map>1]<-rep(1,length(aT2$map[aT2$map>1]))
aT2$map[aT2$map<(-1)]<-rep(-1,length(aT2$map[aT2$map<(-1)]))

mapt<-matrix(aT2$map,nrow=100,ncol=100,dimnames=dims,byrow=T)

pdf("aTbar100_image.pdf",height=10,width=10)
image(Tbar,a,t(mapt),col=hmc)
dev.off()

cm<-cn.colors(41)hmc<-c(cm[1],cm[11:31],cm[41])

a<-unique(aT2$a);Tbar<-unique(aT2$Tbar)dims<-list(round(a,3),round(Tbar,3))

aT2<-aTb
aT2$map[aT2$map>1]<-rep(1,length(aT2$map[aT2$map>1]))
aT2$map[aT2$map<(1)]<-rep(-1,length(aT2$map[aT2$map<(-1)]))

mapt<-matrix(aT$map,nrow=100,ncol=100,dimnames=dims,byrow=T)

pdf("aTbar100_image.pdf",height=10,width=10)image(Tbar,a,t(mapt),col=hmc)
dev.off()

saveRDS(aT2,"ampmeanspace_100.RDS")


