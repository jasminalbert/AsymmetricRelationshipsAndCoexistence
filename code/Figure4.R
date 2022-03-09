source("./decomp2_plotting.R")

numeric_results_loc <- "../results_numeric"
if(dir.exists(numeric_results_loc)==FALSE){
  dir.create(numeric_results_loc)
}
params_loc <- paste(numeric_results_loc, "/params.RData", sep="")
noise_loc <- paste(numeric_results_loc, "/noise.RData", sep = "")
fig4numres_loc <- paste(numeric_results_loc, "/fig4numres_", sep="")

fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
fig4_loc <- paste(fig_loc,"fig4.pdf",sep="")

load(params_loc)

mudif_4 = seq(range(mudif)[2], range(mudif)[1]-0.5, length.out = 100)

#FIGURE 4# 
#plotting IGR and IGR-ATA to find ATA exclusion and ATA rescue in trait space (plot against mudif)
#for different sigma and delta

nd <- length(delta) 
ns <- length(sigma)
ph <- 3 #panel height
pw <- 2.5 #panel widht
ht <- c(1,rep(ph, ns)) #vector of panel heights
wd <- c(rep(pw, nd),2) #vector of panel widths
midx <- (pw*(nd/2))/sum(wd) #middle of figure panels; horizontal
midy <- (ph*(ns/2))/sum(ht) #middle of figure panels; vertical 

#########################################################################################
pdf(fig4_loc)

par(mgp=c(3,0.5,0), mar = c(0.5,1,1,1), oma=c(4,4,2,2), xpd=NA)

layout(matrix(c(2,3,4,5,1,
                6,7,8,9,1,	
                10,11,12,13,1,
                14,15,16,17,1,
                18,19,20,21,1), ncol=5, byrow=TRUE), 
       heights=ht, widths=wd)

plot.new() #1
legend("topright", 
       legend=c(expression(GWR),expression(GWR-Delta[i]^"[EC]"), expression(Delta[i]^"[EC]")),
       col = c("orange","navy", "red"),
       lty = c(1,1,1), bty="n", cex=1.8, inset=c(-0.51,-0.05),
       y.intersp = 1.1, x.intersp = 0.05, seg.len=0.8, lwd=1.5)
legend("topright", legend=c("ATA \nrescue", "ATA \nexclusion"), density=50, fill=c("darkgoldenrod2","hotpink2"),
       bty="n", border="white", inset=c(-0.33, 0.09), cex=1.6, y.intersp = 1.6,
       x.intersp = 0.3)
#2-5
for (d in 1:nd){
  plot.new()
}

#6-21
res <- vector(mode='list',length=1)
m <- 1
for (s in 1:ns){
  for (d in 1:nd){
    res_loc <- paste(fig4numres_loc,m,".RDS", sep="")
    res <- append(res, list(dePlot2(noise_loc, res_exist=ifelse(file.exists(res_loc), yes=TRUE, no=FALSE), res_loc=res_loc, mudif_4,sigma[s], delta[d], xaxt="n", cex.axis=1.3)))
    axis(1, labels=ifelse(m>12, yes=TRUE, no=FALSE), tick=TRUE, cex.axis=1.2)
    mtext(paste0("(", letters[m],")"), side=3, line=-1.7, at=-4.6, adj=1, cex=1.3)
    
    if(m<5){
      mtext(paste(delta[d]), side=3, line=0.2, font=2, cex=1)
    }
    if(m%%4==0){
      mtext(paste(sigma[s]), side=4, line=0.6, font=2, cex=1.1)
    }
    if(m==12){
      abline(v=-0.8, col="red", lty=2, xpd=T)
    }
    m <- m+1
  }
}

mtext("contribution to coexistence", side=2, outer=TRUE, line=1.3, font=2, cex=1.3, at=midy-0.005)
mtext(expression(sigma), side=4, outer=TRUE, line=-5.5, font=2, cex=2, at=midy)
mtext(expression(delta), side=3, outer=TRUE, line=-2.5, font=2, cex=2, at=midx)
mtext(expression(mu[1]-mu[2]), outer=TRUE, side=1, line=2, cex=1.5, at=midx)

dev.off()

fig4maxse <- max(unlist(lapply(res[-1], function(Y){lapply(Y, function(X){X$D_se})})), na.rm=TRUE)
cat("\nmaximum standard error in figure four is", fig4maxse, "\n(M=", M, ")\n")
fig4maxse_loc <- paste(numeric_results_loc, "/fig4maxse.RDS", sep="")
saveRDS(fig4maxse, file=fig4maxse_loc)
Sys.time()
