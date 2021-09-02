source("./decomp2_plotting.R")

numeric_results_loc <- "../results_numeric"
if(dir.exists(numeric_results_loc)==FALSE){
  dir.create(numeric_results_loc)
}
params_loc <- paste(numeric_results_loc, "/params.RData", sep="")
noise_loc <- paste(numeric_results_loc, "/noise.RData", sep = "")

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

par(mgp=c(3,0.5,0), mar = c(1,1,1,1), oma=c(3,5,3,2), xpd=TRUE)

layout(matrix(c(2,3,4,5,1,
                6,7,8,9,1,	
                10,11,12,13,1,
                14,15,16,17,1,
                18,19,20,21,1), ncol=5, byrow=TRUE), 
       heights=ht, widths=wd)

plot.new() #1
legend("topright", 
       legend=c(expression(IGR),expression(IGR-Delta["[E||C]"]), expression(Delta["[E||C]"])),
       col = c("orange","navy", "red"),
       lty = c(1,1,1), bty="n", cex=1.5, inset=c(-0.25,-0.042),
       y.intersp = 1.25, x.intersp = 0.1, seg.len=1)
#2-5
for (d in 1:nd){
  plot.new()
}

#6-21
res <- vector(mode='list',length=1)
m <- 1
for (s in 1:ns){
  for (d in 1:nd){
    res <- append(res,dePlot2(noise_loc, mudif_4,sigma[s], delta[d], xaxt="n"))
    axis(1, labels=ifelse(m>12, yes=TRUE, no=FALSE), tick=TRUE)
    mtext(LETTERS[m], side=3, line=-1.45, at=-4.5, adj=1)
    
    if(m<5){
      mtext(paste(delta[d]), side=3, line=0.75, font=2, cex=0.8)
    }
    if(m%%4==0){
      mtext(paste(sigma[s]), side=4, line=0.75, font=2, cex=0.8)
    }
    if(m==12){
      abline(v=-0.8, col="red", lty=2)
    }
    m <- m+1
  }
}

mtext("contribution to coexistence", side=2, outer=TRUE, line=1.5, font=2, cex=1, at=midy)
mtext(expression(sigma), side=4, outer=TRUE, line=-5, font=2, cex=1.5, at=midy)
mtext(expression(delta), side=3, outer=TRUE, line=-2, font=2, cex=1.5, at=midx)
mtext(expression(mu[1]-mu[2]), outer=TRUE, side=1, line=1, cex.lab=1.3, at=midx)

dev.off()
fig4maxse <- max(unlist(lapply(res, function(X){X$D_se})), na.rm = TRUE)
cat("maximum standard error in figure four is", fig4maxse, "\n(M=", M, ")\n")
fig4maxse_loc <- paste(numeric_results_loc, "/fig4maxse.RDS", sep="")
saveRDS(fig4maxse, file=fig4maxse_loc)
Sys.time()
