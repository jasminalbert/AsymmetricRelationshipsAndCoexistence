source("./decomp1_plotting.R")

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
fig3_loc <- paste(fig_loc,"fig3.pdf",sep="")
fig3qij_loc <- paste(fig_loc,"fig3_qij.pdf",sep="")

load(params_loc)

##FIGURE 3##
#plotting decomposed mechanisms contributions against sigma
#for different deltas and mu1-mu2

nd <- length(delta) 
nm <- length(mudif)
ph <- 3 #panel height
pw <- 2.5 #panel widht
ht <- c(1,rep(ph, nm)) #vector of panel heights
wd <- c(rep(pw, nd),2) #vector of panel widths
midx <- (pw*(nd/2))/sum(wd) #middle of figure panels; horizontal
midy <- (ph*(nm/2))/sum(ht) #middle of figure panels; vertical 

#########################################################################################
pdf(fig3_loc)

par(mgp=c(3,0.6,0), mar = c(0.5,1,1,1), oma=c(4,4,2,2), xpd=NA)

layout(matrix(c(2,3,4,5,1,
                6,7,8,9,1,	
                10,11,12,13,1,
                14,15,16,17,1,
                18,19,20,21,1), ncol=5, byrow=TRUE), 
       heights=ht, widths=wd)

plot.new() #1
legend("topright", 
       legend=c(expression(Delta[i]^0), expression(Delta[i]^E),
                expression(Delta[i]^C), expression(Delta[i]^"(E#C)"),
                expression(Delta[i]^"[E||C]"), expression(Delta[i]^"[EC]")
                ,expression(IGR)), 
       col = c("black","black","black","black","blue","red","orange"),
       lty = c(1,2,4,3,1,1,1), bty="n", cex=1.8, inset=c(0,-0.05),
       y.intersp = 1.1, x.intersp = 0.1, seg.len=0.8, lwd=1.5)
#2-5
for (d in 1:nd){
  plot.new()
}

#6-21
res <- vector(mode='list',length=1)
n <- 1
for (m in 1:nm){
  for (d in 1:nd){
    res <- append(res,dePlot1(noise_loc, mudif[m], delta[d], xaxt="n", cex.axis=1.2))
    axis(1, labels=ifelse(n>12, yes=TRUE, no=FALSE), tick=TRUE, cex.axis=1.3)
    mtext(paste0("(", letters[n],")"), side=3, line=-1.7, at=0, cex=1.3, adj=0)
    
    if(n<5){
      mtext(paste(delta[d]), side=3, line=0.2, font=2, cex=1)
    }
    if(n%%4==0){
      mtext(paste(mudif[m]), side=4, line=0.6, font=2, cex=1)
    }
    n <- n+1
  }
}

mtext("contribution to coexistence", side=2, outer=TRUE, line=1.3, font=2, cex=1.3, at=midy-0.005)
mtext(expression(mu[1]-mu[2]), side=4, outer=TRUE, line=-5.5, font=2, cex=2, at=midy)
mtext(expression(delta), side=3, outer=TRUE, line=-2.5, font=2, cex=2, at=midx)
mtext(expression(sigma), outer=TRUE, side=1, line=2, cex=1.5, at=midx)

dev.off()
fig3maxse <- max(unlist(lapply(res, function(X){X$D_se})), na.rm = TRUE)
cat("maximum standard error in figure three is", fig3maxse, "\n(M=", M, ")\n")
fig3maxse_loc <- paste(numeric_results_loc, "/fig3maxse.RDS", sep="")
saveRDS(fig3maxse, file=fig3maxse_loc)
#Sys.time()
#####################################################################################
pdf(fig3qij_loc)

par(mgp=c(3,0.6,0), mar = c(0.5,1,1,1), oma=c(4,4,2,2), xpd=NA)

layout(matrix(c(2,3,4,5,1,
                6,7,8,9,1,	
                10,11,12,13,1,
                14,15,16,17,1,
                18,19,20,21,1), ncol=5, byrow=TRUE), 
       heights=ht, widths=wd)

plot.new() #1
legend("topright", 
       legend=c(expression(Delta[i]^0), expression(Delta[i]^E),
                expression(Delta[i]^C), expression(Delta[i]^"(E#C)"),
                expression(Delta[i]^"[E||C]"), expression(Delta[i]^"[EC]")
                ,expression(IGR)), 
       col = c("black","black","black","black","blue","red","orange"),
       lty = c(1,2,4,3,1,1,1), bty="n", cex=1.8, inset=c(0,-0.05),
       y.intersp = 1.1, x.intersp = 0.1, seg.len=0.8, lwd=1.5)
#2-5
for (d in 1:nd){
  plot.new()
}

#6-21
res <- vector(mode='list',length=1)
n <- 1
for (m in 1:nm){
  for (d in 1:nd){
    res <- append(res,dePlot1(noise_loc, mudif[m], delta[d], qij=TRUE, xaxt="n", cex.axis=1.2))
    axis(1, labels=ifelse(n>12, yes=TRUE, no=FALSE), tick=TRUE, cex.axis=1.3)
    mtext(paste0("(", letters[n],")"), side=3, line=-1.7, at=0, cex=1.3, adj=0)
    
    if(n<5){
      mtext(paste(delta[d]), side=3, line=0.2, font=2, cex=1)
    }
    if(n%%4==0){
      mtext(paste(mudif[m]), side=4, line=0.6, font=2, cex=1)
    }
    n <- n+1
  }
}

mtext("contribution to coexistence", side=2, outer=TRUE, line=1.3, font=2, cex=1.3, at=midy-0.005)
mtext(expression(mu[1]-mu[2]), side=4, outer=TRUE, line=-5.5, font=2, cex=2, at=midy)
mtext(expression(delta), side=3, outer=TRUE, line=-2.5, font=2, cex=2, at=midx)
mtext(expression(sigma), outer=TRUE, side=1, line=2, cex=1.5, at=midx)

dev.off()
fig3qijmaxse <- max(unlist(lapply(res, function(X){X$Dq_se})), na.rm = TRUE)
cat("maximum standard error in figure three (qij) is", fig3qijmaxse, "\n(M=", M, ")\n")
fig3qijmaxse_loc <- paste(numeric_results_loc, "/fig3qijmaxse.RDS", sep="")
saveRDS(fig3qijmaxse, file=fig3qijmaxse_loc)
#Sys.time()