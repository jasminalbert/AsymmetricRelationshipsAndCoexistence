load("results_numeric/params.RData")
source("./decomp1_plotting.R")
mudif <- c(0, -0.5, -2, -4)
#FIGURE 3#
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
pdf("results_figs/fig3.2.pdf")

par(mgp=c(3,0.5,0), mar = c(1,1,1,1), oma=c(3,5,3,2), xpd=TRUE)

layout(matrix(c(2,3,4,5,1,
                6,7,8,9,1,	
                10,11,12,13,1,
                14,15,16,17,1,
                18,19,20,21,1), ncol=5, byrow=TRUE), 
       heights=ht, widths=wd)
#layout.show(n=25)

#plot(0, bty="n")
plot.new() #1
legend("topright", 
       legend=c(expression(Delta[i]^0), expression(Delta[i]^E),
                expression(Delta[i]^C), expression(Delta[i]^"(E#C)"),
                expression(Delta[i]^"[EC]"), expression(Delta[i]^"[E||C]")
                ,expression(IGR)), 
       col = c("black","black","black","black","blue","red","orange"),
       lty = c(1,2,4,3,1,1,1), bty="n", cex=1.5, inset=c(0.1,-0.042),
       y.intersp = 1.25, x.intersp = 0.1, seg.len=1)
#2-5
for (d in 1:nd){
  plot.new()
}

#6-21
res <- vector(mode='list',length=1)
n <- 1
for (m in 1:nm){
  for (d in 1:nd){
    res <- append(res,dePlot1(mudif[m], delta[d], xaxt="n"))
    axis(1, labels=ifelse(n>12, yes=TRUE, no=FALSE), tick=TRUE)
    mtext(LETTERS[n], side=3, line=-1.45, at=0.5)
    
    if(n<5){
      mtext(paste(delta[d]), side=3, line=0.75, font=2, cex=0.8)
    }
    if(n%%4==0){
      mtext(paste(mudif[m]), side=4, line=0.75, font=2, cex=0.8)
    }
    n <- n+1
  }
}

mtext("contribution to coexistence", side=2, outer=TRUE, line=1.5, font=2, cex=1, at=midy)
mtext(expression(mu[1]-mu[2]), side=4, outer=TRUE, line=-5, font=2, cex=1.5, at=midy)
mtext(expression(delta), side=3, outer=TRUE, line=-2, font=2, cex=1.5, at=midx)
mtext(expression(sigma), outer=TRUE, side=1, line=1, cex.lab=1.3, at=midx)

dev.off()
fig3.2maxse <- max(unlist(lapply(res, function(X){X$D_se})), na.rm = TRUE)
cat("maximum standard error in figure three is", fig3.2maxse, "\n(M=", M, ")\n")
saveRDS(fig3.2maxse, file="results_numeric/fig3.2maxse.RDS")
#Sys.time()
#####################################################################################
pdf("results_figs/fig3.2_qij.pdf")

par(mgp=c(3,0.5,0), mar = c(1,1,1,1), oma=c(3,5,3,2))

layout(matrix(c(2,3,4,5,1,
                6,7,8,9,1,	
                10,11,12,13,1,
                14,15,16,17,1,
                18,19,20,21,1), ncol=5, byrow=TRUE), 
       heights=ht, widths=wd)
#layout.show(n=25)

#plot(0, bty="n")
plot.new() #1
legend("topright", 
       legend=c(expression(Delta[i]^0), expression(Delta[i]^E),
                expression(Delta[i]^C), expression(Delta[i]^"(E#C)"),
                expression(Delta[i]^"[EC]"), expression(Delta[i]^"[E||C]")
                ,expression(IGR)), 
       col = c("black","black","black","black","blue","red","orange"),
       lty = c(1,2,4,3,1,1,1), bty="n", cex=1.5, inset=c(0.1,-0.025),
       y.intersp = 1.25, x.intersp = 0.1, seg.len=1)
#2-5
for (d in 1:nd){
  plot.new()
}

#6-21
res <- vector(mode='list',length=1)
n <- 1
for (m in 1:nm){
  for (d in 1:nd){
    res <- append(res,dePlot1(mudif[m], delta[d], qij=TRUE, xaxt="n"))
    axis(1, labels=ifelse(n>12, yes=TRUE, no=FALSE), tick=TRUE)
    mtext(LETTERS[n], side=3, line=-1.45, at=0.5)
    
    if(n<5){
      mtext(paste(delta[d]), side=3, line=0.75, font=2, cex=0.8)
    }
    if(n%%4==0){
      mtext(paste(mudif[m]), side=4, line=0.75, font=2, cex=0.8)
    }
    n <- n+1
  }
}

mtext("contribution to coexistence", side=2, outer=TRUE, line=1.5, font=2, cex=1, at=midy)
mtext(expression(mu[1]-mu[2]), side=4, outer=TRUE, line=-5, font=2, cex=1.5, at=midy)
mtext(expression(delta), side=3, outer=TRUE, line=-2, font=2, cex=1.5, at=midx)
mtext(expression(sigma), outer=TRUE, side=1, line=1, cex.lab=1.3, at=midx)

dev.off()
fig3.2qijmaxse <- max(unlist(lapply(res, function(X){X$Dq_se})), na.rm = TRUE)
cat("maximum standard error in figure three (qij) is", fig3.2qijmaxse, "\n(M=", M, ")\n")
saveRDS(fig3.2qijmaxse, file="results_numeric/fig3.2qijmaxse.RDS")
#Sys.time()