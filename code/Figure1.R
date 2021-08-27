load("results_numeric/noise_etc.RData")


P <- c(cor(b_tilde$l, method='pearson')[1,2], cor(b_tilde$s, method='pearson')[1,2], cor(b_tilde$r, method='pearson')[1,2])
S <- c(cor(b_tilde$l, method='spearman')[1,2],cor(b_tilde$s, method='spearman')[1,2],cor(b_tilde$r, method='spearman')[1,2])
PL <- round(P[1],4)

#col2rgb("darkblue")
blue <- rgb(0,0,0.545,0.3)

pdf("results_figs/fig1.pdf", width=7, height=3)
#png("fig1.png", height=560, width=1600,pointsize=22)
par(mar=c(0.5,0.5,0,0), oma=c(5,5,2,2))
layout(matrix(c(4,10,10,6,11,11,8,12,
                1,5 ,13,2,7 ,14,3,9 ),nrow=2,byrow=T),
       heights=c(0.25,1), widths=c(1,0.25,0.25,1,0.25,0.25,1,0.25))
#1-3
plot(b_tilde$l[1:1000,1],b_tilde$l[1:1000,2], col=blue, ylab=NA, xlab=NA, xlim=c(-4,4), ylim=c(-4,4), pch=16, cex=1.5)
#title(main='Left tail association', line=-1.2, cex.main=1.5)
text(x=4,y=-4,labels=bquote(~ rho == .(round(P[1],5))), adj=1)
mtext("A", side=3, line=-1, at=-4)

plot(b_tilde$s[1:1000,1],b_tilde$s[1:1000,2], col=blue, xlab=NA, ylab=NA, xlim=c(-4,4), ylim=c(-4,4), pch=16, cex=1.5)
#title(main='Symmetric association', line=-1.2, cex.main=1.5)
text(x=4,y=-4,labels=bquote(~ rho == .(round(P[2],5))), adj=1)

plot(b_tilde$r[1:1000,1],b_tilde$r[1:1000,2], col=blue, xlab=NA, ylab=NA, xlim=c(-4,4), ylim=c(-4,4), pch=16, cex=1.5)
#title(main='Right tail association', line=-1.2, cex.main=1.5)
text(x=4,y=-4,labels=bquote(~ rho == .(round(P[3],5))), adj=1)



title(xlab="variable 1", outer=TRUE, line=3, cex.lab=1.5)
title(ylab="variable 2", outer=TRUE, line=3, cex.lab=1.5)
#mtext(expression(rho==0.8258437), 1, outer=TRUE)

#marginals
#4-9
d <- density(b_tilde$l[,1])
plot(d$x, d$y, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n',type='l', xaxt='n', yaxt='n')
d <- density(b_tilde$l[,2])
plot(d$y, d$x, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n', type='l',xaxt='n', yaxt='n')

d <- density(b_tilde$s[,1])
plot(d$x, d$y, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n',type='l', xaxt='n', yaxt='n')
d <- density(b_tilde$s[,2])
plot(d$y, d$x, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n', type='l',xaxt='n', yaxt='n')

d <- density(b_tilde$r[,1])
plot(d$x, d$y, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n',type='l', xaxt='n', yaxt='n')
d <- density(b_tilde$r[,2])
plot(d$y, d$x, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n', type='l',xaxt='n', yaxt='n')

#title(main='Left-tailed Symmetric Right-tailed', outer=TRUE)

dev.off()


pdf("results_figs/fig1_vertical.pdf", width=2.8, height=7)

par(mar=c(0.5,0.5,0,0), oma=c(5,5,2,2))

layout(matrix(c(4,10,
                1,5 ,
                13,13,
                6,11,
                2,7,
                14,14,
                8,12,
                3,9),nrow=8,byrow=T),
       heights=c(0.25,1,0.2,0.25,1,0.2,0.25,1), widths=c(1,0.25))

#1-3
plot(b_tilde$l[1:800,1],b_tilde$l[1:800,2], col=blue, ylab=NA, xlab=NA, xlim=c(-4,4), ylim=c(-4,4), pch=16, cex=1.5)
text(x=4,y=-3.8,labels=bquote(~ rho == .(round(P[1],5))), adj=1)
mtext("A", side=3, line=-1.5, at=-3.7)

plot(b_tilde$s[1:800,1],b_tilde$s[1:800,2], col=blue, xlab=NA, ylab=NA, xlim=c(-4,4), ylim=c(-4,4), pch=16, cex=1.5)
text(x=4,y=-3.8,labels=bquote(~ rho == .(round(P[2],5))), adj=1)
mtext("B", side=3, line=-1.5, at=-3.7)

plot(b_tilde$r[1:800,1],b_tilde$r[1:800,2], col=blue, xlab=NA, ylab=NA, xlim=c(-4,4), ylim=c(-4,4), pch=16, cex=1.5)
text(x=4,y=-3.8,labels=bquote(~ rho == .(round(P[3],5))), adj=1)
mtext("C", side=3, line=-1.5, at=-3.7)

mtext("variable 1", outer=TRUE, side=1, line=2.5, cex.lab=1.25, at=0.43)
mtext("variable 2", outer=TRUE, side=2, line=2.5, cex.lab=1.25, at=0.48)

#marginals
#4-9
d <- density(b_tilde$l[,1])
plot(d$x, d$y, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n',type='l', xaxt='n', yaxt='n')
d <- density(b_tilde$l[,2])
plot(d$y, d$x, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n', type='l',xaxt='n', yaxt='n')

d <- density(b_tilde$s[,1])
plot(d$x, d$y, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n',type='l', xaxt='n', yaxt='n')
d <- density(b_tilde$s[,2])
plot(d$y, d$x, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n', type='l',xaxt='n', yaxt='n')

d <- density(b_tilde$r[,1])
plot(d$x, d$y, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n',type='l', xaxt='n', yaxt='n')
d <- density(b_tilde$r[,2])
plot(d$y, d$x, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n', type='l',xaxt='n', yaxt='n')

dev.off()











