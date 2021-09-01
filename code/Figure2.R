source("./simsPlot_fxn.R")




#### ---- IMPEDING ---- ####
#( r w/o ATA is positive and r with ATA is negative - ATA causes r to become negative )

# *** use this one for paper *** #
#source("./makenoise.R")
pdf("../results_figs/figure2.pdf") 
par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(2,1,1,1))
fig2 <- simsPlot(mudif=-0.8, delta=0.8, sigma=5, end=1000) # r_woATA = 0.098, r = -0.0373, [E||C] = -0.1354
dev.off()
save(fig2, file="../results_numeric/fig2data.RData")


png("../results_figs/figure2.PNG", width=1000, height=600, pointsize=15) 
par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(2,1,1,1))
simsPlot(mudif=-0.8, delta=0.8, sigma=5, end=1000) # r_woATA = 0.098, r = -0.0373, [E||C] = -0.1354
dev.off()

png("../results_figs/figure2.1.2.PNG", width=1000, height=600, pointsize=15)
par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(2,1,1,1))
simsPlot(mudif=-0.6, delta=0.8, sigma=4) # r_woATA = 0.0185, r = -0.0602, [E||C] = -0.0787
dev.off()

#### ---- FACILITATING ---- ####
#( r w/o ATA is negative and r with ATA is positive - ATA causes r to become positive )
pdf("../results_figs/figure2.2.pdf")
#png("figure2.2.1.PNG", width=1000, height=600, pointsize=15)
par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(2,1,1,1))
simsPlot(mudif=0, delta=1, sigma=7) # r_woATA = -0.0071, r = 0.00036, [E||C] = 0.0074
#dev.off()
#png("figure2.2.2.PNG", width=1000, height=600, pointsize=15)
#par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(2,1,1,1))
simsPlot(mudif=0, delta=1, sigma=6) # r_woATA = -0.0061, r = 0.00031, [E||C] = 0.0064
#dev.off()
#png("figure2.2.3.PNG", width=1000, height=600, pointsize=15)
#par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(2,1,1,1))
simsPlot(mudif=0, delta=1, sigma=5) # r_woATA = -0.0051, r = 0.00026, [E||C] = 0.0053
#dev.off()
#png("figure2.2.4.PNG", width=1000, height=600, pointsize=15)
#par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(2,1,1,1))
simsPlot(mudif=0, delta=1, sigma=3) # r_woATA = -0.0030, r = 0.00016, [E||C] = 0.0030
#dev.off()
#png("figure2.2.5.PNG", width=1000, height=600, pointsize=15)
#par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(2,1,1,1))
simsPlot(mudif=0, delta=1, sigma=2) # r_woATA = -0.0020, r = 0.00010, [E||C] = 0.0020
#dev.off()
#png("figure2.2.6.PNG", width=1000, height=600, pointsize=15)
#par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(2,1,1,1))
simsPlot(mudif=0, delta=1, sigma=1) # r_woATA = -0.0010, r = 0.00005, [E||C] = 0.0011

dev.off()

#most negative [E||C] - criteria for old sim figure
#also max r
png("../results_figs/figure2.3.1.PNG", width=1000, height=600, pointsize=15)
par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(2,1,1,1))
simsPlot(mudif=0, delta=0.5, sigma=7) # r_woATA = 1.131, r = 0.865, [E||C] = -0.266
dev.off()



#min r (negative)
png("../results_figs/figure2.4.1.PNG", width=1000, height=600, pointsize=15)
par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(2,1,1,1))
simsPlot(mudif=-0.8, delta=1, sigma=1)
dev.off()



