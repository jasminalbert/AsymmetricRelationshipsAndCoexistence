#This script makes figure 3 in main text: 
#plotting decomposed mechanisms contributions against sigma
#for different deltas and mu1-mu2

##libraries used (invoked with ::): graphics, grDevices

### source function ###
#lottery model decomposition
source("./plotMechanisms.R")

### define plotting function ###
#### dePlot1 ###
# description
#ARGS:
  #noise_loc    location of noise .RData object in directory, 
  #				      contains b_tilde and u_tilde
  #sigma        *vector* of sigmas to use (x-axis)
  #mudif        mean difference btwn species, mu1-mu2, length 1
  #delta        death rate, length 1
  #qij          use scaling factor != 1? default FALSE
  #cols         colors of lines (7)
  #ltys         types of lines (7)
  #lwds         widths of lines (7)
#OUT:
  #plot with 7 lines for each mechanism/term and 
  #list of decomposition data.frame for each sigma value 
dePlot1 <- function(noise_loc, sigma, mudif, delta, cols, ltys,lwds,qij=FALSE,...){
  load(noise_loc)
  
  #empty list to store decomposition data.frames
  store <- vector(mode='list', length=length(sigma))
  
  #iterate across sigmas to compute decomposition (see sourced script)
  for (i in 1:length(sigma)){
    store[[i]] <- decompose(mudif,sigma[i],delta,b_tilde,u_tilde)
  }
  #get range of Deltas for plotting
  range <- range(unlist(lapply(store, function(X){(X$D)})))
  
  #empty plot
  graphics::plot(0, xlab="", ylab="", ylim=range*1.1, xlim=c(0,7), col="white",...)
  #line at zero
  graphics::lines(range(sigma), rep(0,2), col="gray",lwd=0.7)
  
  #lines for each mechanism/ Delta term
  for (i in 1:length(sigma)){
    #use qij=1
    if(i>=2){
      dat <- cbind(store[[i-1]]$D, store[[i]]$D)
      #or qij!=1
      if(qij==TRUE){
        dat <- cbind(store[[i-1]]$Dq, store[[i]]$Dq)
      }
      for(lin in 1:7){ #lines for each mechanism/term
        graphics::lines(c(sigma[i-1], sigma[i]), dat[lin, 1:2], 
              col=cols[lin], lty=ltys[lin], lwd=lwds[lin])
      }
    }
  }
  return(store)
}


### location to save results ###
numeric_results_loc <- "../results_numeric"
if(dir.exists(numeric_results_loc)==FALSE){
  dir.create(numeric_results_loc)
}
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
fig3_loc <- paste(fig_loc,"decompFigLLN.pdf",sep="")
fig3qij_loc <- paste(fig_loc,"fig3_qij.pdf",sep="")

### files to load ###
noise_loc <- paste0(numeric_results_loc, "/noise.RData")
params_loc <- paste0(numeric_results_loc, "/params.RData")
res_loc <- paste0(numeric_results_loc,"/fig3dat/")
if(dir.exists(res_loc)==FALSE){
  dir.create(res_loc)
}
fig3_dat_loc <- paste0(res_loc,"Deltas_d")

load(params_loc)
sigma <- seq(0,7,len=50)
delta <- delta[-4]
### setting dimensions ###
nd <- length(delta) 
nm <- length(mudif)
ph <- 3 #panel height
pw <- 2.5 #panel widht
ht <- c(1,rep(ph, nm)) #vector of panel heights
wd <- c(rep(pw, nd),2) #vector of panel widths
midx <- (pw*(nd/2))/sum(wd) #middle of figure panels; horizontal
midy <- (ph*(nm/2))/sum(ht) #middle of figure panels; vertical 

### setting up legend ###
terms <- c(expression(Delta[i]^0), expression(Delta[i]^E),
           expression(Delta[i]^C), expression(Delta[i]^"(E#C)"),
           expression(Delta[i]^"[E||C]"), expression(Delta[i]^"[EC]")
           ,expression(GWR))
cols <- c("black","black","black","black","blue","red","orange")
ltys <- c(1,2,4,3,1,1,1)
lwds <- c(rep(1,6),2)

### make fig ###
#qij=1
grDevices::pdf(fig3_loc, height=8, width=8)
graphics::par(mfrow=c(4,3), oma=c(3,4,4,7), mar=c(1,1,1,1), mgp=c(3,1,0), xpd=NA)
ylims <- list(c(-0.5,1.25),c(-0.5,1.25),c(-2,1.5), c(-4,2.5))
#graphics::legend("topright", legend=terms, col = cols, lty = ltys, bty="n", cex=1.8, inset=c(0,-0.05), y.intersp = 1.1, x.intersp = 0.1, seg.len=0.8, lwd=1.5)
n <- 1
for (m in seq_along(mudif)){ #iterates across mudif values
  for (d in seq_along(delta)){ #iterates across delta values
  	
  	file <- paste0(fig3_dat_loc,delta[d],"_md",mudif[m],".RDS")
  	#make files?
  	if(file.exists(file)==FALSE){
  		Delts <- makeDeltas(noise_loc, sigma, mudif[m], delta[d], qij=FALSE)
  		saveRDS(Delts, file=file)
  	}
  	#plot
  	plotco(file,sigma,ylim=ylims[[m]])
  	if (d==1){
  		graphics::axis(2, cex.axis=1.8, tck=-0.035, lwd.ticks=2)
  	}
  	if (n>9){
  		graphics::axis(1, cex.axis=1.8,tck=-0.028, lwd.ticks=2)
  	}
    graphics::mtext(paste0("(", letters[n],")"), side=3, line=-1.7, at=0, cex=1.3, adj=0)
    
    if (n==2){
    	graphics::mtext("adult death rate, ", line=2.5, font=2, cex=1.5, col="gray40", at=-0.75, adj=0)
    	graphics::mtext(expression(delta), line=2.5, font=2, cex=2.3, col="gray40", side=3, at=8, adj=1)
    }
    if (n<4){
      graphics::mtext(paste(delta[d]), side=3, line=0.2, font=2, cex=1.25, col="gray30")
    }
    if (n%%3==0){
      graphics::text(paste(mudif[m]), srt=-90, x=7.8,y=mean(ylims[[m]]),  font=2, cex=2, col="gray30")
    }
    if (n==6){
    	graphics::text("mean log fecundity difference, ", x=10,y=2.4, font=2, cex=2.5, srt=-90, adj=0, col="gray40")
    	graphics::text(expression(mu[1]-mu[2]), x=10,y=-4, font=2, cex=3, srt=-90, adj=1, col="gray40")
    }
    if (n==11){
    	graphics::mtext("log fecundity standard deviation, ",side=1, line=3, cex=1.25, at=-2.75, adj=0,col="gray40")
    	graphics::mtext(expression(sigma),side=1, line=3, cex=1.5, at=9.75, adj=1,col="gray40")
    }
    n <- n+1
  }
} 
graphics::mtext("contribution to coexistence", side=2, outer=TRUE, line=2, font=2, cex=1.7, col="gray40")

grDevices::dev.off() #finish plotting

### get standard error and save ###
#fig3maxse <- max(unlist(lapply(res, function(X){X$D_se})), na.rm = TRUE)
#cat("maximum standard error in figure three is", fig3maxse, "\n(M=", M, ")\n")
#fig3maxse_loc <- paste(numeric_results_loc, "/fig3maxse.RDS", sep="")
#saveRDS(fig3maxse, file=fig3maxse_loc)



### make fig ###
#qij!=1; same as above but qij=TRUE in plotting function
#grDevices::pdf(fig3qij_loc)
#graphics::par(mgp=c(3,0.6,0), mar = c(0.5,1,1,1), oma=c(4,4,2,2))
#graphics::layout(matrix(c(2,3,4,5,1,
#                6,7,8,9,1,	
#                10,11,12,13,1,
#                14,15,16,17,1,
#                18,19,20,21,1), ncol=5, byrow=TRUE), 
#       heights=ht, widths=wd)
#graphics::plot.new() #1
#graphics::legend("topright", legend=terms, col = cols, lty = ltys, bty="n", cex=1.8, 
#       inset=c(0,-0.05),y.intersp = 1.1, x.intersp = 0.1, #seg.len=0.8, lwd=1.5, xpd=NA)
#2-5
#for (d in 1:nd){
#  graphics::plot.new()
#}
#6-21
#res <- vector(mode='list',length=1)
#n <- 1
#for (m in 1:nm){
#  for (d in 1:nd){
#    res <- append(res,dePlot1(noise_loc, sigma, mudif[m], #delta[d], cols, ltys, lwds, qij=TRUE, xaxt="n", cex.axis=1.2))
#    graphics::axis(1, labels=ifelse(n>12, yes=TRUE, no=FALSE), tick=TRUE, cex.axis=1.3)
#    graphics::mtext(paste0("(", letters[n],")"), side=3, line=-1.7, at=0, cex=1.3, adj=0)
    
#    if(n<5){
#      graphics::mtext(paste(delta[d]), side=3, line=0.2, font=2, cex=1)
#    }
#    if(n%%4==0){
#      graphics::mtext(paste(mudif[m]), side=4, line=0.6, #font=2, cex=1)
#    }
#    n <- n+1
#  }
#}
#graphics::mtext("contribution to coexistence", side=2, outer=TRUE, line=1.3, font=2, cex=1.3, at=midy-0.005)
#graphics::mtext(expression(mu[1]-mu[2]), side=4, outer=TRUE, line=-5.5, font=2, cex=2, at=midy)
#graphics::mtext(expression(delta), side=3, outer=TRUE, line=-2.5, font=2, cex=2, at=midx)
#graphics::mtext(expression(sigma), outer=TRUE, side=1, line=2, cex=1.5, at=midx)

#grDevices::dev.off() #finish plotting

### get standard error and save ###
#fig3qijmaxse <- max(unlist(lapply(res, function(X){X$Dq_se})), na.rm = TRUE)
#cat("maximum standard error in figure three (qij) is", fig3qijmaxse, "\n(M=", M, ")\n")
#fig3qijmaxse_loc <- paste(numeric_results_loc, "/fig3qijmaxse.RDS", sep="")
#saveRDS(fig3qijmaxse, file=fig3qijmaxse_loc)
