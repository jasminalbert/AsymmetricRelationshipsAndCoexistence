#This script make figure 4 in main text:
#plotting IGR and IGR-ATA to find ATA exclusion and ATA rescue in trait space (plot against mudif)
#for different sigma and delta

##libraries used (invoked with ::): graphics, grDevices

### source function ###
#lottery model decomposition
source("./lotteryDecomp_fxns.R")

### define plotting function ###
### dePlot2 ###
# function to plot decomposition against mudif to find "differential coexistence" and ATA effects
# the degree of difference in mu1-mu2 that can result in ATAs impeding or facilitating coexistence
#ARGS:
  #noise_loc  location of noise .RData object in directory, 
  #				      contains b_tilde and u_tilde
  #res_exist  do .RDS numeric results for fig 4 exists? default to TRUE
  #res_loc    location of numeric results in directory
  #mudif      a descending vector of mu1-mu2 to plot against and compute decompostion for
  #sigma      an integer value of standard deviation
  #delta      an integer value of death rate
  #qij        use qij!=1? deafult to FALSE
#OUT:

dePlot2 <- function(noise_loc, res_exist=TRUE, res_loc, mudif, sigma, delta, qij=FALSE, legend=FALSE,...){
  load(noise_loc)
  
  if (res_exist==FALSE){ #create the data
    cat("\n\ncomputing results with M=", length(u_tilde),"\nthis make take awhile...")
    #empty list
    store <- vector(mode='list', length=length(mudif))
    
    #iterate across mudif to compute decomposition (see sourced script)
    for (i in 1:length(mudif)){
      store[[i]] <- decompose(mudif[i],sigma,delta,b_tilde,u=u_tilde)
    }
    saveRDS(store, res_loc)
    cat("...results computed and saved as", res_loc)
  } else {
    cat("\nresults loaded from", res_loc)
    store <- readRDS(res_loc)
  }
  #get range of Deltas for plotting
  range <- range(unlist(lapply(store, function(X){(X$D[6:8])})))
  
  #empty plot
  graphics::plot(0, xlab="", ylab="", ylim=range*1.2, xlim=c(max(mudif),min(mudif)), col="white",... )
  graphics::lines(range(mudif), rep(0,2), col="gray",lwd=0.7)#line at zero
  
  for (i in 1:length(mudif)){
    if(i>=2){ #use Deltas computed with qij=1
      dat <- cbind(store[[i-1]]$D, store[[i]]$D)
      if(qij==TRUE){ #qij!=1
        dat <- cbind(store[[i-1]]$Dq, store[[i]]$Dq)
      }
      graphics::lines(c(mudif[i-1], mudif[i]), dat[6,1:2], col="red")#[E||C]
      graphics::lines(c(mudif[i-1], mudif[i]), dat[7,1:2], col="orange",lwd=2)#r
      graphics::lines(c(mudif[i-1], mudif[i]), dat[8,1:2], col="navy") #growth without ATA
      
      #ATA effect shading 
      if(all(dat[8,]>0) & all(dat[7,]<0) ){ #impeding (ATA exclusion)
        graphics::rect(mudif[i-1], range[1]*1.2, mudif[i], range[2]*1.2, col="hotpink2", border=NA, 
             density=40, lty=3)
      }
      if(all(dat[8,]<0) & all(dat[7,1]>0) ){ #facilitating (ATA rescue)
        graphics::rect(mudif[i-1], range[1]*1.2, mudif[i], range[2]*1.2, col="darkgoldenrod2", border=NA,
             density=40, lty=3)
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
fig4_loc <- paste(fig_loc,"fig4.pdf",sep="")

### location to load required numerics from ###
params_loc <- paste(numeric_results_loc, "/params.RData", sep="")
noise_loc <- paste(numeric_results_loc, "/noise.RData", sep = "")
fig4numres_loc <- paste(numeric_results_loc, "/fig4numres_", sep="")

load(params_loc)
#set mudif
mudif_4 = seq(range(mudif)[2], range(mudif)[1]-0.5, length.out = 100)

### setting up dimensions ###
nd <- length(delta) 
ns <- length(sigma)
ph <- 3 #panel height
pw <- 2.5 #panel widht
ht <- c(1,rep(ph, ns)) #vector of panel heights
wd <- c(rep(pw, nd),2) #vector of panel widths
midx <- (pw*(nd/2))/sum(wd) #middle of figure panels; horizontal
midy <- (ph*(ns/2))/sum(ht) #middle of figure panels; vertical 

### make fig ###
grDevices::pdf(fig4_loc)
graphics::par(mgp=c(3,0.5,0), mar = c(0.5,1,1,1), oma=c(4,4,2,2), xpd=NA)
graphics::layout(matrix(c(2,3,4,5,1,
                6,7,8,9,1,	
                10,11,12,13,1,
                14,15,16,17,1,
                18,19,20,21,1), ncol=5, byrow=TRUE), 
       heights=ht, widths=wd)

graphics::plot.new() #1: legend panel
graphics::legend("topright", 
       legend=c(expression(GWR),expression(GWR-Delta[i]^"[EC]"), expression(Delta[i]^"[EC]")),
       col = c("orange","navy", "red"),
       lty = c(1,1,1), bty="n", cex=1.8, inset=c(-0.51,-0.05),
       y.intersp = 1.1, x.intersp = 0.05, seg.len=0.8, lwd=1.5)
graphics::legend("topright", legend=c("ATA \nrescue", "ATA \nexclusion"), density=50, fill=c("darkgoldenrod2","hotpink2"),
       bty="n", border="white", inset=c(-0.33, 0.09), cex=1.6, y.intersp = 1.6,
       x.intersp = 0.3)
#2-5: empty space panels
for (d in 1:nd){
  graphics::plot.new()
}
#6-21: plots
res <- vector(mode='list',length=1)
m <- 1
for (s in 1:ns){ #iterate through sigmas
  for (d in 1:nd){ #iterate through deltas
    
    #fig 4 numeric results location
    res_loc <- paste(fig4numres_loc,m,".RDS", sep="")
    
    #plot and get results
    res <- append(res, list(dePlot2(noise_loc, res_exist=ifelse(file.exists(res_loc), yes=TRUE, no=FALSE), res_loc=res_loc, mudif_4,sigma[s], delta[d], xaxt="n", cex.axis=1.3)))
    graphics::axis(1, labels=ifelse(m>12, yes=TRUE, no=FALSE), tick=TRUE, cex.axis=1.2)
    graphics::mtext(paste0("(", letters[m],")"), side=3, line=-1.7, at=-4.6, adj=1, cex=1.3)
    
    if(m<5){
      graphics::mtext(paste(delta[d]), side=3, line=0.2, font=2, cex=1)
    }
    if(m%%4==0){
      graphics::mtext(paste(sigma[s]), side=4, line=0.6, font=2, cex=1.1)
    }
    if(m==12){
      graphics::abline(v=-0.8, col="red", lty=2, xpd=T)
    }
    m <- m+1
  }
}
graphics::mtext("contribution to coexistence", side=2, outer=TRUE, line=1.3, font=2, cex=1.3, at=midy-0.005)
graphics::mtext(expression(sigma), side=4, outer=TRUE, line=-5.5, font=2, cex=2, at=midy)
graphics::mtext(expression(delta), side=3, outer=TRUE, line=-2.5, font=2, cex=2, at=midx)
graphics::mtext(expression(mu[1]-mu[2]), outer=TRUE, side=1, line=2, cex=1.5, at=midx)

grDevices::dev.off() #finish plotting

### get and save standard error ###
fig4maxse <- max(unlist(lapply(res[-1], function(Y){lapply(Y, function(X){X$D_se})})), na.rm=TRUE)
cat("\nmaximum standard error in figure four is", fig4maxse, "\n(M=", M, ")\n")
fig4maxse_loc <- paste(numeric_results_loc, "/fig4maxse.RDS", sep="")
saveRDS(fig4maxse, file=fig4maxse_loc)

