source("./decomposition_fxn.R")

#function to plot decomposition against mudif to find "differential coexistence"
  #the degree of difference in mu1-mu2 that can result in ATAs impeding or facilitating coexistence

#ARGS:
#mudif  a descending vector of mu1-mu2 to plot against and compute decompostion for
#sigma  an integer value of sigma
#delta  an integer value of delta

dePlot2 <- function(noise_loc, data_exist=TRUE, data_loc, mudif, sigma, delta, qij=FALSE, legend=FALSE,...){
  load(noise_loc)
  
  if (data_exist==FALSE){ #if data_loc =0, create the data
    cat("\ncomputing results with M=", length(u_tilde),"\nthis make take awhile...")
    store <- vector(mode='list', length=length(mudif))
    for (i in 1:length(mudif)){
      store[[i]] <- decompose(mudif[i],sigma,delta,b_tilde,u=u_tilde)
    }
    saveRDS(store, data_loc)
    cat("\n\nresults computed and saved as", data_loc)
  } else {
    cat("\nresults loaded from", data_loc)
    store <- readRDS(data_loc)
  }
  
  range <- range(unlist(lapply(store, function(X){(X$D[6:8])})))
  
  plot(0, xlab="", ylab="", ylim=range*1.2, xlim=c(max(mudif),min(mudif)), col="white",... )
  lines(range(mudif), rep(0,2), col="gray",lwd=0.7)
  
  for (i in 1:length(mudif)){
    if(i>=2){
      dat <- cbind(store[[i-1]]$D, store[[i]]$D)
      if(qij==TRUE){
        dat <- cbind(store[[i-1]]$Dq, store[[i]]$Dq)
      }
      lines(c(mudif[i-1], mudif[i]), dat[6,1:2], col="red")#[E||C]
      lines(c(mudif[i-1], mudif[i]), dat[7,1:2], col="orange",lwd=2)#r
      lines(c(mudif[i-1], mudif[i]), dat[8,1:2], col="navy")
      
      if(all(dat[8,]>0) & all(dat[7,]<0) ){ #impeding
        rect(mudif[i-1], range[1]*1.2, mudif[i], range[2]*1.2, col="hotpink2", border=NA, 
             density=40, lty=3)
      }
      if(all(dat[8,]<0) & all(dat[7,1]>0) ){ #facilitating 
        rect(mudif[i-1], range[1]*1.2, mudif[i], range[2]*1.2, col="darkgoldenrod2", border=NA,
             density=40, lty=3)
      }
      
    }
  }
  return(store)
}


#dePlot2(mudif_list[[1]], 6,1)
#dePlot2(mudif, 5,1)
#dePlot2(mudif, 4,1)
#dePlot2(mudif, 2,1)

#dePlot2(mudif, 2,0.8)
#dePlot2(mudif, 4,0.8)
#dePlot2(mudif, 5,0.8)
#dePlot2(mudif, 6,0.8)

#dePlot2(mudif, 2,0.5)
#dePlot2(mudif, 4,0.5)
#dePlot2(mudif, 5,0.5)
#dePlot2(mudif, 6,0.5)







#mudif = seq(-0.2,-0.4, length.out = 60)
#dePlot2(mudif, 4,0.95)
#dePlot2(mudif, 4,0.94)
#dePlot2(mudif, 4,0.935)
#dePlot2(mudif, 4,0.93)





