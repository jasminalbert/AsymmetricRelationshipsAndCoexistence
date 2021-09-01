#THIS IS A FUNCTION
source("./decomposition_fxn.R")

dePlot1 <- function(noise_loc, mudif, delta, qij=FALSE, legend=FALSE,...){
  load(noise_loc)
  store <- vector(mode='list', length=7)
  sigma <- seq(0,7,1)
  
  for (i in 1:length(sigma)){
    store[[i]] <- decompose(mudif,sigma[i],delta,b_tilde,u_tilde)
  }
  range <- range(unlist(lapply(store, function(X){(X$D)})))
  
  plot(0, xlab="", ylab="", ylim=range*1.1, xlim=c(0,7), col="white",...)
  lines(0:7, rep(0,8), col="gray",lwd=0.7)
  
  for (i in 1:length(sigma)){
    if(i>=2){
      dat <- cbind(store[[i-1]]$D, store[[i]]$D)
      if(qij==TRUE){
        dat <- cbind(store[[i-1]]$Dq, store[[i]]$Dq)
      }
      lines(c(sigma[i-1], sigma[i]), dat[1,1:2])#null
      lines(c(sigma[i-1], sigma[i]), dat[2,1:2], lty=2)#E
      lines(c(sigma[i-1], sigma[i]), dat[3,1:2], lty=4)#C
      lines(c(sigma[i-1], sigma[i]), dat[4,1:2], lty=3)#(E#C)
      lines(c(sigma[i-1], sigma[i]), dat[5,1:2], col="blue")#[EC]
      lines(c(sigma[i-1], sigma[i]), dat[6,1:2], col="red")#[E||C]
      lines(c(sigma[i-1], sigma[i]), dat[7,1:2], col="orange",lwd=2)#r
    }
  }
  
  if(legend==TRUE){
    legend("topleft", 
           legend=c(expression(Delta[0]), expression(Delta[E]),
                    expression(Delta[C]), expression(Delta[("E#C")]),
                    expression(Delta["[EC]"]), expression(Delta["[E||C]"])
                    ,expression(r)), 
           col = c("black","black","black","black","blue","red","orange"),
           lty = c(1,2,4,3,1,1,1), bty="n", cex=0.8)
  }
  
  
  return(store)
}


