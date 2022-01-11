#function to plot simulations side by side to compare effect of ATA contributions
source("./pop_sim.R")
source("./transform2.R")
source("./co.periods.R")
#end must be an even number

co.pPlot <- function(co.p, pop, start, end,...){
  
  if(end %%2 != 0){
    stop("end must be set to an even number")
  }
  plot(0, xlab="", ylab="", ylim=c(0,50), xlim=c(start,end), 
       col="white",bty="l",...)
  ybottom <- 0
  ytop <- 50
    
  if (co.p$values[1] == FALSE){
    xleft <- sum(co.p$lengths[1]+1)
    c <- 1
  }
  c <- 0
  odds <- seq(1,end-1,2)
  evens <- seq(2,end,2)
  xleft <- 0
    
  for (p in 1:length(co.p$lengths)){
      
    xright <- sum(co.p$lengths[1:(odds[p]+c)])
      
    rect(xleft, ybottom, xright, ytop, col="grey90", border=NA)
      
    xleft <- sum(co.p$lengths[1:(evens[p]+c)]) + 1
     
  }
  lines(pop$N1[start:end])
}

simsPlot <- function(noise_loc, mudif, delta, sigma, start=1, end=500){
  load(noise_loc)
  #transform noise
  n <- transform(b_tilde, u_tilde, rho, sigma=sigma, mudif=mudif, b_s=T)
  
  #asymmetric pop
  popA <- popsim(n$b_l1,n$b_l2,N=50,N1=25,delta=delta,M)
  copA <- co.periods(popA[start:end,], dom=0.95, N=50)
  #symmetric pop
  popS <- popsim(n$b_s1,n$b_s2,N=50,N1=25,delta=delta,M)
  copS <- co.periods(popS[start:end,], dom=0.95, N=50)
  
  comeanA <- round(mean(copA$length[copA$values==TRUE]),3)
  comeanS <- round(mean(copS$length[copS$values==TRUE]),3)
  dommeanA <- round(mean(copA$length[copA$values==FALSE]),3)
  dommeanS <- round(mean(copS$length[copS$values==FALSE]),3)
    
  co.pPlot(copS, popS, start, end)
  title(main="without ATA", line=0)   
  mtext("(a)", side=3, line = -.85, at=-16)
  
  co.pPlot(copA, popA, start, end)
  title(main="with ATA", line=0)
  mtext("(b)", side=3, line = -.85, at=-16)
  
  title(xlab="time", ylab="population of species 1", outer=T, line=-0.5)
  
  return(list(coSym = comeanS, domSym = dommeanS, coATA = comeanA, domATA = dommeanA, params = c(mudif, sigma, delta)))

}

#simsPlot(0,0.5,6)

