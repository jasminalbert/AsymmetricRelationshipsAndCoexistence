#set of 5 different functions - all used to make figure 2:
  #transform
  #popsim
  #co.periods
  #co.pPlot
  #simsPlot

### transform ###
# transform noise by sigma and mu as defined by lottery model
#ARGS:
  #b_tilde      list-3 sets bivariate noise 
  #sigma        common sd of bivariate noise sets
  #mudif        mu1-mu2
#OUT:
  #list with each element as the environmental noise of one of the two species, 3 different types, =length 6
transform <- function(b_tilde, sigma, mudif){
  #multiply by sigma, shift second noise by mudif    
  #left tail ATA
  b_l1 <- sigma*b_tilde$l[,1] 
  b_l2 <- sigma*b_tilde$l[,2] - mudif
  
  #right tail ATA
  b_r1 <- sigma*b_tilde$r[,1] 
  b_r2 <- sigma*b_tilde$r[,2] - mudif
  
  #symmetric
  b_s1 <- sigma*b_tilde$s[,1] 
  b_s2 <- sigma*b_tilde$s[,2] - mudif
  
  #return as list length 6  
  return(list(b_l1=b_l1, b_l2=b_l2, b_r1=b_r1, b_r2=b_r2, 
                b_s1=b_s1, b_s2=b_s2))
}

### popsim ###
#ARGS:
  #b1     species 1 environmental noise, vector length M
  #b2     species 2 environmental noise, vector length M
  #N      total abundance (N1 + N2), numeric
  #N1     initial abundance of species 1, numeric
  #delta  death rate (0-1), numeric
  #M      length of noise, numeric
#OUT:
  #data.frame 2xM+1 with columns N1 and N2 that are the abundance of sp1 and sp2 over time
popsim <- function(b1,b2,N,N1,delta,M){
  #Fecundity, B_i=exp(b_i)
  B1 <- exp(b1)
  B2 <- exp(b2)
  
  #species 2 initial abundance
  N2 <- N-N1
  
  #simulate competitive lottery model to get abundances at each time step
  for (t in 1:M){
    #total new juveniles is sum of fecundity times abundance for each species
    tot_new_juvs <- (B1[t]*N1[t]) + (B2[t]*N2[t])
    
    #lottery model equation
    N1[t+1] <- (1-delta)*N1[t] + delta*N*((B1[t]*N1[t])/tot_new_juvs)	
    N2[t+1] <- (1-delta)*N2[t] + delta*N*((B2[t]*N2[t])/tot_new_juvs)
  }
  return(data.frame(N1=N1, N2=N2))
}

### co.periods ###
#computes lengths of periods of coexisting community states: 
  #noticeable coexistence or single species dominance 
  #when community is in one statw, how long is that length of time? (before it changes to the other) 
  #"noticeable coexistence" = both species abundance above a dominance threshold 
#ARGS:
  #pop    M x 6 data.frame from popsim (see above)
  #dom    dominance threshold that defines noticable coexistence, fraction of total abundance, 0-1
  #N      total abundance
#OUT
  #list of three, co.periods (length of noticeable coexistence) for each noise type
co.periods <- function(pop, dom, N){
  N1 <- pop[,1]
  
  #get logical vector of noticeable coexistence
  #noticably coexisting => TRUE
  coexist <- N1 < dom*N & N1 > (1-dom)*N
  
  #get lengths of runs (runs of TRUE, runs of FALSE)
  RLE <- rle(coexist)
  
  #return RLE object: list with length element and value (T/F) element 
  return(RLE)
}

### co.pPlot ###
#function to plot pop simulation with highlighted grey areas of noticeable coexistence
#ARGS:
  #co.p   RLE object (list length 2 with $length and $value)
  #pop    data.frame of N1 and N2
  #start  first time step to plot
  #end    last time step to plot, **must be even number**
#OUT:
  #plot
co.pPlot <- function(co.p, pop, start, end,...){
  if(end %%2 != 0){
    stop("end must be set to an even number")
  }
  #total abundance
  N <- sum(pop[1,])
  #empty plot
  graphics::plot(0, xlab="", ylab="", ylim=c(0,N), xlim=c(start,end), 
       col="white",bty="l",...)
  
  #rectangle highlights
  #set up borders
  ybottom <- 0
  ytop <- sum(N)
    
  if (co.p$values[1] == FALSE){
    xleft <- sum(co.p$lengths[1]+1)
    c <- 1
  }
  c <- 0
  odds <- seq(1,end-1,2)
  evens <- seq(2,end,2)
  xleft <- 0
  #loop through to plot rectangles  
  for (p in 1:length(co.p$lengths)){
    xright <- sum(co.p$lengths[1:(odds[p]+c)])
    graphics::rect(xleft, ybottom, xright, ytop, col="grey90", border=NA)
    xleft <- sum(co.p$lengths[1:(evens[p]+c)]) + 1
  }
  
  #plot lines
  graphics::lines(pop$N1[start:end])
}

### simsPlot ###
#runs population simulation with and without ATA and
#computes noticeable coexistence and dominance metrics
#plots simulations using co.pPlot 
#ARGS:
	#noise_loc		location of noise .RData object in directory, 
	#				contains b_tilde and u_tilde
	#mudif			mean difference between species; mu1-mu2
	#delta			death rate
	#sigma			standard deviation of noise
	#start			first time step to plot	
	#end				last time step to plot, **must be even**
#OUT:
	#plot
	#list of metrics, length 5

simsPlot <- function(noise_loc, mudif, delta, sigma, start=1, end=500){
	#plotting titles
	titles <- list(main=c("without ATA", "with ATA"), lab=c("(a)","(b)"))
	
  	load(noise_loc)
  	#transform standard noise by mu and sigma
  	n <- transform(b_tilde, sigma=sigma, mudif=mudif)
  	#pop noise; only need left and symmetric; repackage for loop
  	pn <- list(s=list(n$b_s1, n$b_s2), l=list(n$b_l1, n$b_l2))
  	
  	coMets <- list(Sym=0, Asy=0)
  
  	for (n in 1:length(pn)){
  		#population simulation
  		pop <- popsim(pn[[n]][[1]], pn[[n]][[2]],N=50,N1=25,delta=delta,M)
  		#coexisting periods
  		cop <- co.periods(pop[start:end,], dom=.95, N=50)
  		#mean noticeable coexistence period length
  		comean <- round(mean(cop$length[cop$values==TRUE]),3)
  		#mean dominance period length
  		dommean <- round(mean(cop$length[cop$values==FALSE]),3) 
  		coMetssub <- list(co=comean, dom=dommean)	
  		coMets[[n]] <- coMetssub
  	
  		#plot
  		co.pPlot(cop, pop, start, end)
  		graphics::title(main=titles$main[[n]], line=-0.3, cex.main=1.8)   
  		graphics::mtext(titles$lab[[n]], side=3, line = -1, at=1000, cex=1.5)
  	}
	
	return(list(coMets=coMets, params = c(mu1=0, mu2=abs(mudif), sigma, delta)))

}


