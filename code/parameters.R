delta <- seq(0.1,1,0.1)
mu1 <- seq(.1,.9,.1) 
mu2 <- mu1 
sigma <- c(1.6, 3.2,6.4)

#mu1>mu2 
mu <- c(NA,NA)
for (i in 1:length(mu1)){
  for (j in 1:length(mu2)){
    if(mu1[i]<=mu2[j]){
      mu <- rbind(mu, cbind(mu1[i],mu2[j]))
    }
  }
}
mu <- mu[-1,]
#print(mu) #45x2

delta <- c(0.2, 0.4, 0.6, 0.8)
sigma <- c(2, 4, 5, 6)
mudif <- c(0.0, -0.3, -0.6, -0.8)

if(dir.exists("../results_numeric")==FALSE){
  dir.create("../results_numeric")
}


save(delta, mu, mudif, sigma, file = "../results_numeric/params.RData")
