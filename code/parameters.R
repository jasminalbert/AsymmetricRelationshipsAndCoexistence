numeric_results_loc <- "../results_numeric"

if(dir.exists(numeric_results_loc)==FALSE){
  dir.create(numeric_results_loc)
}

params_loc <- paste(numeric_results_loc, "/params.RData", sep="")

delta <- c(0.2, 0.4, 0.6, 0.8)
sigma <- c(2, 4, 5, 6)
mudif <- c(0, -0.5, -2, -4)


save(delta, mudif, sigma, file = params_loc)
