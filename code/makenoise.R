source("./makenoise_fxn.R")

M <- 1000000
b_tilde <- makenoise(M)
u_tilde <- rnorm(M)

rho <- cor(b_tilde$l)[1,2]

cat("M=",M, " rho=", rho)

if(dir.exists("../results_numeric")==FALSE){
  dir.create("../results_numeric")
}

save(b_tilde,u_tilde, file = "../results_numeric/noise_etc.RData")
save(M, file="../results_numeric/M.rda")
save(rho, file="../results_numeric/rho.rda")