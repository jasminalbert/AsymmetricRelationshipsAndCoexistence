source("./makenoise_fxn.R")

numeric_results_loc <- "../results_numeric"

if(dir.exists(numeric_results_loc)==FALSE){
  dir.create(numeric_results_loc)
}

noise_loc <- paste(numeric_results_loc, "/noise.RData", sep = "")
M_loc <- paste(numeric_results_loc, "/M.RDS", sep = "")
rho_loc <- paste(numeric_results_loc, "/rho.RDS", sep = "")

M <- 1000000
b_tilde <- makenoise(M)
u_tilde <- rnorm(M)

rho <- cor(b_tilde$l)[1,2]

cat("M=",M, " rho=", rho)


save(b_tilde,u_tilde, file = noise_loc)
save(M, file=M_loc)
save(rho, file=rho_loc)

