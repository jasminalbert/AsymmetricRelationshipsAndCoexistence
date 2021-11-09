
#data
Tvals <-  c(6,12,18,24); 
V1data <-  c(0.25,0.4,0.37,0.42);  # Fragiliaria data
V2data <-  c(0.32,0.36,0.39,0.02); # Cyclotella data 

K1data <-  c(0.06,0.30,0.20,0.25)
K2data <-  c(0.07,0.13,0.15,NA)

Q1data <-  c(3.21,1.44,2.03,2.02)
Q2data <-  c(0.301,0.862,0.619,NA)

#functions
V1fun <-  approxfun(Tvals,V1data,rule = 2) 
V2fun <-  approxfun(Tvals,V2data,rule = 2)

K1fun <-  approxfun(Tvals,K1data,rule = 2) 
K2fun <-  approxfun(Tvals,K2data,rule = 2) 

Q1fun <-  approxfun(Tvals,Q1data,rule = 2)
Q2fun <-  approxfun(Tvals,Q2data,rule = 2)

#modifications
fitV1 <-  lm(V1data~Tvals+I(Tvals^2)); c1 <-  coef(fitV1); 
V1quad <-  function(temp) c1[1] + c1[2]*temp + c1[3]*temp^2

K1flat <-  c(0.06,0.25,0.25,0.25); 
K1flatfun<- approxfun(Tvals,K1flat,rule = 2); 

K2flat <-  c(0.07,0.14,0.14,0.14); 
K2flatfun<- approxfun(Tvals,K2flat,rule = 2);

V2dataMod <-  c(0.32,0.36,0.39,0.06)   # Cyclotella data with last value modified 
V2modfun <-  approxfun(Tvals,V2dataMod,rule = 2)