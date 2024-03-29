################################################### 
# Temperature-dependent growth rate data
# Species 1 is Fragillaria, species 2 is Cyclotella 
################################################### 

Tvals <-  c(6,12,18,24); 
V1data <-  c(0.25,0.4,0.37,0.42);  # Fragiliaria data
V2data <-  c(0.32,0.36,0.39,0.02); # Cyclotella data 

K1data <-  c(0.06,0.30,0.20,0.25)
K2data <-  c(0.07,0.13,0.15,NA)

Q1data <-  c(3.21,1.44,2.03,2.02)
Q2data <-  c(0.301,0.862,0.619,NA)

############### Fit the data by connect-the-dots 
V1fun <-  approxfun(Tvals,V1data,rule = 2); 
V2fun <-  approxfun(Tvals,V2data,rule = 2);

K1fun <-  approxfun(Tvals,K1data,rule = 2); 
K2fun <-  approxfun(Tvals,K2data,rule = 2); 

Q1fun <-  approxfun(Tvals,Q1data,rule = 2); 
Q2fun <-  approxfun(Tvals,Q2data,rule = 2); 

############################################## 
# Break some fits away from connect-the-dots
##############################################
fitV1 <-  lm(V1data~Tvals+I(Tvals^2)); c1 <-  coef(fitV1); 
V1fun <-  function(temp) c1[1] + c1[2]*temp + c1[3]*temp^2

K1flat <-  c(0.06,0.25,0.25,0.25); 
K1fun<- approxfun(Tvals,K1flat,rule = 2); 

K2flat <-  c(0.07,0.14,0.14,0.14); 
K2fun<- approxfun(Tvals,K2flat,rule = 2);

V2dataMod <-  c(0.32,0.36,0.39,0.06)   # Cyclotella data with last value modified 
V2fun <-  approxfun(Tvals,V2dataMod,rule = 2); 

######################################
# Model 
######################################
forceChemo <- function(t,y,parms) {
	temp <- parms["Tbar"] + parms["a"]*sin(2*pi*t/parms["P"]); 
	V1 <- V1fun(temp); V2<- V2fun(temp); K1<- K1fun(temp); K2 <-  K2fun(temp) 
	Q1 <- Q1fun(temp); Q2 <-  Q2fun(temp);
	R <- y[1]; x1<- y[2]; x2<- y[3]; 
	up1 <- V1*R/(K1 + R); 
	up2 <- V2*R/(K2 + R)  
	dR <- parms["D"]*(parms["S"]-R) - Q1*x1*up1 - Q2*x2*up2; 
 	D <- parms["D"]; names(D) <- NULL;
	dx1 <- x1*(up1-D); 
	dx2 <- x2*(up2-D); 
	names(V1) <- NULL; names(up1) <- NULL; names(up2) <- NULL; names(R) <- NULL; 
	return( list(dx = c(dR,dx1,dx2), temp = temp) ) 
}	

r1 <-  function(temp,R,parms) {
    V1 <- V1fun(temp); K1<- K1fun(temp); 
 	up1 <- V1*R/(K1 + R); 
    return(up1 - parms["D"]); 
}

r2 <- function(temp,R,parms) {
    V2 <- V2fun(temp); K2<- K2fun(temp); 
 	up2 <- V2*R/(K2 + R); 
    return(up2 - parms["D"]); 
}

########################################################################### 
# Functions for doing our partition using Chesson's C where  r = V(E)/C - D 
###########################################################################
r1C <-  function(E,C,parms) {
    V1 <- V1fun(E); up1 <- V1/C;  
    return(up1 - parms["D"]); 
}

r2C <- function(E,C,parms) {
    V2 <- V2fun(E); up2 <- V2/C; 
    return(up2 - parms["D"]); 
}


######################################################################### 
# Function for doing our partition using Chesson's E and C 
# where  r = E/C - D 
#########################################################################
rEC <-  function(E,C,parms) {E/C - parms["D"]} 

