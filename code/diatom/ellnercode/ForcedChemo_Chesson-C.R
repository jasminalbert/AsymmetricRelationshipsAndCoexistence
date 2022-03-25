##################################################################################### 
# Forced chemostat model: do the simulations to calculate various \bar{r}'s 
# and put them together into the various deltas, epsilons, and Deltas. 
#
# Here the partitioning uses E = temperature, C =  Chesson's C as in our 2016 paper 
# Following our dictum, we do not use q_ir.  
#
# This script does the analysis with a no-fluctuations simulation as the baseline, 
# rather than mean values from an invasion simulation with varying temperature.  
#
# WARNING: 
# For historical reasons, variables called delta here correspond to $\varepsilon$
# terms in the manuscript. We figure it's better to leave it that way, than to risk
# introducing errors by trying to make this consistent with the text. 
#
# Species 1 is Fragillaria, species 2 is Cyclotella 
#####################################################################################
#rm(list=ls(all=TRUE))
#graphics.off(); 

#root=ifelse(.Platform$OS.type=="windows","c:/repos","~/ipm"); # modify as needed
#setwd(paste0(root,"/drivers/Manuscripts/Partitioning/Chemostat")); # modify as needed 

require(deSolve); 
source("ForcedChemoSubs.R"); rm(r1); rm(r2); 
# species 1 is Fragillaria, species 2 is Cyclotella 

# q12= 1.099705;  q21= 0.8933319; # Chesson q_ir
q12 = q21 = 1; 
#######################################################################	
# Parameters. Units are: 
#          time in days, resource in mu-mol/L, populations in 10^6/L 
#######################################################################
parms=c(Tbar=18, a=6, P=60, D=0.09, S=35); 

############################################### 
# simulate with 1 as invader, 2 as resident 
###############################################
y0=c(R=.1,x1=0,x2=20); times=seq(0,3600,by=.1); 
out1i = ode(y0,times,func=forceChemo,parms=parms); 
e = times > max(times-1200); out1i = out1i[e,]; 
minus1 = data.frame(out1i); 
names(minus1)=c("time","R","x1","x2","temp"); 
    
#### repeat with no fluctuations     
parms2=c(Tbar=18, a=0, P=60, D=0.09, S=35);    
y0=c(R=.1,x1=0,x2=20); times=seq(0,3600,by=.1); 
out1i = ode(y0,times,func=forceChemo,parms=parms2); 
e = times > max(times-1200); out1i = out1i[e,]; 
minus1.noVar = data.frame(out1i); 
names(minus1.noVar)=c("time","R","x1","x2","temp"); 

############################################### 
# simulate with 2 as invader, 1 as resident 
###############################################
y0=c(R=.1,x1=20,x2=0); times=seq(0,3600,by=.1); 
out2i = ode(y0,times,func=forceChemo,parms=parms); 
e = times > max(times-1200); out2i = out2i[e,];
minus2 = data.frame(out2i); 
names(minus2)=c("time","R","x1","x2","temp"); 
    
#### repeat with no fluctuations     
parms2=c(Tbar=18, a=0, P=60, D=0.09, S=35);       
out2i = ode(y0,times,func=forceChemo,parms=parms2); 
e = times > max(times-1200); out2i = out2i[e,];
minus2.noVar = data.frame(out2i); 
names(minus2.noVar)=c("time","R","x1","x2","temp");     
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculations for 1 invading 2 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 
#################################################
# r's for 1 invading 2 
#################################################
E = minus1$temp;  E.noVar = minus1.noVar$temp; 
C1 = (K1fun(E)+minus1$R)/minus1$R 
C1.noVar = (K1fun(E.noVar) + minus1.noVar$R)/minus1.noVar$R   

C2 = (K2fun(E)+minus1$R)/minus1$R 
C2.noVar = (K2fun(E.noVar) + minus1.noVar$R)/minus1.noVar$R   

Estar = mean(E.noVar); C1star=mean(C1.noVar); C2star=mean(C2.noVar); 
E0 = mean(E); C10=mean(C1); C20=mean(C2); 

### "exact sharping" - all possible combinations 
EC1=expand.grid(E,C1); E1sharp=EC1[,1]; C1sharp=EC1[,2]; rm(EC1); 
EC2=expand.grid(E,C2); E2sharp=EC2[,1]; C2sharp=EC2[,2]; rm(EC2); 

r1star = r1C(Estar,C1star,parms2); # these are r1(E*,C*), no-variance baseline 
r2star = r2C(Estar,C2star,parms2); 

r10 = r1C(E0,C10,parms);  # r1(E0,C0), mean (E,C) in presence of varying temperature 
r20 = r2C(E0,C20,parms); 

r1sharp = mean(r1C(E1sharp,C1sharp,parms)) # these are rbar(E-sharp,C-sharp)
r2sharp = mean(r2C(E2sharp,C2sharp,parms)) 

r1bar = mean(r1C(E,C1,parms)); # these are rbar(E,C) 
r2bar = mean(r2C(E,C2,parms)); 

#################################################
# delta's for 1 invading 2 
#################################################
delta1.0 = r10-r1star; # this is epsilon_1^{prime} in the paper's notation 
delta2.0 = r20-r2star;  

delta1.sigmaE = mean(r1C(E,C10,parms))-r10; 
delta2.sigmaE = mean(r2C(E,C20,parms))-r20;

delta1.sigmaC = mean(r1C(E0,C1,parms))-r10; 
delta2.sigmaC = mean(r2C(E0,C2,parms))-r20;

epsilon1.EC = r1sharp - r10 - (delta1.sigmaE + delta1.sigmaC); 
epsilon2.EC = r2sharp - r20 - (delta2.sigmaE + delta2.sigmaC); 

delta1.rhoEC = r1bar - r1sharp;
delta2.rhoEC = r2bar - r2sharp; 

#################################################
# Delta's for 1 invading 2 
#################################################
Delta1.star = r1star - q12*r2star; 
Delta1.0 = delta1.0 - q12*delta2.0; 
Delta1.sigmaE = delta1.sigmaE - q12*delta2.sigmaE;
Delta1.sigmaC = delta1.sigmaC - q12*delta2.sigmaC; 
Delta1.rhoEC = delta1.rhoEC - q12*delta2.rhoEC; 
Delta1.epsilonEC = epsilon1.EC - q12*epsilon2.EC; 

delta1I = c(r1star,r10,delta1.0,delta1.sigmaE,delta1.sigmaC,epsilon1.EC,delta1.rhoEC); 
delta2R = c(r2star,r20,delta2.0,delta2.sigmaE,delta2.sigmaC,epsilon2.EC,delta2.rhoEC); 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculations for 2 invading 1 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#################################################
# r's for 2 invading 1 
#################################################
E = minus2$temp;  E.noVar = minus2.noVar$temp; 
C1 = (K1fun(E)+minus2$R)/minus2$R 
C1.noVar = (K1fun(E.noVar)+minus2.noVar$R)/minus2.noVar$R   

C2 = (K2fun(E)+minus2$R)/minus2$R 
C2.noVar = (K2fun(E.noVar)+minus2.noVar$R)/minus2.noVar$R   

Estar = mean(E.noVar); C1star=mean(C1.noVar); C2star=mean(C2.noVar); 
E0 = mean(E); C10=mean(C1); C20=mean(C2); 

### "exact sharping" - all possible combinations 
EC1=expand.grid(E,C1); E1sharp=EC1[,1]; C1sharp=EC1[,2]; 
EC2=expand.grid(E,C2); E2sharp=EC2[,1]; C2sharp=EC2[,2]; 

r1star = r1C(Estar,C1star,parms2); # these are r1(E*,C*), no-variance baseline 
r2star = r2C(Estar,C2star,parms2); 

r10 = r1C(E0,C10,parms);  # r1(E0,C0), mean (E,C) in presence of varying temperature 
r20 = r2C(E0,C20,parms); 

r1sharp = mean(r1C(E1sharp,C1sharp,parms)) # these are rbar(E-sharp,C-sharp)
r2sharp = mean(r2C(E2sharp,C2sharp,parms)) 

r1bar = mean(r1C(E,C1,parms)); # these are rbar(E,C) 
r2bar = mean(r2C(E,C2,parms)); 

#################################################
# delta's for 2 invading 1 
#################################################
delta1.0 = r10-r1star; # this is epsilon_1^{prime} in the paper's notation 
delta2.0 = r20-r2star;  

delta1.sigmaE = mean(r1C(E,C10,parms))-r10; 
delta2.sigmaE = mean(r2C(E,C20,parms))-r20;

delta1.sigmaC = mean(r1C(E0,C1,parms))-r10; 
delta2.sigmaC = mean(r2C(E0,C2,parms))-r20;

epsilon1.EC = r1sharp - r10 - (delta1.sigmaE + delta1.sigmaC); 
epsilon2.EC = r2sharp - r20 - (delta2.sigmaE + delta2.sigmaC); 

delta1.rhoEC = r1bar - r1sharp;
delta2.rhoEC = r2bar - r2sharp; 

#################################################
# Delta's for 2 invading 1 
#################################################
Delta2.star = r2star - q21*r1star; 
Delta2.0 = delta2.0 - q21*delta1.0; 
Delta2.sigmaE = delta2.sigmaE - q21*delta1.sigmaE;
Delta2.sigmaC = delta2.sigmaC - q21*delta1.sigmaC; 
Delta2.rhoEC = delta2.rhoEC - q21*delta1.rhoEC; 
Delta2.epsilonEC = epsilon2.EC - q21*epsilon1.EC; 

delta1I = c(r1star,r10,delta1.0,delta1.sigmaE,delta1.sigmaC,epsilon1.EC,delta1.rhoEC); 
delta2R = c(r2star,r20,delta2.0,delta2.sigmaE,delta2.sigmaC,epsilon2.EC,delta2.rhoEC);  
 
##################################################################### 
# Results
#####################################################################

# re-compute invasion growth rates of both species 
# Uses the nonChesson growth rate functions, which is OK  
r1.inv = mean(r1(minus1$temp,minus1$R,parms)); 
r2.inv = mean(r2(minus2$temp,minus2$R,parms)); 

cat("Fragillaria","\n"); 
u1 <- c(Delta1.star,Delta1.0,Delta1.sigmaE,Delta1.sigmaC,Delta1.rhoEC,Delta1.epsilonEC); 
names(u1)<-c("base","mean","sigmaE","sigmaC","rhoEC","epsilonEC"); 
cat(sum(u1),r1.inv,"should be equal", "\n");   
 
cat("Cyclotella","\n"); 
u2 <- c(Delta2.star,Delta2.0,Delta2.sigmaE,Delta2.sigmaC,Delta2.rhoEC,Delta2.epsilonEC); 
names(u2)<-c("base","mean","sigmaE","sigmaC","rhoEC","epsilonEC"); 

cat(sum(u2),r2.inv,"should be equal", "\n");  

# Stabilizing and equalizing components
Delta.bar = 0.5*(u1+u2); 
zeta.1 = u1 - Delta.bar
zeta.2 = u2 - Delta.bar; 


round(u1,digits=3); 
round(u2,digits=3); 
round(Delta.bar,digits=3); 
round(zeta.1,digits=3); 
round(zeta.2,digits=3); 





