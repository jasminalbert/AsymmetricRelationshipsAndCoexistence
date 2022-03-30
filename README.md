# Code repository for the paper ``Relationships in the extremes and their influence on competition and coexistence''

Pimsupa Jasmin Albert, University of Kansas  
Daniel C. Reuman, University of Kansas  

## Introduction

This repository can be used to reproduce the complete analyses behind the paper "Relationships in the extremes and their 
influence on competition and coexistence" and to recompile the paper itself.

## How to run analyses and compile paper

To reproduce the analyses: 1) make your R working directory the ``code'' directory; 2) run MASTER.R. If all dependencies 
are in place (see next section) then all results supporting the paper will be saved to one of the results directories. 
This will take some time (a few days, on the first run, depending on your computer speed). 

To compile the paper latex, open paper/Paper.Rnw and paper/SupMat.Rnw and compile them using Sweave. If you are using
Rstudio and if it is configured to work with Sweave, this will just be a button in Rstudio entitled "Compile PDF".

## Dependencies

### Core dependencies



To mention: this was done on ubuntu linux, you need X numbers of cores available




