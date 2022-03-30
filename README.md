# Code repository for the paper ``Relationships in the extremes and their influence on competition and coexistence''

Pimsupa Jasmin Albert, University of Kansas  
Daniel C. Reuman, University of Kansas  

## Introduction

This repository can be used to reproduce the complete analyses behind the paper "Relationships in the extremes and their 
influence on competition and coexistence" and to recompile the latex for the paper itself.

## How to run analyses and compile paper

To reproduce the analyses: 1) make your R working directory the ``code'' directory; 2) run MASTER.R. If all dependencies 
are in place (see next section) then all results supporting the paper will be saved to the results directories. This 
will take some time (a few days, on the first run, depending on your computer speed). 

To compile the paper latex, open paper/Paper.Rnw and paper/SupMat.Rnw and compile them using Sweave. If you are using
Rstudio and if it is configured to work with Sweave, this will just be a button in Rstudio entitled "Compile PDF".

## Dependencies

### Core dependencies

R, R studio configured to compile Sweave documents, latex and bibtex, checkpoint package in R. 

We used R version 3.6.3 running on 

The codes uses the 
R `checkpoint` package. This is set up in the MASTER.R file. The `checkpoint` package then automatically scans through 
other files looking for other required R packages. It then downloads and installs the newest versions of those packages 
available on the given date. This helps ensure that re-compiling the document uses the same code that was originally 
used. This can take some time on first run (10 minutes or so) but it is faster on subsequent runs because the packages 
are already installed. This also means that R package dependencies should only be the `checkpoint` package, since that 
package should scan for other packages and install them locally. 


To mention: this was done on ubuntu linux, you need X numbers of cores available




