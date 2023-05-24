#library(latex2exp)

nx<-seq(from=-4,to=4,by=0.02)
ny<-seq(from=-4,to=4,by=0.02)

bx<-seq(from=0,to=1.0,by=0.01)
by<-seq(from=0,to=1.0,by=0.01)

# A BETA MARGS CLAYTON COP
sampsA <- theoryfigpanel(dbmarg, pbmarg, qbmarg, ccop, bx, by, gnrt_var=TRUE,col=T,log_points=T)
text(x=lblc["E","x"], y=lblc["E","y"], labels=expression('E'['i']), cex=2)
text(x=lblc["C","x"], y=lblc["C","y"], labels=expression('C'['i\\i']), srt=90, cex=2)
#TeX version: labels = 
	# TeX('$E_i$')	# TeX('($C_{i\\setminus i}$)')
text(x=lblc["lab","x"], y=lblc["lab","y"],labels=paste0("(",letters[1],")"),cex=2)

# B NORM MARGS CLAY COP
sampsB <- theoryfigpanel(dnmarg, pnmarg, qnmarg, ccop, nx, ny, samps=sampsA, oldp=pbmarg,xpd=NA,col=T, log_points=FALSE)
text(x=lblc["E","x"], y=lblc["E","y"], labels=c(expression(varphi^-1*('F'['E'['i']]("E"["i"])) )), cex=2)
text(x=lblc["C","x"], y=lblc["C","y"], labels=c(expression(varphi^-1*('F'['C'['i\\i']]("C"["i\\i"])) )), srt=90, cex=2)
#TeX('$\\varphi^{-1} ( F_{E_i}(E_i) )$')
#TeX('$\\varphi^{-1} ( F_{C_{i\\setminus i}}(C_{i\\setminus i}) )$')
text(x=lblc["lab","x"], y=lblc["lab","y"],labels=paste0("(",letters[2],")"), cex=2)

# C NORM MARGS NORM COP
sampsC <- theoryfigpanel(dnmarg, pnmarg, qnmarg, ncop, nx, ny, gnrt_var=TRUE, xpd=NA,col=T, log_points=FALSE)
text(x=lblc["E","x"], y=lblc["E","y"], labels=c(expression('e'['i'])), cex=2)
text(x=lblc["C","x"], y=lblc["C","y"], labels=c(expression('c'['i\\i'])), srt=90, cex=2)
#TeX('$e_i$')
#TeX('($c_{i\\setminus i}$)')
text(x=lblc["lab","x"], y=lblc["lab","y"],labels=paste0("(",letters[3],")"), cex=2)

# D BETA MARGS NORM COP
sampsD<-theoryfigpanel(dbmarg, pbmarg, qbmarg, ncop, bx, by,samps=sampsC, oldp=pnmarg,col=T,log_points=T)
text(x=lblc["E","x"], y=lblc["E","y"], labels=c(expression('E'['i']^'||'=='F'['E'['i']]^-1*(varphi('e'['i'])) )), cex=2)
text(x=lblc["C","x"], y=lblc["C","y"], labels=c(expression('C'['i\\i']^'||'=='F'['C'['i\\i']]^-1*(varphi('c'['i/i'])) )), srt=90, cex=2)
#TeX('($E_i^{||} = F_{E_i}^{-1}\circ \varphi(e_i)$)')
#TeX(r'($C_{i\setminus i}^{||} = F_{C_{i\setminus i}}^{-1}\circ \varphi(c_{i\setminus i})$)')
text(x=lblc["lab","x"], y=lblc["lab","y"],labels=paste0("(",letters[4],")"), cex=2)



