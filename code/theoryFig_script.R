# this script makes the panels for fig 2 in paper. calls plotting function from theory_panel_fxn.R for each panel

#library(latex2exp)

nx<-seq(from=-4,to=4,by=0.02)
ny<-seq(from=-4,to=4,by=0.02)

bx<-seq(from=0,to=1.0,by=0.01)
by<-seq(from=0,to=1.0,by=0.01) 

# A) BETA MARGS CLAYTON COP
# actual EC
sampsA <- theoryfigpanel(dbmarg, pbmarg, qbmarg, ccop, bx, by, method="gnrt_var",col=T,log_points=T)
graphics::text(x=lblc["E","x"], y=lblc["E","y"], labels=expression('E'['i']), cex=2)
graphics::text(x=lblc["C","x"], y=lblc["C","y"], labels=expression('C'['i\\i']), srt=90, cex=2)
#TeX version: labels = 
	# TeX('$E_i$')	# TeX('($C_{i\\setminus i}$)')
graphics::text(x=lblc["lab","x"], y=lblc["lab","y"],labels=paste0("(",letters[1],")"),cex=2)


# B) NORM MARGS CLAY COP
# EC marginals individially normalized
sampsB <- theoryfigpanel(dnmarg, pnmarg, qnmarg, ccop, nx, ny, samps=sampsA,oldp=pbmarg,xpd=NA,col=T, log_points=FALSE)
graphics::text(x=lblc["E","x"], y=lblc["E","y"], labels=c(expression(varphi^-1*('F'['E'['i']]("E"["i"])) )), cex=2)
graphics::text(x=lblc["C","x"]-.05, y=lblc["C","y"], labels=c(expression(varphi^-1*('F'['C'['i\\i']]("C"["i\\i"])) )), srt=90, cex=2)
#TeX('$\\varphi^{-1} ( F_{E_i}(E_i) )$')
#TeX('$\\varphi^{-1} ( F_{C_{i\\setminus i}}(C_{i\\setminus i}) )$')	
graphics::text(x=lblc["lab","x"], y=lblc["lab","y"],labels=paste0("(",letters[2],")"), cex=2)

# C) NORM MARGS NORM COP
# normal with rho same as actual EC
sampsC <- theoryfigpanel(dnmarg, pnmarg, qnmarg, ncop, nx, ny, method="gnrt_var", xpd=NA,col=T, log_points=FALSE)
graphics::text(x=lblc["E","x"], y=lblc["E","y"], labels=c(expression('e'['i'])), cex=2)
graphics::text(x=lblc["C","x"]+.06, y=lblc["C","y"], labels=c(expression('c'['i\\i'])), srt=90, cex=2)
#TeX('$e_i$')
#TeX('($c_{i\\setminus i}$)')
graphics::text(x=lblc["lab","x"], y=lblc["lab","y"],labels=paste0("(",letters[3],")"), cex=2)

# D) BETA MARGS NORM COP
#same marginals and rho as EC but association structure is symmetric (correlation per se)
sampsD<-theoryfigpanel(dbmarg, pbmarg, qbmarg, ncop, bx, by,samps=sampsC, oldp=pnmarg,col=T,log_points=T)
graphics::text(x=lblc["E","x"], y=lblc["E","y"], labels=c(expression('E'['i']^'||')), cex=2)
#=='F'['E'['i']]^-1*(varphi('e'['i'])) ))
graphics::text(x=lblc["C","x"], y=lblc["C","y"], labels=c(expression('C'['i\\i']^'||')), srt=90, cex=2)
#=='F'['C'['i\\i']]^-1*(varphi('c'['i/i'])) ))
#TeX('($E_i^{||} = F_{E_i}^{-1}\circ \varphi(e_i)$)')
#TeX(r'($C_{i\setminus i}^{||} = F_{C_{i\setminus i}}^{-1}\circ \varphi(c_{i\setminus i})$)')
graphics::text(x=lblc["lab","x"], y=lblc["lab","y"],labels=paste0("(",letters[4],")"), cex=2)

# E) BETA MARGS NO COV CLAY COP
# same marginals - no correlation (variance per se)
if(new){
graphics::plot.new();
sampsB2 <- theoryfigpanel(dbmarg, pbmarg, qbmarg, shcop, bx, by, samps=sampsA, method="var_perse",col=T,log_points=T)
graphics::text(x=lblc["E","x"], y=lblc["E","y"], labels=c(expression('E'['i']^'#')), cex=2)
graphics::text(x=lblc["C","x"], y=lblc["C","y"], labels=c(expression('C'['i\\i']^'#')), srt=90, cex=2)
graphics::text(x=lblc["lab","x"], y=lblc["lab","y"],labels=paste0("(",letters[5],")"), cex=2)
}

