# A BETA MARGS CLAYTON COP
sampsA <- theoryfigpanel(dbmarg, pbmarg, qbmarg, ccop, x, y, gnrt_var=TRUE)
text(x=lblc["E","x"], y=lblc["E","y"], labels=c(expression('E'['i'])), cex=2)
text(x=lblc["C","x"], y=lblc["C","y"], labels=c(expression('C'['i/i'])), srt=90, cex=2)
text(x=lblc["lab","x"], y=lblc["lab","y"],labels=paste("(",letters[1],")"),cex=2)

# B NORM MARGS CLAY COP
sampsB <- theoryfigpanel(dnmarg, pnmarg, qnmarg, ccop, x, y, samps=sampsA, oldp=pbmarg,xpd=NA)
text(x=lblc["E","x"], y=lblc["E","y"], labels=c(expression(phi^-1*'o'*'F'['E'['i']]("E"["i"]))), cex=2)
text(x=lblc["C","x"], y=lblc["C","y"], labels=c(expression(phi^-1*'o'*'F'['C'['i/i']]("C"["i/i"]))), srt=90, cex=2)
text(x=lblc["lab","x"], y=lblc["lab","y"],labels=paste("(",letters[2],")"), cex=2)

# C NORM MARGS NORM COP
sampsC <- theoryfigpanel(dnmarg, pnmarg, qnmarg, ncop, x, y, gnrt_var=TRUE, xpd=NA)
text(x=lblc["E","x"], y=lblc["E","y"], labels=c(expression('e'['i'])), cex=2)
text(x=lblc["C","x"], y=lblc["C","y"], labels=c(expression('c'['i/i'])), srt=90, cex=2)
text(x=lblc["lab","x"], y=lblc["lab","y"],labels=paste("(",letters[3],")"), cex=2)

# D BETA MARGS NORM COP
sampsD<-theoryfigpanel(dbmarg, pbmarg, qbmarg, ncop, x, y,samps=sampsC, oldp=pnmarg)
text(x=lblc["E","x"], y=lblc["E","y"], labels=c(expression('E'['i']^'||'=='F'['E'['i']]^-1*'o'*phi('e'['i']))), cex=2)
text(x=lblc["C","x"], y=lblc["C","y"], labels=c(expression('C'['i/i']^'||'=='F'['C'['i/i']]^-1*'o'*phi('c'['i/i']))), srt=90, cex=2)
text(x=lblc["lab","x"], y=lblc["lab","y"],labels=paste("(",letters[4],")"), cex=2)