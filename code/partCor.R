#Compute the partial Spearman correlation in the left halves and right halves 
#of the distributions, and the difference. Based on the Ghosh et al work.

partCor <- function(x,y, b, ties){
  Tot <- length(x)
  u <- rank(x, ties=ties)/(Tot+1); v <- rank(y, ties=ties)/(Tot+1)
  left <- (u+v) > 0 & (u+v) < b*2
  right <- (u+v) > b*2 & (u+v) < 2 
  corl <- sum((u[left]-mean(u))*(v[left]-mean(v)))/((Tot-1)*sqrt(var(u)*var(v)))
  corr <- sum((u[right]-mean(u))*(v[right]-mean(v)))/((Tot-1)*sqrt(var(u)*var(v)))	
  return(list(left=corl, right=corr, diff=corl-corr, Sh=corl+corr, S=cor(u,v, method="spearman"), uv=cbind(u,v), bounds=cbind(left, right)))
}