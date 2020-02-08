# Restricted Quadratic Spline 
#
# Steven Xu, 2020 (roughly based on SAS code by CJ Howe et al., 2011 Epidemiology)
#
# Arguments:
#
#    data         The data containing the variable(s) for which restricted
#                 quadratic spline is expanded.
#    x.name       The variables for which restircted quadratic spline is 
#                 expanded.
#    equal        If equidistant knots should be placed, otherwise percen-
#                 tile knots are placed.
#    probs        Percentile knots. If not specified default locations will
#                 be used.
#    slice        Wether restricted quadratic spline should be expanded 
#                 separatedly for each factor of a variable
#
# The added variables are named var1, var2, ...
#
# The return value is the original dataset appended with k-1 restricted 
# quadratic spline variables.

rqspline = function(data,x.name,k=4,equal=FALSE,probs=NULL,slice=NULL){
  if(k<3 || k>7) stop("k has to be between 3 and 7")
  if(!all(x.name%in%colnames(data))) stop("One of the variables specified is not in the dataset")
  if(!is.null(slice)){
    if(length(which(colnames(data)==slice))==0) stop(paste("There is no variable named",slice,"in the dataset"))
  }
  if(is.null(probs)){
    if(!equal){
      probs.list=list(
        q1 = c(0.05,0.5,0.95),
        q2 = c(0.05,0.35,0.65,0.95),
        q3 = c(0.05,0.28,0.5,0.73,0.95),
        q4 = c(0.05,0.23,0.41,0.59,0.77,0.95),
        q5 = c(0.03,0.18,0.34,0.5,0.66,0.82,0.98)
      )
      probs = probs.list[[k-2]]
    }else{
      probs = seq(0,1,length.out = k+2)[-c(1,k+2)]
    }
  }
  x.pos = which(colnames(data)%in%x.name)
  x.ori = as.matrix(data[,x.pos])
  p = ncol(x.ori)
  x.out = rep(list(matrix(0,nrow=nrow(x.ori),ncol=k-1)),p)
  if(is.null(slice)){
    x.qt = matrix(0,ncol=p,nrow=length(probs))
    for(m in 1:p){
      x.qt[,m] = quantile(as.numeric(x.ori[,m]),probs=probs,na.rm = T)
      for(i in 1:nrow(x.ori)){
        for(j in 1:(k-1))
          x.out[[m]][i,j] = (max(0,x.ori[i,m]-x.qt[j,m])^2-max(0,x.ori[i,m]-x.qt[k,m])^2)/(x.qt[k,m]-x.qt[1,m])
      }
    }
  }else{
    s.pos = which(colnames(data)==slice)
    s.vec = data[,s.pos]
    level = unique(s.vec)
    n.level = length(level)
    for(i in 1:n.level){
      cur.pos = which(s.vec==level[i])
      x.cur.qt = matrix(0,ncol=p,nrow=length(probs))
      cur.length = length(cur.pos)
      for(m in 1:p){
        x.cur.qt[,m] = quantile(as.numeric(x.ori[cur.pos,m]),probs=probs,na.rm = T) 
        for(j in 1:cur.length){
          for(l in 1:(k-1)){
            x.out[[m]][cur.pos[j],l] = (max(0,x.ori[cur.pos[j],m]-x.cur.qt[l,m])^2-max(0,x.ori[cur.pos[j],m]-x.cur.qt[k,m])^2)/(x.cur.qt[k,m]-x.cur.qt[1,m])
          }
        }
      }
    }
  }
  for(m in 1:p){
    colnames(x.out[[m]])<-paste(x.name[m],1:(k-1),sep="")
  }
  x.out.mat = do.call("cbind",x.out)
  x.final = data.frame(data,x.out.mat)
  return(x.final)
}
