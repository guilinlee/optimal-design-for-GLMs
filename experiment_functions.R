#experiment=cbind(1,E[,1],E[,2],E[,3],E[,1]^2,E[,1]*E[,2],E[,2]*E[,2])

ufun=function(design,beta)
{ 
  x1=design[1];x2=design[2];x3=design[3]
  b0=beta[1];b1=beta[2];b2=beta[3];b3=beta[4];b4=beta[5];b5=beta[6];b6=beta[7]
  eta=b0+b1*x1+b2*x2+b3*x3+b4*x1^2+b5*x1*x2+b6*x2*x2
  u=exp(eta) 
}
optim_ufun=function(beta)
{
  est=optim(c(0.1,0.1,0.1),ufun,lower=c(-1,-1,-1),upper=c(1,1,1),
            method="L-BFGS-B",beta=beta)
  optimdesign=est$par
  return(optimdesign)
}
aggdpois=function(lambday)
{prob=dpois(x=lambday[-1],lambda=lambday[1])}
#########################################################1 prepare the prior
#1)sample beta from gaussian
discrete_beta=function(mean,sd,betaN)
{
  u_prior_mean=mean; u_prior_sd=sd
  #1+beta*width+beta2*width^2+..
  #get a sample of beta
  betasample=matrix(0,betaN,dimcov)
  for (i in 1:dimcov)
  {betasample[,i]=rnorm(betaN,u_prior_mean[i],u_prior_sd[i])}
  return(betasample)
}
###### 2) get the X* table
Xstarfun=function(betasample)
{ 
  Xstar=c(1:betaN)
  optimdesign_sim=t(apply(betasample,1,optim_ufun))
  ###round up xstar up to two digits
  optimdesign_round=round(optimdesign_sim,digits = 2)
  #####re-incoding x*=(width,length)  
  ##adjust the range from 0-200 * 0-200
  re_optimdesign=(optimdesign_round+1)*100
  re_optimdesign=formatC(re_optimdesign,width=3,flag = "0") 
  Xstar=paste0(re_optimdesign[,1],re_optimdesign[,2], re_optimdesign[,3])
  Xstar=as.numeric(Xstar)
  ##look up table
  Xstartable=data.frame(cbind(betasample,Xstar,optimdesign_round))
  ##this result into a table of 4k list
  splitXstartable=split(Xstartable,Xstartable$Xstar)
  return(splitXstartable) 
}
missing_xstari=function(eleXstartable,texperiment,unitnum)
{
  ni=dim(eleXstartable)[1]#number of theta
  ##### calculte theta_j j=1..ni for xstar=xstar_i
  beta_xstari=eleXstartable[,(1:dimcov)]; beta_xstari=unname(data.matrix(beta_xstari))
  ##calculate poisson means j=1,..ni for experiment
  lambda=exp(beta_xstari%*%texperiment)  #row is experiment=(e1,...e_E)
  simy=apply(lambda,c(1,2),rpois,n=unitnum)
  logprob_ratio=c(1:ni)
  for (j in 1:ni)
  {
    batchratio=1
    for (e in 1: numE)
    { 
      ##for fixed experiment e, the mixlambda is lambda[,e], the corresponding sample is simy[,,e], where row is p(y|beta_i), column is unitn samples 
      lambdaymatrix=t(rbind(lambda[j,e],simy[,j,e]))
      trueprob=apply(lambdaymatrix,1,aggdpois)## p(yij）i for sample index, j for lambda index
      mixprob=sapply(lambda[,e],dpois,x=as.vector(simy[,j,e]))
      ##mixprob[[1]] is assum all from beta1
      ##each sample need mixprob/true prob
      mixtrueratio=mixprob/replicate(ni,trueprob)
      ##the ratio need to be a product of P(y1..yE|xstar)/P(y1...yE|theta) so it's a product of E dimension ratios
      batchratio=batchratio*mixtrueratio
    }
    logprob_ratio_j=log(rowMeans(batchratio)) # 
    logprob_ratio_j[logprob_ratio_j=-Inf]=0
    logprob_ratio[j]=mean(logprob_ratio_j)
  }
  diffHy_i=mean(logprob_ratio)
  return(diffHy_i)
}
# missing_xstar=function(splitXuse,weight,texperiment,unitnum)
# {
#   lmissing=sapply(splitXuse,missing_xstari,texperiment=texperiment,unitnum=unitnum)
#   H_xstari=unname(lmissing)
#   H_xstar=sum(weight*H_xstari)
#   return(H_xstar)
# }
missing_xstar=function(splitXuse,weight,texperiment,unitnum)
{ 
  n=length(splitXuse)
  lmissing=c(1:n)
  for(i in 1:n)
  { eleXstartable=splitXuse[[i]]
  lmissing[i]=missing_xstari(eleXstartable,texperiment,unitnum)
  }
  H_xstari=unname(lmissing)
  H_xstar=sum(weight*H_xstari)
  return(H_xstar)
}
Doptim_beta=function(elelambda,texperiment)
{
  elelambda=as.vector(elelambda)
  R=diag(1/u_prior_sd)
  ZQZ=texperiment%*%diag(elelambda)%*%t(texperiment)+R
  return(det(ZQZ))
}
Doptim=function(betasample,texperiment)
{
  lambda=exp(betasample%*%texperiment)  #row is experiment=(e1,...e_E)
  detI=c(1:betaN)
  for (i in 1:betaN)
  {
    detI[i]=Doptim_beta(lambda[i,],texperiment)
  }
  Dvalue=1/2*mean(log(detI))-dimcov/2*log(2*pi)-dimcov/2
  return(Dvalue)
  #-dimcov/2*log(2*pi)-dimcov/2
}
###I(y,xstar)=I(theta,y)+H(y|theta)-H(y|xstar) maximize I(y,xtar)
Uystar=function(designstring)
{
  E=matrix(designstring,nrow=numE,ncol=dimX,byrow = T)
  experiment=cbind(1,E[,1],E[,2],E[,3],E[,1]^2,E[,1]*E[,2],E[,2]*E[,2])
  texperiment=t(experiment)
  XTX=det(texperiment%*%experiment)
  if (XTX<=0)
  {return(-Inf)}
  else
  {
    Uystar=Doptim(betasample,texperiment)+missing_xstar(splitXuse,weight,texperiment,unitnum)
    return(Uystar)}
}
Doptimstring=function(designstring)
{
  E=matrix(designstring,nrow=numE,ncol=dimX,byrow = T)
  experiment=cbind(1,E[,1],E[,2],E[,3],E[,1]^2,E[,1]*E[,2],E[,2]*E[,2])
  texperiment=t(experiment)
  XTX=det(texperiment%*%experiment)
  if (XTX<=0)
  {return(-Inf)}
  else
  {
    lambda=exp(betasample%*%texperiment)  #row is experiment=(e1,...e_E)
    detI=apply(lambda,1,Doptim_beta,texperiment)
    Dvalue=1/2*mean(log(detI))-dimcov/2*log(2*pi)-dimcov/2
    return(Dvalue)}
}

