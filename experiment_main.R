##note: 3e+06 sample use 200s
#library(doParallel)
args=commandArgs()
install.packages("GA",repos="http://cran.r-project.org")
install.packages("doParallel",repos="http://cran.r-project.org")
require(doParallel)
require(GA)
#
#
source('experiment_functions.R')
source('GA_functions.R')
load("betasample.Rdata")
numE=15
dimX=3
betaN=2000
unitnum=50
## intercetp+x1+x2+x3+x1^2+x1*x2+x2^2
u_prior_mean=c(1.895,0.631,0.73,-0.367,0.266, -0.546, -0.04)
u_prior_sd=c(0.045,0.038,0.040,0.14,0.037,0.039,0.023)
dimcov=length(u_prior_mean)
betasample=discrete_beta(u_prior_mean,u_prior_sd,betaN=betaN)
splitXstartable=Xstarfun(betasample)
splitXuse=splitXstartable[which(lapply(splitXstartable,function(x) { nrow(x[1])})>1)]
weight=do.call(rbind,lapply(splitXuse,nrow))/betaN ##weight is used in missing_star
minv=rep(-1,dimX*numE);maxv=rep(1,dimX*numE)
##
popSize=20
GADoptimal <- ga(type = "real-valued", fitness = Doptimstring,
                  #population =gamy_Population,
                  #mutation=gareal_rsMutation,
                  elitism = base::max(1, round(popSize*0.1)),
                  pmutation = 0.2,
                    min = minv, max = maxv,
                    popSize = popSize, maxiter = 400,parallel ="snow",monitor=F,seed=12345)
########################
popSize=20
GA <- ga(type = "real-valued", fitness = Uystar,
         #mutation=gareal_rsMutation,
         elitism = base::max(1, round(popSize*0.1)),
         #population =gamy_Population,
         pmutation = 0.2,
         min = minv, max = maxv,
         popSize = popSize, maxiter = 400,parallel ="snow",monitor=F,seed=12345)
resultlist=list(GADoptimal,GA)
save(resultlist,file=args[5])

###
GA@fitnessValue
designstring=GA@solution
E=matrix(designstring,nrow=numE,ncol=dimX,byrow = T)
Ddesignstring=GADoptimal@solution
ED=matrix(Ddesignstring,nrow=numE,ncol=dimX,byrow = T)
plot(E[1,],E[3,])
points(ED[1,],ED[3,],col="red")
points(ED,col="red")
value=function(x1,x2,x3)
  {
  value=1.895+0.631*x1+0.73*x2-0.367*x3+0.266*x1^2-0.546*x1*x2-0.04*x2^2
}
  x1=seq(-1,1,0.1)
x2=seq(-1,1,0.1)
x3=0
grid=expand.grid(x1,x2)
surface=value(grid$Var1,grid$Var2,x3)
install.packages("plotly")
library(plotly)
packageVersion('plotly')
plot_ly(x=grid$Var1,y=grid$Var2,z=surface);add_surface()
plot(x=grid$Var1,surface)
