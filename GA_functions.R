load("elismpop.Rdata")
gamy_Population=function (object) 
{
  nvars=numE*dimX
  population <- matrix(as.double(NA), nrow = object@popSize-10, 
                       ncol = nvars)
  unirange=seq(-0.9,0.9,0.1)
  for (j in 1:nvars) {
    population[, j] <- sample(unirange,(object@popSize-10),replace=T)
  }
  population=rbind(population,elismpop)
  return(population)
}