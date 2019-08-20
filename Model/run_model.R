#### Run community occupancy model in JAGS
#### by Jesse F. Abrams

library(snowfall)
library(rjags)
library(rlecuyer)
library(ggplot2)

wd <- "~/Dropbox (ScreenForBio)/Projects/Species_comparison/Model reruns/Sabah"
setwd(wd)

mdata<-dget("occ.in_sabah.R")
species <- read.csv("~/Dropbox (ScreenForBio)/Projects/Species_comparison/species_list/Sabah_species_list.csv") 

##save only BP values 
params<-c("alpha",
          "mu.b1","sigma.b1",
          "mu.b2","sigma.b2",
          "mu.br", "sigma.br",
          "beta1","beta2",
          "b.eff")

data2<-mdata
data2$y<-data2$y[1:15,,]  
data2$M<-15              ### change the munber of recorded species and undetected species

zin<-apply(data2$y,1:2, sum)
zin[zin>1]<-1

inits2<-function(){list( mu.a=runif(3, -.5, .5), sig.a=runif(3, 0.5, 1.5), 
                         mu.bp=runif(1, 0, 2), 
                         tau.bp=runif(1, 0.5, 1),
                         mu.br=runif(1), tau.br=1, 
                         mu.b1=runif(1),
                         tau.b1=runif(1, 0.5, 1),
                         mu.b2=runif(1),
                         tau.b2=runif(1, 0.5, 1),
                         b.eff=runif(1), z=zin  )}


filename <- "com_model.txt"
modelFile2<-filename

wrapper2<-function(a){
  mod<-jags.model(modelFile2, data2, inits2, n.chain=1, n.adapt=500)
  out<-coda.samples(mod,params, n.iter=200000, thin=20 )
  return(out)
}

sfInit( parallel=TRUE, cpus=3 )
sfLibrary(rjags)
sfExportAll()
sfClusterSetupRNG()
time.2 <- system.time(
  (out.sp=sfLapply(1:3,wrapper2))
)

sfStop()

dput(out.sp, "Sabah_results_CAN_HEAT.R")
