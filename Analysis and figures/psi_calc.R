##### extraction of mean psi for each degraded site, change input files for hunted sites
##### by Jesse F. Abrams

library(coda)
library(secr)
library(raster)
library(snowfall)
library(rjags)
library(rlecuyer)

wd <- ""
setwd(wd)

out.sp <-dget("Sabah_results_CAN_HEAT.R")
out.f<-mcmc.list(mcmc(out.sp[[1]][[1]]), mcmc(out.sp[[2]][[1]]), mcmc(out.sp[[3]][[1]]))

parms<-dimnames(out.f[[1]])[[2]]

inm<-seq(2501,10000)
mod<-rbind(out.f[[1]][inm,], out.f[[2]][inm,], out.f[[3]][inm,])
nsim<-dim(mod)[1]

spin<-1:15 ###change for listed species only

###extract parameters beforehand
a1<-mod[ ,pmatch(paste('alpha[1,',spin, ']', sep=''), parms)]
a2<- mod[ ,  pmatch(paste('alpha[2,', spin, ']', sep=''), parms)]
a3<-  mod[ ,  pmatch(paste('alpha[3,', spin, ']', sep=''), parms)]
a.mat<-array(NA, c(nsim, length(spin), 3))
a.mat[,,1]<-a1
a.mat[,,2]<-a2
a.mat[,,3]<-a3

b1<- mod[,  pmatch(paste('beta1[', spin, ']', sep=''), parms)]
b2<- mod[,  pmatch(paste('beta2[', spin, ']', sep=''), parms)]

####re-analysis of comm. occ model, to get at model fit

mdata<-dget("occ.in_sabah.R")
covs<-mdata

lpsi_dfr <- a.mat[,,1] + b1 * mean(covs$CAN) + b2 * mean(covs$HEAT)
lpsi_kfr <- a.mat[,,2] + b1 * mean(covs$CAN) + b2 * mean(covs$HEAT)
lpsi_tfr <- a.mat[,,3] + b1 * mean(covs$CAN) + b2 * mean(covs$HEAT)
psi_dfr <- exp(lpsi_dfr)/(1+exp(lpsi_dfr))
psi_kfr <- exp(lpsi_kfr)/(1+exp(lpsi_kfr))
psi_tfr <- exp(lpsi_tfr)/(1+exp(lpsi_tfr))

psi_dfr2 <- data.frame(psi_dfr[,1],
                       psi_dfr[,2],  
                       psi_dfr[,3],
                       psi_dfr[,4],
                       psi_dfr[,5],
                       psi_dfr[,6],
                       psi_dfr[,7],
                       psi_dfr[,8],
                       psi_dfr[,9],
                       psi_dfr[,10],
                       psi_dfr[,11],
                       psi_dfr[,12],
                       psi_dfr[,13],
                       psi_dfr[,14],
                       psi_dfr[,15]
)

psi_kfr2 <- data.frame(psi_kfr[,1],
                       psi_kfr[,2],  
                       psi_kfr[,3],
                       psi_kfr[,4],
                       psi_kfr[,5],
                       psi_kfr[,6],
                       psi_kfr[,7],
                       psi_kfr[,8],
                       psi_kfr[,9],
                       psi_kfr[,10],
                       psi_kfr[,11],
                       psi_kfr[,12],
                       psi_kfr[,13],
                       psi_kfr[,14],
                       psi_kfr[,15]
)

psi_tfr2 <- data.frame(psi_tfr[,1],
                       psi_tfr[,2],  
                       psi_tfr[,3],
                       psi_tfr[,4],
                       psi_tfr[,5],
                       psi_tfr[,6],
                       psi_tfr[,7],
                       psi_tfr[,8],
                       psi_tfr[,9],
                       psi_tfr[,10],
                       psi_tfr[,11],
                       psi_tfr[,12],
                       psi_tfr[,13],
                       psi_tfr[,14],
                       psi_tfr[,15]
)

species <- read.csv("~/Dropbox (ScreenForBio)/Projects/Species_comparison/species_list/Sabah_species_list.csv") 
colnames(psi_dfr2) <- as.character(species$Species)
colnames(psi_kfr2) <- as.character(species$Species)
colnames(psi_tfr2) <- as.character(species$Species)

dput(psi_dfr2, paste("psi_DFR_FINAL.R"))
dput(psi_kfr2, paste("psi_KFR_FINAL.R"))
dput(psi_tfr2, paste("psi_TFR_FINAL.R"))

