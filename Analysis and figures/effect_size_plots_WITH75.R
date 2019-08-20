#### analysis and plotting of effect sizes of covariates
#### by Jesse F. Abrams

library(rjags)
library(snowfall)
library(rlecuyer)
library(ggplot2)
library(reshape)
library(ggpubr)
library(grid)
library(gridExtra)
library(lattice)
library(cowplot)
library(egg) 

wd <- ""
setwd(wd)

#--------------------------------------
# Vietnam
#--------------------------------------


# load model output and summarize
##################################################
out.sp <-dget("Vietnam/Vietnam_results_CAN_HEAT.R")
out.f<-mcmc.list(mcmc(out.sp[[1]][[1]]), mcmc(out.sp[[2]][[1]]), mcmc(out.sp[[3]][[1]]))
sm<-summary(window(out.f, start=2500))
##################################################


species <- read.csv("~/Dropbox (ScreenForBio)/Projects/Species_comparison/species_list/VN_species_list.csv")

subset<-1:15
nsp2<-15
ssp<-as.character(species$Species)
newsp<-ssp
best<-sm

##CANOPY
betf<-grep( "beta1", dimnames(best$statistics)[[1]]  )
bet<-betf[subset]
resmat<-cbind(best[[1]][,1:2], best[[2]][,c(1,5)])
nums<-resmat[bet, c(1,3,4)]
lim<-c(min(nums[,2]), max(nums[,3]))
mu<-grep( "mu.b1", dimnames(best$statistics)[[1]]  )
resmat[mu, c(1,3,4)]
muv<-resmat[mu,1]
muv_CAN<- muv
viet_CAN <- data.frame(nums, newsp)

##HEAT
betf<-grep( "beta2", dimnames(best$statistics)[[1]]  )
bet<-betf[subset]
resmat<-cbind(best[[1]][,1:2], best[[2]][,c(1,5)])
nums<-resmat[bet, c(1,3,4)]
lim<-c(min(nums[,2]), max(nums[,3]))
mu<-grep( "mu.b2", dimnames(best$statistics)[[1]]  )
muv<-resmat[mu,1]
muv_HEAT <- muv
viet_HEAT <- data.frame(nums, newsp)

#--------------------------------------
# SABAH
#--------------------------------------

# load model output and summarize
##################################################
out.sp <-dget("Sabah/Sabah_results_CAN_HEAT.R")
out.f<-mcmc.list(mcmc(out.sp[[1]][[1]]), mcmc(out.sp[[2]][[1]]), mcmc(out.sp[[3]][[1]]))
sm<-summary(window(out.f, start=2500))
##################################################

species <- read.csv("~/Dropbox (ScreenForBio)/Projects/Species_comparison/species_list/Sabah_species_list.csv") 
subset<-1:15
nsp2<-15
ssp<-as.character(species$Species)
newsp<-ssp
best<-sm

##CANOPY
betf<-grep( "beta1", dimnames(best$statistics)[[1]]  )
bet<-betf[subset]
resmat<-cbind(best[[1]][,1:2], best[[2]][,c(1,5)])
nums<-resmat[bet, c(1,3,4)]
lim<-c(min(nums[,2]), max(nums[,3]))
mu<-grep( "mu.b1", dimnames(best$statistics)[[1]]  )
resmat[mu, c(1,3,4)]
muv<-resmat[mu,1]
mus_CAN <- muv
sabah_CAN <- data.frame(nums, newsp)

##HEAT
betf<-grep( "beta2", dimnames(best$statistics)[[1]]  )
bet<-betf[subset]
resmat<-cbind(best[[1]][,1:2], best[[2]][,c(1,5)])
nums<-resmat[bet, c(1,3,4)]
lim<-c(min(nums[,2]), max(nums[,3]))
mu<-grep( "mu.b2", dimnames(best$statistics)[[1]]  )
muv<-resmat[mu,1]
mus_HEAT <- muv
sabah_HEAT <- data.frame(nums, newsp)



colnames(viet_CAN) <- colnames(sabah_CAN)
CAN <-rbind(viet_CAN,sabah_CAN)
order <- seq(1,30)
Site <- c("Vietnam","Vietnam","Vietnam","Vietnam","Vietnam","Vietnam","Vietnam",
          "Vietnam","Vietnam","Vietnam","Vietnam","Vietnam","Vietnam","Vietnam","Vietnam",
          "Sabah","Sabah","Sabah","Sabah","Sabah","Sabah","Sabah",
          "Sabah","Sabah","Sabah","Sabah","Sabah","Sabah","Sabah","Sabah")
CAN <- cbind(CAN, Site, order)
CAN$Site <- factor(CAN$Site, levels = c("Vietnam", "Sabah"))
levels(CAN$Site)[levels(CAN$Site)=="Sabah"] <- "Degraded"
levels(CAN$Site)[levels(CAN$Site)=="Vietnam"] <- "Hunted"

colnames(viet_HEAT) <- colnames(sabah_HEAT)
HEAT <-rbind(viet_HEAT,sabah_HEAT)
HEAT <- cbind(HEAT, Site, order)
HEAT$Site <- factor(HEAT$Site, levels = c("Vietnam", "Sabah"))
levels(HEAT$Site)[levels(HEAT$Site)=="Sabah"] <- "Degraded"
levels(HEAT$Site)[levels(HEAT$Site)=="Vietnam"] <- "Hunted"


viet_CAN <- viet_CAN[ order(row.names(viet_CAN)), ]
viet_HEAT <- viet_HEAT[ order(row.names(viet_HEAT)), ]
viet_all <- data.frame(species=viet_CAN$newsp,
                       "Village heat map"=viet_HEAT$Mean,
                       "Canopy closure"=viet_CAN$Mean)
sabah_CAN <- sabah_CAN[ order(row.names(sabah_CAN)), ]
sabah_HEAT <- sabah_HEAT[ order(row.names(sabah_HEAT)), ]
sabah_all <- data.frame(species=sabah_CAN$newsp,
                       "Village heat map"=sabah_HEAT$Mean,
                       "Canopy closure"=sabah_CAN$Mean)
all_combined <- rbind(viet_all, sabah_all)
all_combined <- cbind(Site, order,all_combined)
all_combined$Site <- factor(all_combined$Site, levels = c("Vietnam", "Sabah"))
levels(all_combined$Site)[levels(all_combined$Site)=="Sabah"] <- "Degraded"
levels(all_combined$Site)[levels(all_combined$Site)=="Vietnam"] <- "Hunted"

all_combined2 <- all_combined
all_combined2$order <- NULL
all_combined2$species <- NULL
colnames(all_combined2) <- c("Site","Village heat map","Canopy closure")

heat_sorted <- HEAT
heat_sorted$nam.ord <- as.character(heat_sorted$newsp)
heat_sorted$which <- "Village heat map"
canopy_sorted <- CAN
canopy_sorted$nam.ord <- as.character(canopy_sorted$newsp)
canopy_sorted$which <- "Canopy closure"

allcombined3 <- rbind(canopy_sorted,heat_sorted)
allcombined3$ord <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
         1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)

write.csv(allcombined3,"effect_size_table.csv")


cols1 <- c("#005580",
          "#b37700",
          "#4dc3ff",
          "#99ddff",
          "#0088cc",
          "#ffaa00",
          "#ffc34d",
          "#ffd480")


allcombined3$species <- allcombined3$newsp
levels(allcombined3$species) <- allcombined3$nam.ord[1:30]

p <- ggplot(data=allcombined3, aes(x = Mean,y = order,colour=Site)) +
  geom_point(size=2) +
  geom_errorbarh(aes(xmin = X2.5.,xmax = X97.5.),size=0.6,height=0,color="gray48") +
  geom_errorbarh(aes(xmin = X25.,xmax = X75.),size=1.2,height=0,color="brown3") +
  scale_color_manual(breaks=c("Degraded","Hunted"),values = c(cols1[1],cols1[2])) +
  scale_color_manual(breaks=c("Degraded","Hunted"),values = c(cols1[1],cols1[2])) +
  geom_vline(xintercept=0, linetype="dashed", color = "black",size=0.5) +
  scale_y_continuous(breaks=order,labels=(as.character(allcombined3$nam.ord[1:30]))) +
  geom_errorbarh(aes(xmin = X2.5.,xmax = X97.5.),size=0.6,height=0,color="black") +
  geom_errorbarh(aes(xmin = X25.,xmax = X75.),size=1.2,height=0,color="brown3") +
  geom_point(data=allcombined3, aes(x = Mean,y = order,colour=Site),size=2) +
  ylab("Species") +
  coord_flip() +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        axis.title.x=element_text(size=16, face="bold"),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=14,angle = 90, hjust = 1,vjust=0.3),
        axis.text.y=element_text(size=14),
        legend.position='none',
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

p2_eff <- p + labs(tag = "A")+ facet_grid(which ~ Site , switch="y", scales = "free", space = "free")  +
  theme(strip.text.y = element_text(size=16, face="bold"),
        strip.text.x = element_text(size=16, face="bold"),
        strip.placement = "outside",
        plot.tag=element_text(size=18, face="bold"))

p_scat <- ggplot(data=all_combined2, aes(y = `Village heat map`,x = `Canopy closure`,color=Site,group=Site)) +
  geom_point(size=3) +
  stat_ellipse(aes(color = Site), type = "t") +
  scale_color_manual(breaks=c("Degraded","Hunted"),values = c(cols1[1],cols1[2])) +
  xlab("Canopy closure") +
  ylab("Village heat map") +
  labs(tag = "B") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white",
                                         colour = "white",
                                         size = 0.5, linetype = "solid"),
        axis.title.x=element_text(size=16, face="bold"),
        axis.title.y=element_text(size=16, face="bold"),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        legend.position='bottom',
        legend.text=element_text(size=14),
        legend.title=element_text(size=16),
        plot.tag=element_text(size=18, face="bold"))

right <- grid.arrange(p_scat, nrow = 2,heights=c(3,1))
right2 <- right
p2_eff2 <- p2_eff
fig <- grid.arrange(p2_eff2, rectGrob(gp=gpar(col=NA)), right2, ncol = 3,widths=c(2.7,0.3,2))

png('figure3.png', width = 16, height = 8, units = 'in', res = 300)
grid.draw(fig)
dev.off()
