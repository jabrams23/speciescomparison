#### analysis of site specific occupancies and community occupancy defaunation index
#### by Jesse F. Abrams


##################################################
# load required libraries
library(snowfall)
library(rjags)
library(rlecuyer)
library(matrixStats)
library(readxl)
library(stringr)
library(reshape)
library(ggpubr)
library(grid)
library(gridExtra)
library(lattice)
library(cowplot)
library(egg) 
library(ggplot2)
##################################################

##################################################
#set working directory
wd <- " "
setwd(wd)
##################################################


##################################################
# load occupancies calculated in file "psi_calc.R"
dfr <- dget("Sabah/psi_DFR_FINAL.R")
tfr <- dget("Sabah/psi_TFR_FINAL.R")
kfr <- dget("Sabah/psi_KFR_FINAL.R")
vbc <- dget("Vietnam/psi_VBC_FINAL.R")
vsc <- dget("Vietnam/psi_VSC_FINAL.R")
xs <- dget("Vietnam/psi_XS_FINAL.R")
##################################################

##################################################
# get statistics for each site
mean_dfr <- colMeans(dfr)
sd_dfr <- colSds(as.matrix(dfr))
quant_dfr <- colQuantiles(as.matrix(dfr), probs = c(2.5, 97.5)/100)
dfr_out <- data.frame(mean=mean_dfr,
                      sd=sd_dfr,
                      lower=quant_dfr[,1],
                      upper=quant_dfr[,2])
mean_tfr <- colMeans(tfr)
sd_tfr <- colSds(as.matrix(tfr))
quant_tfr <- colQuantiles(as.matrix(tfr), probs = c(2.5, 97.5)/100)
tfr_out <- data.frame(mean=mean_tfr,
                      sd=sd_tfr,
                      lower=quant_tfr[,1],
                      upper=quant_tfr[,2])
mean_kfr <- colMeans(kfr)
sd_kfr <- colSds(as.matrix(kfr))
quant_kfr <- colQuantiles(as.matrix(kfr), probs = c(2.5, 97.5)/100)
kfr_out <- data.frame(mean=mean_kfr,
                      sd=sd_kfr,
                      lower=quant_kfr[,1],
                      upper=quant_kfr[,2])


mean_bm <- colMeans(vbc)
sd_bm <- colSds(as.matrix(vbc))
quant_bm <- colQuantiles(as.matrix(vbc), probs = c(2.5, 97.5)/100)
bm_out <- data.frame(mean=mean_bm,
                      sd=sd_bm,
                      lower=quant_bm[,1],
                      upper=quant_bm[,2])
mean_sn <- colMeans(vsc)
sd_sn <- colSds(as.matrix(vsc))
quant_sn <- colQuantiles(as.matrix(vsc), probs = c(2.5, 97.5)/100)
sn_out <- data.frame(mean=mean_sn,
                      sd=sd_sn,
                      lower=quant_sn[,1],
                      upper=quant_sn[,2])
mean_xs <- colMeans(xs)
sd_xs <- colSds(as.matrix(xs))
quant_xs <- colQuantiles(as.matrix(xs), probs = c(2.5, 97.5)/100)
xs_out <- data.frame(mean=mean_xs,
                      sd=sd_xs,
                      lower=quant_xs[,1],
                      upper=quant_xs[,2])
##################################################

##################################################
# lazy but easy way to reorder factors
colnames(dfr) <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
colnames(tfr) <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
colnames(kfr) <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
colnames(vbc) <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
colnames(vsc) <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
colnames(xs) <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
##################################################



##################################################
# collapse dataframes
dfr_collapse <- melt(dfr)
tfr_collapse <- melt(tfr)
kfr_collapse <- melt(kfr)
vbc_collapse <- melt(vbc)
vsc_collapse <- melt(vsc)
xs_collapse <- melt(xs)
##################################################

##################################################

species_pairs <- read.csv("~/Dropbox (ScreenForBio)/Projects/Species_comparison/species_list/species_pairs_combined2.csv")
colnames(species_pairs)[1]<-"order"
species_pairs$species <- as.character(species_pairs$pairs)
##################################################

##################################################
# rename factors
dfr_collapse$Site <- "Deramakot FR"
tfr_collapse$Site <- "Tangkulap FR"
kfr_collapse$Site <- "Kuamut FR"
vbc_collapse$Site <- "Bach Ma NP"
vsc_collapse$Site <- "Saola NR"
xs_collapse$Site <- "Xe Sap NPA"
##################################################

##################################################
# create dataframe for plotting
dat4plot <- rbind(dfr_collapse,
                  tfr_collapse,
                  kfr_collapse,
                  vbc_collapse,
                  vsc_collapse,
                  xs_collapse)
colnames(dat4plot)[1] <- "order"
dat4plot$Site <- as.factor(dat4plot$Site)
##################################################

##################################################
sabah_list <- str_extract(species_pairs$species, "/.*")
vietnam_list<- str_extract(species_pairs$species, ".*/")
sabah_list <- gsub(" /","",sabah_list) 
vietnam_list <- gsub("/ ","",vietnam_list) 
species_axis <- data.frame(vietnam=vietnam_list,sabah=sabah_list)
##################################################


##################################################
cols <- c("#99ddff",
          "#4dc3ff",
          "#0088cc",
          "#ffdc8f",
          "#edaf13",
          "#a87602"
          )
##################################################

##################################################
dat4plot$Site <- factor(dat4plot$Site, levels = c("Xe Sap NPA","Saola NR","Bach Ma NP",
                                                  "Deramakot FR","Kuamut FR","Tangkulap FR"))

test <- as.character(species_pairs$species)
test <- gsub('/','\n',test)

dat4plot$species <- dat4plot$order
dat4plot$species <- as.character(dat4plot$species)

dat4plot$species[dat4plot$species == 1] <- test[1]
dat4plot$species[dat4plot$species == 2] <- test[2]
dat4plot$species[dat4plot$species == 3] <- test[3]
dat4plot$species[dat4plot$species == 4] <- test[4]
dat4plot$species[dat4plot$species == 5] <- test[5]
dat4plot$species[dat4plot$species == 6] <- test[6]
dat4plot$species[dat4plot$species == 7] <- test[7]
dat4plot$species[dat4plot$species == 8] <- test[8]
dat4plot$species[dat4plot$species == 9] <- test[9]
dat4plot$species[dat4plot$species == 10] <- test[10]
dat4plot$species[dat4plot$species == 11] <- test[11]
dat4plot$species[dat4plot$species == 12] <- test[12]
dat4plot$species[dat4plot$species == 13] <- test[13]
dat4plot$species[dat4plot$species == 14] <- test[14]
dat4plot$species[dat4plot$species == 15] <- test[15]

dat4plot <- as.data.frame(dat4plot)

three <- dat4plot[dat4plot$order == c(4,6,8,10,11,13,14,15),]
three_species <- test[c(4,6,8,10,11,13,14,15)]
two <- dat4plot[dat4plot$order == 9,]
two_species <- test[c(9)]
one <- dat4plot[dat4plot$order == c(3,7,12),]
one_species <- test[c(3,7,12)]
zero <- dat4plot[dat4plot$order == c(1,2,5),]
zero_species <- test[c(1,2,5)]
##################################################

##################################################
p_three <- ggplot(three, aes(species, value, fill = Site)) +
  geom_boxplot(outlier.shape = NA, position=position_dodge(1)) +
  scale_fill_manual(breaks=list_order,values = cols) +
  coord_flip() +
  ylab("") +
  geom_vline(xintercept=1.5, linetype="solid", color = "grey",size=1) +
  geom_vline(xintercept=2.5, linetype="solid", color = "grey",size=1) +
  geom_vline(xintercept=3.5, linetype="solid", color = "grey",size=1) +
  geom_vline(xintercept=4.5, linetype="solid", color = "grey",size=1) +
  geom_vline(xintercept=5.5, linetype="solid", color = "grey",size=1) +
  geom_vline(xintercept=6.5, linetype="solid", color = "grey",size=1) +
  geom_vline(xintercept=7.5, linetype="solid", color = "grey",size=1) +
  labs(y=expression(paste("Occupancy probability (", psi, ")",sep=""))) +
  guides(fill=guide_legend(ncol=2)) +
  theme_bw() +
  theme(axis.title.x=element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.position="bottom",
        legend.text=element_text(size=22),
        legend.title=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())


p_two <- ggplot(two, aes(species, value, fill = Site)) +
  geom_boxplot(outlier.shape = NA, position=position_dodge(1)) +
  scale_fill_manual(breaks=list_order,values = cols) +
  coord_flip() +
  ylab("") +
  labs(y=expression(paste("Occupancy probability (", psi, ")",sep=""))) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_text(size=22),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

p_one <- ggplot(one, aes(species, value, fill = Site)) +
  geom_boxplot(outlier.shape = NA, position=position_dodge(1)) +
  scale_fill_manual(breaks=list_order,values = cols) +
  coord_flip() +
  ylab("") +
  geom_vline(xintercept=1.5, linetype="solid", color = "grey",size=1) +
  geom_vline(xintercept=2.5, linetype="solid", color = "grey",size=1) +
  labs(y=expression(paste("Occupancy probability (", psi, ")",sep=""))) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_text(size=22),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())


p_zero <- ggplot(zero, aes(species, value, fill = Site)) +
  geom_boxplot(outlier.shape = NA, position=position_dodge(1)) +
  scale_fill_manual(breaks=list_order,values = cols) +
  coord_flip() +
  ylab("") +
  geom_vline(xintercept=1.5, linetype="solid", color = "grey",size=1) +
  geom_vline(xintercept=2.5, linetype="solid", color = "grey",size=1) +
  labs(y=expression(paste("Occupancy probability (", psi, ")",sep=""))) +
  guides(fill=guide_legend(ncol=2)) +
  theme_bw() +
  theme(axis.title.x=element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.position="bottom",
        legend.text=element_text(size=22),
        legend.title=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

right <- ggarrange(p_two, p_one, p_zero, ncol=1, nrow=3, heights=c(1,3,3))
right <- ggarrange(p_two, p_one, p_zero, ncol=1, nrow=3, heights=c(1,3,3))
top <- grid.arrange(p_three, right, ncol = 2)
top <- grid.arrange(p_three, right, ncol = 2)

#####################################################
## MONTE CARLO DRAWS FOR OCCUPANCY DEFAUNATION INDEX
#####################################################

dfr <- dget("Sabah/psi_DFR_FINAL.R")
tfr <- dget("Sabah/psi_TFR_FINAL.R")
kfr <- dget("Sabah/psi_KFR_FINAL.R")
vbc <- dget("Vietnam/psi_VBC_FINAL.R")
vsc <- dget("Vietnam/psi_VSC_FINAL.R")
xs <- dget("Vietnam/psi_XS_FINAL.R")

species_pairs <- read.csv("~/Dropbox (ScreenForBio)/Projects/Species_comparison/species_list/species_pairs_combined.csv")
colnames(species_pairs)[1]<-"order"

omegas <- read.csv("~/Dropbox (ScreenForBio)/Projects/Species_comparison/defaunation_data/comparison_omegas2.csv")
omegas<-as.data.frame(omegas)
rownames(omegas)<-omegas$`Common.name`

omegas_sort <-omegas[c(3,5,4,1,8,9,11,10,12,2,6,14,13,7,15),]

omega_listed <- omegas_sort$`Conservation.wk`
omega_size <- omegas_sort$`Size.wk`

dfr <- t(dfr)
tfr <- t(tfr)
kfr <- t(kfr)
vbc <- t(vbc)
vsc <- t(vsc)
xs <- t(xs)
nsim <- ncol(dfr)
nspecies <- nrow(dfr)

# or use a species importance where all species are treated the same
omega_equal <- rep(1, nspecies)

defaunation_equal <- array(0, c(nsim, 5))
defaunation_size <- array(0, c(nsim, 5))
defaunation_listed <- array(0, c(nsim, 5))

focus <- array(0, c(15, 22500, 5))
focus[,,1] <- tfr
focus[,,2] <- kfr
focus[,,3] <- vbc
focus[,,4] <- vsc
focus[,,5] <- xs

for (x in 1:nsim) {
  for (i in 1:5) {
    
    d_top <- omega_equal*(dfr[,x] - focus[,x,i])
    d_bottom <- omega_equal*(dfr[,x] + focus[,x,i])
    defaunation_equal[x,i] <- (sum(d_top))/(sum(d_bottom))
    
    d_top <- omega_size*(dfr[,x] - focus[,x,i])
    d_bottom <- omega_size*(dfr[,x] + focus[,x,i])
    defaunation_size[x,i] <- (sum(d_top))/(sum(d_bottom))
    
    d_top <- omega_listed*(dfr[,x] - focus[,x,i])
    d_bottom <- omega_listed*(dfr[,x] + focus[,x,i])
    defaunation_listed[x,i] <- (sum(d_top))/(sum(d_bottom))
    
  }
}

hold <- 0


defaunation_equal <- as.data.frame(defaunation_equal)
defaunation_size <- as.data.frame(defaunation_size)
defaunation_listed <- as.data.frame(defaunation_listed)

colnames(defaunation_equal) <- c("TFR",
                                 "KFR",
                                 "BMNP",
                                 "SNR",
                                 "XS")
colnames(defaunation_size) <- c("TFR",
                                "KFR",
                                "BMNP",
                                "SNR",
                                "XS")
colnames(defaunation_listed) <- c("TFR",
                                  "KFR",
                                  "BMNP",
                                  "SNR",
                                  "XS")

defaunation_equal2 <- melt(defaunation_equal)
defaunation_size2 <- melt(defaunation_size)
defaunation_listed2 <- melt(defaunation_listed)

defaunation_equal2$weight <- "equal"
defaunation_size2$weight <- "size"
defaunation_listed2$weight <- "conservation"

defaunation <- defaunation_equal2

dat2plot <- defaunation


mean_defaunation_equal <- colMeans(defaunation_equal)
sd_equal <- colSds(as.matrix(defaunation_equal))
quant_equal <- colQuantiles(as.matrix(defaunation_equal), probs = c(2.5, 97.5)/100)

mean_defaunation_size <- colMeans(defaunation_size)
sd_size <- colSds(as.matrix(defaunation_size))
quant_size <- colQuantiles(as.matrix(defaunation_size), probs = c(2.5, 97.5)/100)

mean_defaunation_listed <- colMeans(defaunation_listed)
sd_listed <- colSds(as.matrix(defaunation_listed))
quant_listed <- colQuantiles(as.matrix(defaunation_listed), probs = c(2.5, 97.5)/100)

quant_equal <- t(quant_equal)
quant_size <- t(quant_size)
quant_listed <- t(quant_listed)

stats_equal <- rbind(mean_defaunation_equal,sd_equal,quant_equal)
stats_size <- rbind(mean_defaunation_size,sd_size,quant_size)
stats_listed <- rbind(mean_defaunation_listed,sd_listed,quant_listed)

stats_equal <- as.data.frame(stats_equal)
stats_size <- as.data.frame(stats_size)
stats_listed <- as.data.frame(stats_listed)

rownames(stats_equal)<-c("mean","sd","lower","upper")
rownames(stats_size)<-c("mean","sd","lower","upper")
rownames(stats_listed)<-c("mean","sd","lower","upper")

stats_equal2 <- as.data.frame(t(stats_equal))
stats_size2 <- as.data.frame(t(stats_size))
stats_listed2 <- as.data.frame(t(stats_listed))

stats_equal2$Weight <- "equal"
stats_equal2$Site <- rownames(stats_equal2)
stats_size2$Weight <- "size"
stats_size2$Site <- rownames(stats_size2)
stats_listed2$Weight <- "conservation"
stats_listed2$Site <- rownames(stats_listed2)

dat4plot <- stats_equal2

dat4plot$Site <- factor(dat4plot$Site,levels=c("BMNP","SNR","XS","TFR","KFR"))
dat2plot$Site <- factor(dat2plot$Site,levels=c("BMNP","SNR","XS","TFR","KFR"))

colnames(dat2plot) <- c("Site","value","Weight")

dat4plot$cols <- as.character(dat4plot$Site)
dat2plot$cols <- as.character(dat2plot$Site)


cols <- c("#99ddff",
          "#4dc3ff",
          "#0088cc",
          "#edaf13",
          "#a87602"
)



dat4plot$cols <- gsub("TFR", "#a87602", dat4plot$cols)
dat4plot$cols <- gsub("KFR", "#edaf13", dat4plot$cols)
dat4plot$cols <- gsub("BMNP", "#0088cc", dat4plot$cols)
dat4plot$cols <- gsub("SNR", "#4dc3ff", dat4plot$cols)
dat4plot$cols <- gsub("XS", "#99ddff", dat4plot$cols)

dat2plot$cols <- gsub("TFR", "#a87602", dat2plot$cols)
dat2plot$cols <- gsub("KFR", "#edaf13", dat2plot$cols)
dat2plot$cols <- gsub("BMNP", "#0088cc", dat2plot$cols)
dat2plot$cols <- gsub("SNR", "#4dc3ff", dat2plot$cols)
dat2plot$cols <- gsub("XS", "#99ddff", dat2plot$cols)


rects <- data.frame(xstart = 0, xend = 0.5, ymin = 0, ymax = 120, col = letters[1])


dat4plot$Weight <- as.factor(dat4plot$Weight)
dat4plot$Weight <- factor(dat4plot$Weight,levels=c("equal","conservation","size"))

dat2plot$Weight <- as.factor(dat2plot$Weight)
dat2plot$Weight <- factor(dat2plot$Weight,levels=c("equal","conservation","size"))

cols <- rev(cols)

p_def <- ggplot(dat2plot, aes(x=value, fill=Site, color=Site)) +
  geom_histogram(data=dat2plot, aes(x=value, fill=Site, color=Site), position="dodge",binwidth=0.001,alpha=0.5,na.rm=TRUE) +
  geom_vline(data=dat4plot, aes(xintercept=mean, color=Site), color = dat4plot$cols) + 
  geom_vline(data=dat4plot, aes(xintercept=(mean-sd), color=Site),linetype="longdash", color = dat4plot$cols) + 
  geom_vline(data=dat4plot, aes(xintercept=(mean+sd), color=Site),linetype="longdash", color = dat4plot$cols) + 
  xlim(-0.25,0.8) +
  ylab("") +
  xlab("Defaunation Index") +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) +
  annotate("rect",xmin=-0.02,xmax=0.02, ymin=-Inf, ymax=Inf, fill="darkgrey", alpha=0.7) +
  annotate("text", x=0, y=100, size=7, label="reference site, DFR", color="white", angle=90) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="bottom",
        axis.text.x=element_text(size=20),
        axis.title.x=element_text(size=22),
        legend.text=element_text(size=20),
        legend.title=element_text(size=22),
        plot.tag=element_text(size=22, face="bold"))


fig <- grid.arrange(top, p_def, nrow = 2,heights=c(2,1))
fig <- grid.arrange(top, p_def, nrow = 2,heights=c(2,1))
fig <- grid.arrange(top, p_def, nrow = 2,heights=c(2,1))

png('Figure2_ACCEPTED.png', width = 15, height = 15, units = 'in', res = 300)
grid.draw(fig)
dev.off()

setEPS()
postscript("Figure2_ACCEPTED.eps")
grid.draw(fig)
dev.off()
