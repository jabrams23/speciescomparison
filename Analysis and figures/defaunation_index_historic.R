#####################################################
# Code for defaunation index
# by Jesse F. Abrams
#####################################################

library(readxl)
library(ggplot2)
library(reshape2)
library(grid)

setwd("")

# read in
# VIETNAM
current_vietnam <- read_excel("~/Dropbox (ScreenForBio)/Projects/Species_comparison/defaunation_data/vietnam_laos_pres_updated.xlsx")
historic_vietnam <- read_excel("~/Dropbox (ScreenForBio)/Projects/Species_comparison/defaunation_data/vietnam_laos_historical_updated2.xls")

# SABAH
current_sabah <- read_excel("~/Dropbox (ScreenForBio)/Projects/Species_comparison/defaunation_data/sabah_pres_funct_updated.xlsx")
historic_sabah <- read_excel("~/Dropbox (ScreenForBio)/Projects/Species_comparison/defaunation_data/sabah_pres_funct_historical_updated2.xls")


current_bmnp_complete <- read_excel("~/Dropbox (ScreenForBio)/Projects/Species_comparison/defaunation_data/BMNP.xlsx")
current_snr_complete <- read_excel("~/Dropbox (ScreenForBio)/Projects/Species_comparison/defaunation_data/SNRs.xlsx")
current_laos_complete <- read_excel("~/Dropbox (ScreenForBio)/Projects/Species_comparison/defaunation_data/Laos.xlsx")
current_dfr_complete <- read_excel("~/Dropbox (ScreenForBio)/Projects/Species_comparison/defaunation_data/Deramakot.xlsx")
current_tfr_complete <- read_excel("~/Dropbox (ScreenForBio)/Projects/Species_comparison/defaunation_data/Tangkulap.xlsx")
current_kfr_complete <- read_excel("~/Dropbox (ScreenForBio)/Projects/Species_comparison/defaunation_data/Kuamut.xlsx")


current_occupancy_vietnam <- current_vietnam$`Present / absent`[1:45]
historic_occupancy_vietnam <- historic_vietnam$`Present / absent`[1:45]
current_bmnp <- current_bmnp_complete$`Present / absent`[1:45]
current_snr <- current_snr_complete$`Present / absent`[1:45]
current_laos <- current_laos_complete$`Present / absent`[1:45]

current_occupancy_sabah <- current_sabah$`Present / absent`[1:36]
historic_occupancy_sabah <- historic_sabah$`Present / absent`[1:36]
current_dfr <- current_dfr_complete$`Present / absent`[1:36]
current_tfr <- current_tfr_complete$`Present / absent`[1:36]
current_kfr <- current_kfr_complete$`Present / absent`[1:36]

# assign species importance

omega_vietnam_con <- historic_vietnam$`Conservation wk`[1:45]
omega_vietnam_size <- historic_vietnam$`Average body size`[1:45]
omega_vietnam_equal <- rep(1, length(historic_occupancy_vietnam))

omega_sabah_con <- historic_sabah$`Conservation wk`[1:36]
omega_sabah_size <- historic_sabah$`Average body size`[1:36]
omega_sabah_equal <- rep(1, length(historic_occupancy_sabah))



vietnam_pres <- data.frame(historic=historic_occupancy_vietnam,
                           vietnam_all=current_occupancy_vietnam,
                           bmnp=current_bmnp,
                           snr=current_snr,
                           laos=current_laos)


def_vietnam <- matrix(0L, nrow = 4, ncol = 3)
def_vietnam <- as.data.frame(def_vietnam)
colnames(def_vietnam) <- c("equal","conservation","size")
rownames(def_vietnam) <- colnames(vietnam_pres)[2:5]

for (i in 1:4) {
  d_top_con_vietnam <- omega_vietnam_con*(vietnam_pres$historic - vietnam_pres[,(i+1)])
  d_bottom_con_vietnam <- omega_vietnam_con*(vietnam_pres$historic + vietnam_pres[,(i+1)])
  d_con_vietnam <- (sum(d_top_con_vietnam))/(sum(d_bottom_con_vietnam))
  
  d_top_size_vietnam <- omega_vietnam_size*(vietnam_pres$historic - vietnam_pres[,(i+1)])
  d_bottom_size_vietnam <- omega_vietnam_size*(vietnam_pres$historic + vietnam_pres[,(i+1)])
  d_size_vietnam <- (sum(d_top_size_vietnam))/(sum(d_bottom_size_vietnam))
  
  d_top_equal_vietnam <- omega_vietnam_equal*(vietnam_pres$historic - vietnam_pres[,(i+1)])
  d_bottom_equal_vietnam <- omega_vietnam_equal*(vietnam_pres$historic + vietnam_pres[,(i+1)])
  d_equal_vietnam <- (sum(d_top_equal_vietnam))/(sum(d_bottom_equal_vietnam))
  
  def_vietnam[i,2] <- d_con_vietnam
  def_vietnam[i,3] <- d_size_vietnam
  def_vietnam[i,1] <- d_equal_vietnam
}


sabah_pres <- data.frame(historic=historic_occupancy_sabah,
                           sabah_all=current_occupancy_sabah,
                           dfr=current_dfr,
                           tfr=current_tfr,
                           kfr=current_kfr)

def_sabah <- matrix(0L, nrow = 4, ncol = 3)
def_sabah <- as.data.frame(def_sabah)
colnames(def_sabah) <- c("equal","conservation","size")
rownames(def_sabah) <- colnames(sabah_pres)[2:5]

for (i in 1:4) {
  d_top_con_sabah <- omega_sabah_con*(sabah_pres$historic - sabah_pres[,(i+1)])
  d_bottom_con_sabah <- omega_sabah_con*(sabah_pres$historic + sabah_pres[,(i+1)])
  d_con_sabah <- (sum(d_top_con_sabah))/(sum(d_bottom_con_sabah))
  
  d_top_size_sabah <- omega_sabah_size*(sabah_pres$historic - sabah_pres[,(i+1)])
  d_bottom_size_sabah <- omega_sabah_size*(sabah_pres$historic + sabah_pres[,(i+1)])
  d_size_sabah <- (sum(d_top_size_sabah))/(sum(d_bottom_size_sabah))
  
  d_top_equal_sabah <- omega_sabah_equal*(sabah_pres$historic - sabah_pres[,(i+1)])
  d_bottom_equal_sabah <- omega_sabah_equal*(sabah_pres$historic + sabah_pres[,(i+1)])
  d_equal_sabah <- (sum(d_top_equal_sabah))/(sum(d_bottom_equal_sabah))
  
  def_sabah[i,2] <- d_con_sabah
  def_sabah[i,3] <- d_size_sabah
  def_sabah[i,1] <- d_equal_sabah
}


defaunation <- rbind(def_vietnam,def_sabah)
rownames(defaunation) <- c("Hunted","Bach Ma NP","Saola NR","Xe Sap NPA","Degraded","Deramakot FR","Tangkulap FR","Kuamut FR")

df.long<-melt(defaunation)

df.long$site <- rownames(defaunation)

cols <- c("#005580",
          "#b37700",
          "#0088cc",
          "#4dc3ff",
          "#99ddff",
          "#a87602",
          "#edaf13",
          "#ffdc8f")
  
df.long$site <- factor(df.long$site, levels = c("Hunted", "Degraded", 
                                                "Bach Ma NP","Saola NR","Xe Sap NPA",
                                                "Tangkulap FR","Kuamut FR","Deramakot FR"))

df.long$what <- "site"
df.long$what[df.long$site=="Vietnam"] <- "region"
df.long$what[df.long$site=="Sabah"] <- "region"

p <- ggplot(df.long,aes(site,value,fill=site, width=0.5)) +
  geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
  geom_vline(xintercept=2.5, linetype="solid", color = "grey",size=1) +
  facet_grid(~variable, scales="free_x", space="free") +
  ylab(expression("Defaunation index")) +
  scale_fill_manual(values = cols) +
  labs(x = "", fill = "Extinction") +
  ylim(0,1) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.text.x=element_text(angle = 90, vjust = 0.4, hjust = 1),
        axis.title=element_text(size=14,face="bold"),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
        legend.position = "none",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold")
  )

pp <- p + facet_grid(~what, switch = "y")+
  theme(strip.text.y = element_text(size=14,angle = 180),
        strip.background = element_blank())


png(filename = "defaunation_fig_notranked.png", width=10*600, height=4*600, 
    res=600, bg="white")
p
dev.off()

