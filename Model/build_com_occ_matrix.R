#### setup detection history
#### by Jesse F. Abrams

library(gdata)
library(readxl)
library(camtrapR)
library(optimbase)
library(plyr)

wd <-("~/Dropbox (ScreenForBio)/Projects/Species_comparison")
setwd(wd)

site <- "Sabah_coarse"
#### load record table ####
rec.table <- read.csv("C:/Users/tilker/Dropbox (ScreenForBio)/Projects/Vietnam camera + leeches occupancy analysis/test_Jan_script/Sabah_Record_table_60min_deltaT_Combine_Coarse_modified.csv")
rec.table <- as.data.frame(rec.table)

#### load species list ####
species_list <- read.csv("C:/Users/tilker/Dropbox (ScreenForBio)/Projects/Vietnam camera + leeches occupancy analysis/test_Jan_script/Sabah_species_list.csv") 
View(species_list)

#### load CT table ####
CTtable <- read.csv("C:/Users/tilker/Dropbox (ScreenForBio)/Projects/Vietnam camera + leeches occupancy analysis/test_Jan_script/Sabah_DataTable_All_Coarse.csv")
CTtable <- as.data.frame(CTtable)
View(CTtable)


species_found <- unique(rec.table$Species)
num.species_found <- length(species_found)

species_found <- as.data.frame(species_found)

species.match <- (species_found$species_found %in% species_list$`Species`)
listed_species_found <- species_found[species.match,]
num.species <- length(listed_species_found)

camop <-  cameraOperation(CTtable = CTtable, 
                                stationCol = "Station", 
                                setupCol = "Date_setup",
                                retrievalCol = "Date_retrieval", 
                                hasProblems = FALSE,
                                dateFormat = "%m/%d/%Y", 
                                writecsv = TRUE, 
                                outDir = "C:/Users/tilker/Dropbox (ScreenForBio)/Projects/Vietnam camera + leeches occupancy analysis/test_Jan_script"
)


occasion <- 10

data.dummy <- rep(0, nrow(camop)*15*10)
site.det <- array(data.dummy, c(nrow(camop), 10, 15))  ## occasion length, number of species 

DateTime <- paste(rec.table$Date,rec.table$Time,sep=" ")
DateTime <- as.data.frame(DateTime)
new_record <- cbind(DateTime,rec.table)

for(i in 1:num.species) {
  spec <- as.character(listed_species_found[i])
  history_effort <- detectionHistory(recordTable          = rec.table,
                                     camOp                = camop,
                                     stationCol           = "Station",
                                     speciesCol           = "Species",
                                     recordDateTimeCol    = "DateTimeOriginal",
                                     recordDateTimeFormat = "%d/%m/%Y %H:%M",
                                     species              = spec,
                                     occasionLength       = occasion,
                                     day1                 = "station",
                                     maxNumberDays        = 100,
                                     occasionStartTime = 0,
                                     datesAsOccasionNames = FALSE,
                                     includeEffort        = TRUE,
                                     scaleEffort          = FALSE,
                                     timeZone             = "Asia/Saigon"#,
                                     #writecsv = TRUE
  )
  site.det[ , ,i] <- (history_effort$detection_history)
  
}

site.effort <- history_effort$effort
site.det <- aperm(site.det,c(3,1,2))

write.table(site.effort, file = paste(site,"_effort_matrix_mammal.csv", sep=""), sep = ",")

dput(site.det, paste(site,"_com_occ_detection_mammal.R", sep=""))
dput(site.effort, paste(site,"_com_occ_effort_mammal.R", sep=""))

species_out <- as.character(listed_species_found)
write.csv(species_out, paste(site,"_species_found_com_occ.csv", sep=""))



