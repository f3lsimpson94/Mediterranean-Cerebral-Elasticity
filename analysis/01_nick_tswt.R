### --------------------------------------------------------------
### Behavioral preprocessing pipeline for TSWT, PCC, EEG&OPT SESSION
### --------------------------------------------------------------

### --------------------------------------------------------------
### Load libraries
### --------------------------------------------------------------

library(data.table)
library(ggplot2)
library(readr)
library(plyr)
library(dplyr)
library(doBy)
library(writexl)
library(readxl)

rm(list = ls())
options(scipen=999)

### --------------------------------------------------------------
### Check IDs for EEG session 
### --------------------------------------------------------------
setwd("F:\\TSWT\\NC_36\\") #set working directory to the data folder
paths <- dir(pattern = "\\DATA.dat$") #import names of the files
names(paths) <- basename(paths)
IDtable <- as.data.frame(unique(paths))

paths1 <- dir(pattern = "\\DATA_1.dat$") #some participants have two output files -- check them with Binderinfo to select the correct one
names(paths1) <- basename(paths1)
IDtable1 <- as.data.frame(unique(paths1))

#1007 has one file with "1" only

### --------------------------------------------------------------
### Import/load the data for EEG
### --------------------------------------------------------------

rm(list = ls())
options(scipen=999)

# ## Formula to load new data
setwd("F:\\TSWT\\NC_36\\")
paths <- dir(pattern = "\\.dat$")
names(paths) <- basename(paths)
all.eeg <- ldply(paths, read.delim)
all.eeg <- all.eeg[grepl("DATA", all.eeg$`.id`),]
colnames(all.eeg)[1] = "ID"
str(all.eeg)
head(all.eeg)
names(all.eeg)
save(all.eeg, file = "F:\\TSWT\\NC_36\\all.eeg.RDat")
load("F:\\TSWT\\NC_36\\all.eeg.RDat")
setwd("F:\\TSWT\\NC_36\\") #set your working directory

eeg = all.eeg

### Fix ID issues (here I'm removing invalid files)
unique(eeg$ID)
eeg$ID[eeg$ID == "1007_Phase1_T_DATA_1.dat"] <-  "1007_Phase1_T_DATA.dat"

unique(eeg$SubjectID)
eeg$SubjectID[eeg$SubjectID == "1002-1"] <-  "1002"
eeg$SubjectID[eeg$SubjectID == "1006_1"] <-  "1006"

### --------------------------------------------------------------
### Reduce the columns and decode events
### --------------------------------------------------------------
eeg = eeg[,c(5,13:18,11,24:25,2:4,19:23)]
eeg = eeg[,c(1,12,13,11,9,2,3,4,5,6,7,8,10,14:18)]

# Block type
colnames(eeg)[12] <- "BlockType"
eeg$BlockType[eeg$BlockType == 5] <- "single"
eeg$BlockType[eeg$BlockType == 6] <- "single"
eeg$BlockType[eeg$BlockType == 7] <- "mixed"

# Trial type
eeg$TT[eeg$BlockType == "mixed" & eeg$TrialType == "S"] <- "MS"
eeg$TT[eeg$BlockType == "mixed" & eeg$TrialType == "R"] <- "MR"
eeg$TT[eeg$BlockType == "single" & eeg$TrialType == "R"] <- "AR"
eeg$TT[eeg$TrialNumPerBlock == 1] <- "warmup"
eeg$TT[eeg$TrialNumPerBlock == 2] <- "warmup"

# Task type
eeg$TaskName[eeg$TaskName == "Vowel/Consonant"] <- "Letter"
eeg$TaskName[eeg$TaskName == "Odd/Even"] <- "Digit"

# Task type transition
eeg$ABTrans[eeg$ABTrans == "AA"] <- "DD" #digit-digit
eeg$ABTrans[eeg$ABTrans == "BB"] <- "LL" #letter-letter
eeg$ABTrans[eeg$ABTrans == "AB"] <- "DL" #digit-LETTER
eeg$ABTrans[eeg$ABTrans == "BA"] <- "LD" #letter-DIGIT (digit in the current trial)
eeg$ABTrans[eeg$ABTrans == "~A"] <- "~D" #first digit trial in a block
eeg$ABTrans[eeg$ABTrans == "~B"] <- "~L" #first letter trial in a block

# Response hand
eeg$Resp[eeg$Resp == "1"] <- "Left"
eeg$Resp[eeg$Resp == "2"] <- "Right"
eeg$Resp[eeg$Resp == "3"] <- "NoResponse"

# Congruency
eeg$Congruency[eeg$Congruency == "N"] <- "Neu"
eeg$Congruency[eeg$Congruency == "I"] <- "Inc"

# #Equalize number of trials 
# singleLetter <- 
#   subset(eeg, TT %in% c('AR') & TaskName == "letter") %>%
#   group_by(SubjectID) %>%
#   mutate(Seq=1:n())
# singleLetter <- singleLetter[singleLetter$Seq < 121,]
# singleDigit <- 
#   subset(eeg, TT %in% c('AR') & TaskName == "digit") %>%
#   group_by(SubjectID) %>%
#   mutate(Seq=1:n())
# singleDigit <- singleDigit[singleDigit$Seq < 121,]
# eeg <- eeg[!eeg$TT == "AR",]
# eeg <- rbind(eeg, singleDigit[,-c(20)])
# eeg <- rbind(eeg, singleLetter[,-c(20)])
# eeg <- eeg %>% arrange(SubjectID, TrialNum)

# ID as a factor
eeg$SubjectID = as.factor(eeg$SubjectID)


#Previous congruency
eeg$TT1 = lag(eeg$Congruency)
eeg$TT1 <- ifelse(eeg$Congruency == eeg$TT1, "R", "S")

### --------------------------------------------------------------
### Data trimming
### --------------------------------------------------------------

# Compute total number of trials
total <- reshape2::dcast(eeg, SubjectID ~ TT, length, value.var = c("TT"), drop = FALSE)

# Exclude warm-ups -- 10 for every participant, see "total" above
eeg <- eeg[!eeg$TT == "warmup",] 
eeg$TT = as.factor(eeg$TT)


# Exclude misses 
misses <-  reshape2::dcast(eeg[eeg$Accuracy == 3,], SubjectID ~ TT, length, value.var = c("TT"), drop = FALSE) #misses are coded as "3"
eeg <- eeg[!eeg$Accuracy == 3,]

# Exclude post-error trials and post-misses
perr <- reshape2::dcast(eeg, SubjectID ~ LastTrialErr + TT, length, value.var = c("TT"), drop = FALSE) 
eeg <- eeg[!eeg$LastTrialErr == 1,]

# Too fast exclusion
fast <- reshape2::dcast(eeg[eeg$ReactionTime <= 200,], SubjectID ~ TT, length, value.var = c("TT"), drop = FALSE)
eeg <- eeg[!eeg$ReactionTime <= 200,]

# Too slow -- +-3 SD from mean (TT x Conruency)
meanRT = reshape2::dcast(eeg, SubjectID + TT + Congruency ~ ., mean, value.var = c("ReactionTime"))
colnames(meanRT)[4] <- "meanRT"
sdRT = reshape2::dcast(eeg, SubjectID + TT + Congruency ~ ., sd, value.var = c("ReactionTime"))
colnames(sdRT)[4] <- "sdRT"
slow <- merge(meanRT, sdRT, by = c("SubjectID", "TT", "Congruency"), all = TRUE)
eeg = merge(eeg, slow, by = c("SubjectID", "TT", "Congruency"))
slow_summary <- eeg[eeg$ReactionTime >= eeg$meanRT+3*eeg$sdRT,]
slow <-  reshape2::dcast(slow_summary, SubjectID ~ TT + Congruency, length, value.var = c("TT"), drop = FALSE)
eeg <- eeg[!eeg$ReactionTime >= eeg$meanRT+3*eeg$sdRT,]
eeg <- eeg[,-c(21:22)]
rm(sdRT, meanRT, slow_summary, all.eeg)

# Valid trials
valid = reshape2::dcast(SubjectID ~ TT + Congruency, length, data = eeg, value.var = c("TT"), drop = FALSE)

### --------------------------------------------------------------
### Summarize data trimming 
### --------------------------------------------------------------

colnames(total) <- c("SubjectID", "Total_AR", "Total_MR", "Total_MS", "Total_warmup")
colnames(misses) <- c("SubjectID", "Misses_AR", "Misses_MR", "Misses_MS")
perr = perr[,c(1,5:7)]
colnames(perr) <- c("SubjectID", "Perr_AR", "Perr_MR", "Perr_MS")
colnames(fast) <- c("SubjectID", "Fast_AR", "Fast_MR", "Fast_MS")
colnames(slow) <- c("SubjectID", "Slow_AR_Inc", "Slow_AR_Neu", "Slow_MR_Inc", "Slow_MR_Neu", "Slow_MS_Inc", "Slow_MS_Neu")
colnames(valid) <- c("SubjectID", "Valid_AR_Inc", "Valid_AR_Neu", "Valid_MR_Inc", "Valid_MR_Neu", "Valid_MS_Inc", "Valid_MS_Neu")

preprocInfo = join_all(list(total, misses, perr, slow, valid), by = "SubjectID", type = 'full')

total_long = reshape2::melt(total, id.vars = "SubjectID")
summaryBy(value ~ variable , data = total_long,
          FUN = function(x) { c(n = length(x), mean = mean(x), sd = sd(x), min = min(x), max = max(x)) } )

perr_long = reshape2::melt(perr, id.vars = "SubjectID")
summaryBy(value ~ variable , data = perr_long,
          FUN = function(x) { c(n = length(x), mean = mean(x), sd = sd(x), min = min(x), max = max(x)) } )

misses_long = reshape2::melt(misses, id.vars = "SubjectID")
summaryBy(value ~ variable , data = misses_long,
          FUN = function(x) { c(n = length(x), mean = mean(x), sd = sd(x), min = min(x), max = max(x)) } )

slow_long = reshape2::melt(slow, id.vars = "SubjectID")
summaryBy(value ~ variable , data = slow_long,
          FUN = function(x) { c(n = length(x), mean = mean(x), sd = sd(x), min = min(x), max = max(x)) } )

valid_long = reshape2::melt(valid, id.vars = "SubjectID")
summaryBy(value ~ variable , data = valid_long,
          FUN = function(x) { c(n = length(x), mean = mean(x), sd = sd(x), min = min(x), max = max(x)) } )


ggplot(perr_long[perr_long$variable == "Perr_AR" | perr_long$variable == "Perr_MR" | perr_long$variable == "Perr_MS",], aes(x = value, fill = variable)) + geom_histogram(color='#e9ecef', alpha=0.6, position='identity') + facet_grid(variable ~ .) + labs(title = "Counting post-incorrect-response", x = "number of trials", y = "number of participants") + theme_classic(base_size = 14)
ggplot(misses_long, aes(x = value, fill = variable)) + geom_histogram(color='#e9ecef', alpha=0.6, position='identity')  + facet_grid(variable ~ .) + labs(title = "Counting misses", x = "number of trials", y = "number of participants") + theme_classic(base_size = 14)
ggplot(slow_long, aes(x = value, fill = variable)) + geom_histogram(color='#e9ecef', alpha=0.6, position='identity') + facet_grid(variable ~ .) + labs(title = "Counting slow trials", x = "number of trials", y = "number of participants") + theme_classic(base_size = 14)
ggplot(valid_long, aes(x = value, fill = variable)) + geom_histogram(color='#e9ecef', alpha=0.6, position='identity') + facet_grid(variable ~ .) + labs(title = "Counting valid trials", x = "number of trials", y = "number of participants") + theme_classic(base_size = 14)

notes <- c(
  "Behavioral preprocessing for EEG session",
  "Data trimming: misses (no response), post-errors (including post-misses), too fast (<= 200 ms), too slow (+-3SD from Mean for TrialType x Congruency)",
  "Code written in R by PK")
notes <- as.data.frame(notes)


### --------------------------------------------------------------
### Save the data
### --------------------------------------------------------------

#Histogram
hist(eeg$ReactionTime, labels = T)

#Accuracy
eeg$Accuracy[eeg$Accuracy == 2] <- 0

#Save preprocessed data
sheets <- list("Data" = eeg, "PreprocessingSummary" = preprocInfo, "Notes" = notes)
write_xlsx(sheets, "F:\\TSWT\\NC_36\\TSWT_EEG_NC_36M.xlsx")


### --------------------------------------------------------------
### Compute Costs
### --------------------------------------------------------------


EEG <- read_excel("F:\\TSWT\\NC_36\\TSWT_EEG_NC_36M.xlsx", sheet = "Data")

### RT
eeg.rt.congr = data.table::dcast(setDT(EEG[EEG$Accuracy == 1,]), SubjectID ~ Congruency, value.var = c("ReactionTime"), mean)
eeg.rt.congr$congruencyCost = eeg.rt.congr$Inc - eeg.rt.congr$Neu

eeg.rt.tt = data.table::dcast(setDT(EEG[EEG$Accuracy == 1,]), SubjectID ~ TT, value.var = c("ReactionTime"), mean)
eeg.rt.tt$mixCost = eeg.rt.tt$MR - eeg.rt.tt$AR
eeg.rt.tt$swCost = eeg.rt.tt$MS - eeg.rt.tt$MR

eeg.rt = merge(eeg.rt.congr, eeg.rt.tt, by = c("SubjectID"))
rm(eeg.rt.congr, eeg.rt.tt)

### ERR
eeg.acc.congr = data.table::dcast(setDT(EEG), SubjectID ~ Congruency, value.var = c("Accuracy"), mean)
eeg.err.congr <- eeg.acc.congr %>% mutate(Inc = 1 - Inc, Neu = 1 - Neu)
eeg.err.congr$Inc[ eeg.err.congr$Inc > .50] <- NA
eeg.err.congr$Neu[ eeg.err.congr$Neu > .50] <- NA
eeg.err.congr$congruencyCost = eeg.err.congr$Inc - eeg.err.congr$Neu

eeg.acc.tt = data.table::dcast(setDT(EEG), SubjectID ~ TT, value.var = c("Accuracy"), mean)
eeg.err.tt <- eeg.acc.tt %>% mutate(AR = 1 - AR, 
                                    MR = 1 - MR,
                                    MS = 1 - MS)
eeg.err.tt$AR[ eeg.err.tt$AR > .50] <- NA
eeg.err.tt$MR[ eeg.err.tt$MR > .50] <- NA
eeg.err.tt$MS[ eeg.err.tt$MS > .50] <- NA
eeg.err.tt$mixCost = eeg.err.tt$MR - eeg.err.tt$AR
eeg.err.tt$swCost = eeg.err.tt$MS - eeg.err.tt$MR

eeg.err = merge(eeg.err.congr, eeg.err.tt, by = c("SubjectID"))
rm(eeg.err.congr, eeg.err.tt, eeg.acc.congr, eeg.acc.tt)



# Get all column names
col_names <- colnames(eeg.rt)
# Add prefix "rt_" to each column name except for "SubjectID"
colnames(eeg.rt) <- ifelse(col_names != "SubjectID", paste0("RT_", col_names), col_names)
# Get all column names
col_names <- colnames(eeg.err)
# Add prefix "rt_" to each column name except for "SubjectID"
colnames(eeg.err) <- ifelse(col_names != "SubjectID", paste0("ERR_", col_names), col_names)


task = merge(eeg.rt, eeg.err, id.vars = c("SubjectID"), all = TRUE)
colnames(task)[1] <- "ID"
task$ID = as.character(task$ID)

# ID_vector <- read.csv("//bi-cnl-nas3/data/pcc/matlab_scripts/TSWT_scripts/NEWCASTLE EEG 2023/4ANALYSES/IDlist.txt", sep="")$SubjectID
# task <- task[task$ID %in% ID_vector,]

task$IES_Inc <- task$RT_Inc/(1- task$ERR_Inc)
task$IES_Neu <- task$RT_Neu/(1- task$ERR_Neu)
task$IES_congruencyCost = task$IES_Inc - task$IES_Neu

task$IES_AR <- task$RT_AR/(1- task$ERR_AR)
task$IES_MR <- task$RT_MR/(1- task$ERR_MR)
task$IES_MS <- task$RT_MS/(1- task$ERR_MS)
task$IES_mixCost = task$IES_MR - task$IES_AR
task$IES_swCost = task$IES_MS - task$IES_MR

rm(EEG, eeg.err, eeg.rt)

write_xlsx(task, "F:\\TSWT\\NC_36\\TSWT_EEG_NC_36M_TASK.xlsx")
