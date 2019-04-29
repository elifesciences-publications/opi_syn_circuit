# This r script will read 'rawdataset.csv' and preprocess the data to get to a
# new dataframe that includes the effect amplitude of agonists/antagonists
# relative to baseline.

# Exclusion criteria are used to reject particular experiments. See below for
# explaination.

library(reshape2)
library(plyr)

df = read.csv("data/rawdataset.csv")
df$X = NULL # remove extra column with indexes
# make new column describing the part of the circuit measured.
df$circuit= paste(df$stimSource, df$recordLayer, df$cellType, sep='_')
df = df[, c(1:8,length(df),9:(length(df)-1))]
#rename instances of "eWill..." in animalID
df$animalID = gsub("eWill","Will",df$animalID)

# There is a set of animalIDs that are either called by date of experiment vs
# injection ID. To stratify this we replace all the date's with known injection
# ID to injectionID
df$animalID = gsub("Will20160328", "AV0113",df$animalID)
df$animalID = gsub("Will20160330", "AV0110",df$animalID)
df$animalID = gsub("Will20160830", "AV0165",df$animalID)
df$animalID = gsub("Will20161005", "AV0179",df$animalID)
df$animalID = gsub("Will20161212", "AV0216",df$animalID)
df$animalID = gsub("Will20161213", "AV0217",df$animalID)

# Criteria:
# =========

#1. When onsetTime of EPSC is between 0 and 1 ms: turn onsetTime and onset2Peak
# to NaN. The detection of the onsetTime is corrupted by noise in those
# traces.
#2. Turn all values from EPSC conditions where EPSC amplitude <= -1 nA into NaN.
# Those instances will most likely be biased by failure to properly clamp.
# So calculation of agonist/antagonist effects may not be accurate.
#3. Turn all values from IPSC conditions where IPSC amplitude >= 1 nA & IPSC
# and halfwidth >= thresholdHWIPSC into NaN. thresholdHWIPSC is mean + 3*SD
# determined from PV;Ai32 IPSCs. The latter is a pure population of GABAergic
# synapses. Similar to criteria #2 signals above 1nA will most likely be biased
# by failure to properly clamp, so the calculation of agonist/antagonist effects
# may not be accurate

# Calculate halfwidth treshold for IPSCs based on the IPSC from PV;Ai32
# recordings
hwPvIPSC = subset(df,df$stimSource == 'PV')
meanHw = mean(hwPvIPSC$halfwidth, na.rm = TRUE)
sdHw = sd(hwPvIPSC$halfwidth, na.rm = TRUE)
threshHWIPSC = meanHw + 3 * sdHw

# Clean-up based on criteria
# Criteria #1
df$onset[df$signal == 'EPSC' & df$onset <= 1] = NaN
df$onset2peakTime[df$signal == 'EPSC' & df$onset <= 1] = NaN

# Criteria #2
df[df$signal == 'EPSC' & !is.na(df$amplitude) & df$amplitude <= -1000,
    13:length(df)] = NaN

# Criteria #3
df[df$signal == 'IPSC' & !is.na(df$amplitude) & df$amplitude <= 1000 &
    df$halfwidth >= threshHWIPSC, 13:length(df)] = NaN

#One exception of non-generalized criteria but experiment that should be
#excluded a paired pulse was given for cell BJ0576. The IPSC grew so big that
# the first IPSC invaded the second, which led to an enormous decay.
df$decay[df$cellID == 'BJ0576' & df$conditionName == 'DAMGOCTAP'] = NaN
# will20171027c002 has two baselines, Take out the first baseline
df[(df$cellID=='will20171027c002' & df$waveformSequence==1),] = NA

#Add a variable: charge(area)/amplitude
df$chargeOAmplitude = df$area/df$amplitude
df = df[,c(1:22,ncol(df),24:(ncol(df)-1))]

# make sure that "DPDPE" is correctly written.
df$conditionName= gsub("DPDE","DPDPE", df$conditionName)
df$conditionName= gsub("BASEDPDE","BASEDPDPE", df$conditionName)
df$conditionName= gsub("DPDEICI","DPDPEICI", df$conditionName)
df$conditionName= gsub("DPDPENTD","DPDPENALTRI", df$conditionName)
df$conditionName= gsub("DPDENALTRI","DPDPENALTRI", df$conditionName)
df$conditionName= gsub("DAMGONLX","DAMGONALOX",df$conditionName)

cols = colnames(df)
mdf=melt(df, id.vars = cols[1:12], measure.vars = cols[13:23])
#take the major 3 agonist and append column drugGroup
dfENK = mdf[grep("Enk", mdf$conditionName, perl=TRUE, ignore.case=TRUE),]
dfENK$drugGroup = 1
dfDAMGO = mdf[grep("DAMGO", mdf$conditionName, perl=TRUE, ignore.case=TRUE),]
dfDAMGO$drugGroup = 2
dfDPDPE = mdf[grep("DPDPE", mdf$conditionName, perl=TRUE, ignore.case=TRUE),]
dfDPDPE$drugGroup = 3
dfDNQX = mdf[grep("DNQX", mdf$conditionName, perl=TRUE, ignore.case=TRUE),]
dfDNQX$drugGroup = 4
dfNBQX = mdf[grep("NBQX", mdf$conditionName, perl=TRUE, ignore.case=TRUE),]
dfNBQX$drugGroup = 5

#re-combine the subdataframes
rmdf=rbind(dfENK, dfDAMGO, dfDPDPE, dfDNQX, dfNBQX)
dfBase = rmdf[grep("Base",rmdf$conditionName, perl=TRUE, ignore.case=TRUE),]
dfAgonist = subset(rmdf, conditionName == "ENK"| conditionName == "DAMGO"|
    conditionName == "DPDPE"| conditionName == "DNQX"| conditionName == "NBQX")
dfAntagonist = subset(rmdf, conditionName == "ENKWASH"|
    conditionName == "ENKNALOX"| conditionName == "DAMGOCTAP"|
    conditionName == "DPDPENALTRI"|conditionName == "DPDPEICI"|
    conditionName == "DPDPENALOX"|conditionName == "DAMGONALOX"|
    conditionName == "DPDPENLX"|conditionName == "DAMGONLX"|
    conditionName == "DPDPENTD")
#combine the base, agonist dataframes
dfDrug = merge(dfBase, dfAgonist, by=c("cellID","animalID","genoType",
    "cellType","recordLayer","recordRegion","stimSource","stimChannel",
    "circuit","signal","samplerate","variable","drugGroup"), all.x=TRUE)
    #the merge puts dfBase unique cols as .x and dfAgonist unique cols as .y,
    # so you need to rename those cols
dfDrug = rename(dfDrug,c("conditionName.x"="baseName", "value.x"="baseValue",
    "conditionName.y"="agonistName","value.y"="agonistValue"))
#combine the dfDrug with dfAntagonist
dfDrug = merge(dfDrug, dfAntagonist, by=c("cellID","animalID","genoType",
    "cellType","recordLayer","recordRegion","stimSource","stimChannel",
    "circuit","signal","samplerate","variable","drugGroup"), all.x=TRUE)
#rename added cols
dfDrug = rename(dfDrug,c("conditionName"="antagonistName",
    "value"="antagonistValue"))

# Next edit entries that have a opposite sign in agonist or antagonist value
# relative to baseline. Cases with sign conversion will get the value 0.00001 or
# -0.00001 depending on the sign of the baseline. This when e.g. a IPSC is
# completely blocked. If the voltage clamp command is not perfectly set at
# sodium reversal potential this will result in the appearance of a small inward
# current. This would make the signal switch from positive float to negative
# float
countAll = 0
countReversals = 0
countIPSCAmpl = 0
for (i in 1:nrow(dfDrug)){
  countAll = countAll + 1
    if(is.null(dfDrug$parameter[i])==FALSE && is.na(dfDrug$parameter[i])==FALSE){
      if(dfDrug$signal[i] == "IPSC" && dfDrug$parameter[i]== "amplitude"){
      countIPSCAmpl = countIPSCAmpl + 1}
    }
    if(is.na(dfDrug$baseValue[i])==FALSE){
        if(is.na(dfDrug$agonistValue[i])==FALSE){
            if(dfDrug$baseValue[i] < 0 & dfDrug$agonistValue[i] > 0){
                dfDrug$agonistValue[i] = -0.00001
                if(is.null(dfDrug$parameter[i])==FALSE && is.na(dfDrug$parameter[i])==FALSE){
                  if(dfDrug$signal[i] == "IPSC" && dfDrug$parameter[i]== "amplitude"){
                  countReversals = countReversals +1}
                }
            }
            else if(dfDrug$baseValue[i] > 0 & dfDrug$agonistValue[i] < 0){
                dfDrug$agonistValue[i] = 0.00001
                if(is.null(dfDrug$parameter[i])==FALSE && is.na(dfDrug$parameter[i])==FALSE){
                  if(dfDrug$signal[i] == "IPSC" && dfDrug$parameter[i]== "amplitude"){
                    countReversals = countReversals +1}
                  }
            }
        }
        if(is.na(dfDrug$agonistValue[i]) == TRUE){
            if(dfDrug$baseValue[i] > 0){
                dfDrug$agonistValue[i] = 0.00001
            }
            else if(dfDrug$baseValue[i] < 0){
                dfDrug$agonistValue[i] = -0.00001
            }
        }
        if(is.na(dfDrug$antagonistValue[i])==FALSE){
              if(dfDrug$baseValue[i] < 0 & dfDrug$antagonistValue[i] > 0){
                  dfDrug$antagonistValue[i] = -0.00001
              }
              else if(dfDrug$baseValue[i] > 0 & dfDrug$antagonistValue[i] < 0){
                  dfDrug$antagonistValue[i] = 0.00001
              }

        }
    }
}

# Calculate effect relative to baseline for agonist and antagonist
dfDrug = transform(dfDrug,
    agonistEffect = dfDrug$agonistValue/dfDrug$baseValue*100)
dfDrug = transform(dfDrug,
    antagonistEffect = dfDrug$antagonistValue/dfDrug$baseValue*100)
dfDrug = rename(dfDrug,c("variable"="parameter"))
colsDrug = colnames(dfDrug)
mdfDrug = melt(dfDrug, id.vars = colsDrug[1:19], measure.vars = colsDrug[20:21])
mdfDrug$varName = paste(mdfDrug$parameter, mdfDrug$variable, sep='_')
# Insert extra row containing tailratio information
tDf = read.csv('data/tailQuantification.csv')
colnames(tDf) = c('animalID','tailRatio')
dfDrug$tailRatio = tDf$tailRatio[match(dfDrug$animalID, tDf$animalID)]

#Append latency spike data and some extra experiment performed by Birdsong to
# the dataset
apdf=read.csv('data/latencySpikePeak.csv')
appendDf = read.csv('data/data_appendMSN.csv')
appendDf$X = NULL
appendDf$agonistEffect = appendDf$agonistEffect*100
appendDf$antagonistEffect = appendDf$antagonistEffect*100
#Add the appendDf and save result as .csv
df = rbind(dfDrug, appendDf)
rm(appendDf)

# Save the dataframe
write.csv(df,file= "data/effectdataset.csv" )
