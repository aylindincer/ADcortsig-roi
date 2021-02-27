#Setup Workspace
#===============================================================================

# Clear workspace
rm(list = ls()) 

#Libraries
library(ggplot2)
library(dplyr)
library(data.table)
library(pROC)

#===============================================================================

#Setup Analysis Parameters
#===============================================================================
ADRCcohort <- 1
DIANcohort <- 0

ADRCmask <- 0
DIANmask <- 1

IndividROI <- 0
CortSigROI <- 1

#Setup directory paths based on paramaters selection
Cohort <- ifelse(DIANcohort==1,"DIAN-analysis-qcache","ADRC-analysis-qcache")
Cohortmap <- ifelse(DIANmask==1,"DIAN-mask","ADRC-mask")
ROImap <- ifelse(IndividROI==1,"IndividROI","CortSigROI")
#===============================================================================

# Weighted Average of LH and RH Optimal Monte Carlo Thresholds
#===============================================================================
#Function to calculate waverage and write csv
if(CortSigROI==1){
  ThickSumMeas <- function(MCTH_LH, MCTH_RH, inputpath, outputpath){
    #Read thickness files from specified MC thresholds
    LHThickpath <- list.files(path=inputpath,full.names=TRUE,
                              pattern=(paste0("^lh\\..*th",MCTH_LH,".*\\.csv$")))
    RHThickpath <- list.files(path=inputpath,full.names=TRUE,
                              pattern=(paste0("^rh\\..*th",MCTH_RH,".*\\.csv$")))
    LHThick <- read.csv(LHThickpath,header = TRUE,stringsAsFactors = FALSE)
    RHThick <- read.csv(RHThickpath,header = TRUE,stringsAsFactors = FALSE)
    
    #Diffentiate between LH and RH
    colnames(LHThick)[2:14] <- paste("LH", colnames(LHThick)[2:14], sep = "_")
    colnames(RHThick)[2:14] <- paste("RH", colnames(RHThick)[2:14], sep = "_")
    concatThick <- data.frame(full_join(LHThick, RHThick, by = "Session"))
    
    #Calculate Weighted Average for each participant
    concatThick$ThickAvg <- ((concatThick$LH_ThickAvg*concatThick$LH_NumVert) + 
                             (concatThick$RH_ThickAvg*concatThick$RH_NumVert))/
                            (concatThick$LH_NumVert+concatThick$RH_NumVert)
    filestr <- (gsub("^lh\\.", "", basename(LHThickpath)))
    filestr <- (gsub("\\.th.*", "", basename(filestr)))
    write.csv(concatThick, file = paste0(outputpath,filestr,".SumMeas.csv"), 
              row.names = FALSE)
  }
  
  #DIAN mask optimal Monte Carlo thresholds
  if(DIANmask==1){
    MCTH_LH <- 40
    MCTH_RH <- 30
  }
  #ADRC mask optimal Monte Carlo thresholds
  if(ADRCmask==1){
    MCTH_LH <- 23
    MCTH_RH <- 13
  }
  inputpath <- paste0(Cohort,"/",Cohortmap,"/IndividROI/data/")
  outputpath <- paste0(Cohort,"/",Cohortmap,"/CortSigROI/data/")
  
  #Run ThickSumMeas function
  ThickSumMeas(MCTH_LH, MCTH_RH, inputpath, outputpath)
}
#===============================================================================

#Define groups
#===============================================================================
#DIAN groups
if(DIANcohort==1){
  preclinicalADAD <- read.csv("DIAN.preclinicalScan.csv",
                          header=TRUE,stringsAsFactors=FALSE)
  preclinicalADAD <- preclinicalADAD[c(1:1556)] #keep columns consistent
  preclinicalADADQ4 <- read.csv("DIAN.preclinicalScanQ4.csv",
                            header=TRUE,stringsAsFactors=FALSE)
  preclinicalADADQ4 <- preclinicalADADQ4[c(1:1556)] #keep columns consistent
  controls <- read.csv("DIAN.controls.csv",
                       header=TRUE,stringsAsFactors=FALSE)
  
  #Preclinical + Control Groups
  PreclinicalADADControls <- data.frame(rbind(controls, preclinicalADAD))
  colnames(PreclinicalADADControls)[1] <- "Session"
  PreclinicalADADControls$ABpositive_Scan_num <- ifelse(PreclinicalADADControls$ABpositive_Scan==TRUE,1,0)
  PreclinicalADADControls$Age <- PreclinicalADADControls$visitage
  PreclinicalADADControls$Gender <- PreclinicalADADControls$CLIN_gender_cat
  colnames(PreclinicalADADControls)[colnames(PreclinicalADADControls)=='totalHippocampus'] <- 'tHCV'
  colnames(PreclinicalADADControls)[colnames(PreclinicalADADControls)=='PET_PET_fSUVR_rsf_TOT_CORTMEAN'] <- 'AmyloidCont'
  
  #Preclinical (top 25%) + Control Groups
  PreclinicalADADQ4Controls <- data.frame(rbind(controls, preclinicalADADQ4))
  colnames(PreclinicalADADQ4Controls)[1] <- "Session"
  PreclinicalADADQ4Controls$ABpositive_Scan_num <- ifelse(PreclinicalADADQ4Controls$ABpositive_Scan==TRUE,1,0)
  PreclinicalADADQ4Controls$Age <- PreclinicalADADQ4Controls$visitage
  PreclinicalADADQ4Controls$Gender <- PreclinicalADADQ4Controls$CLIN_gender_cat
  colnames(PreclinicalADADQ4Controls)[colnames(PreclinicalADADQ4Controls)=='totalHippocampus'] <- 'tHCV'
  colnames(PreclinicalADADQ4Controls)[colnames(PreclinicalADADQ4Controls)=='PET_PET_fSUVR_rsf_TOT_CORTMEAN'] <- 'AmyloidCont'
  
}

#ADRC groups
if(ADRCcohort==1){
  preclinicalAD <- read.csv("ADRC.preclinicalScan.csv",
                          header=TRUE,stringsAsFactors=FALSE)
  preclinicalAD <- preclinicalAD[c(2:296)] #remove misc column
  preclinicalADQ4 <- read.csv("ADRC.preclinicalScanQ4Cent.csv",
                            header=TRUE,stringsAsFactors=FALSE)
  preclinicalADQ4 <- preclinicalADQ4[c(2:296)] #remove misc column
  controls <- read.csv("ADRC.controls.csv",
                       header=TRUE,stringsAsFactors=FALSE)
  controls <- controls[c(2:296)] #remove misc column

  #Preclinical + Control Groups
  PreclinicalADControls <- data.frame(rbind(controls, preclinicalAD))
  colnames(PreclinicalADControls)[1] <- "Session"
  PreclinicalADControls$ABpositive_Scan_num <- ifelse(PreclinicalADControls$ABpositive_Scan==TRUE,1,0)
  colnames(PreclinicalADControls)[colnames(PreclinicalADControls)=='totalHippocampus'] <- 'tHCV'
  colnames(PreclinicalADControls)[colnames(PreclinicalADControls)=='centiloid'] <- 'AmyloidCont'
  
  #Preclinical (top 25%) + Control Groups
  PreclinicalADQ4Controls <- data.frame(rbind(controls, preclinicalADQ4))
  colnames(PreclinicalADQ4Controls)[1] <- "Session"
  PreclinicalADQ4Controls$ABpositive_Scan_num <- ifelse(PreclinicalADQ4Controls$ABpositive_Scan==TRUE,1,0)
  colnames(PreclinicalADQ4Controls)[colnames(PreclinicalADQ4Controls)=='totalHippocampus'] <- 'tHCV'
  colnames(PreclinicalADQ4Controls)[colnames(PreclinicalADQ4Controls)=='centiloid'] <- 'AmyloidCont'
  
}
#===============================================================================

#Combine ROI cortical thickness with clinical and PET data. (cleared)
#===============================================================================
#Function to combine data
concatdata <- function(mygroup, thicknessfile, mysavepath){
  mygroupstr <- deparse(substitute(mygroup))
  for(i in 1:length(thicknessfile)){
    thickness <- read.csv(paste0(thicknessfile[i]), header = TRUE, 
                          stringsAsFactors = FALSE)
    thickness <- thickness[thickness$Session %in% mygroup$Session ,]
    thickness <- data.frame(full_join(thickness, mygroup, by = "Session"))
    filestr <- (gsub(".concat.csv", "", basename(thicknessfile[i])))
    write.csv(thickness, file = paste0(mysavepath,mygroupstr,"/data_groups/",
                                       mygroupstr,".",filestr,".csv"),
              row.names = FALSE)
    
  }
}

#Pull thickness data + Setup save directory path
thicknessfile <- list.files(path=paste0(Cohort,"/",Cohortmap,"/",ROImap,"/data/"),
                              full.names=TRUE,pattern=".csv")
mysavepath <- paste0(Cohort,"/",Cohortmap,"/",ROImap,"/")
  
#Run concatdata function
if(DIANcohort==1){
  concatdata(PreclinicalADADControls, thicknessfile, mysavepath)
  concatdata(PreclinicalADADQ4Controls, thicknessfile, mysavepath)
}

if(ADRCcohort==1){
  concatdata(PreclinicalADControls, thicknessfile, mysavepath)
  concatdata(PreclinicalADQ4Controls, thicknessfile, mysavepath)
}
#===============================================================================

#Statistical Analysis Paramaters (cleared)
#===============================================================================
#Set up directory path for analysis
startpath <- paste0(Cohort,"/",Cohortmap,"/",ROImap,"/")

#Linear model functions
lmfxn <- c(AmyloidCont~ThickAvg+Age+Gender,
           AmyloidCont~ThickAvg+Age+Gender+tHCV)

#Logisitic regression functions
logfxn <- c(ABpositive_Scan_num~ThickAvg+Age+Gender,
            ABpositive_Scan_num~ThickAvg+Age+Gender+tHCV)

#ROC curve functions
throcfxn <- ABpositive_Scan_cat~ThickAvg
hcvrocfxn <- ABpositive_Scan_cat~tHCV
 
#Adjusted plot (lm), age plot AES for ggplot, pib plot AES for ggplot
adjfxn1 <- AmyloidCont~Age+Gender
adjfxn2 <- ThickAvg~Age+Gender
ageplotaes <- aes(Age,ThickAvg)

if(DIANcohort==1){
  amyloidplotaes <- aes(AmyloidCont,ThickAvg,color=ABpositive_Scan_cat)
}

if(ADRCcohort==1){
  amyloidplotaes <- aes(PET_fSUVR_rsf_TOT_CORTMEAN,ThickAvg,
                        color=ABpositive_Scan_cat)
  centplotaes <- aes(AmyloidCont,ThickAvg,color=ABpositive_Scan_cat,
                     shape = Tracer)
}
#===============================================================================

#*******************************************************************************
#Statistical Analysis Master Function
#===============================================================================
statanalysis <- function(mygroup,startpath,ADRCcohort,DIANcohort,IndividROI,
                         lmfxn,logfxn,throcfxn,hcvrocfxn,adjfxn1,adjfxn2,
                         ageplotaes,amyloidplotaes,centplotaes=NULL){
  #Convert mygroup variable to string for future directory paths
  groupstr <- deparse(substitute(mygroup))
  grouppath <- list.files(path=paste0(startpath,groupstr,"/data_groups/"),
                          full.names=TRUE,pattern = ".csv")
  
  #Run analysis for each threshold from each hemisphere saved in 'grouppath'
  for(i in 1:length(grouppath)){
    #Pull thickness file
    groupfile <- basename(grouppath[i])
    grouptitle <- gsub(".csv", "", basename(grouppath[i]))
    thickness <- read.csv(grouppath[i], header = TRUE, stringsAsFactors = FALSE)
    
    #Linear Model
    #===========================================================================
    #Create empty data frame for appending results
    lmresults <- data.frame()
    
    #Run linear model for each function listed in the lmfxn variable
    #Save output and the associated residual and qqplots
    lmanalysis <- function(lmfxn, thickness){
      pdf(file=paste0(startpath,groupstr,"/linear_model/lmplot/",
                      grouptitle,".lmplot.pdf"))
      par(mfrow=c(4,2),mar=c(2,2,2,2)+0.1,ps=9)
      for(i in 1:length(lmfxn)){
        grouplm <- lm(lmfxn[[i]],thickness) #[[]] because lmfxn is a list
        grouplmsum <- data.frame(summary(grouplm)$coefficients)
        grouplmsum$call <- paste(deparse(lmfxn[[i]]))
        grouplmsum$file <- grouptitle
        lmresults <- rbind(lmresults,grouplmsum)
        lmplot1 <- plot(grouplm,which=1)
        title(paste(deparse(lmfxn[[i]])),cex.main=1,line=1.3)
        lmplot2 <- plot(grouplm, which=2)
      }
      dev.off()
      write.csv(lmresults, file=paste0(startpath,groupstr,"/linear_model/",
                                       grouptitle,".lm.csv"),row.names=TRUE)
    }
    
    #Run lmanalysis function
    lmanalysis(lmfxn, thickness)
    #===========================================================================
    
    #Logistic Regression
    #===========================================================================
    #Create empty data frame for appending results
    logresults <- data.frame()
    ORresults <- data.frame()
    
    #Run logisitic regression and odds ratio for each function listed in logfxn 
    #Save output
    loganalysis <- function(logfxn, thickness){
      for(i in 1:length(logfxn)){
        grouplog <- glm(logfxn[[i]],data=thickness,family=binomial(link='logit'))
        grouplogsum <- data.frame(summary(grouplog)$coefficients)
        grouplogsum$call <- paste(deparse(logfxn[[i]]))
        grouplogsum$file <- grouptitle
        logresults <- rbind(logresults, grouplogsum)
        
        grouplogOR <- data.frame(exp(cbind("OR"=coef(grouplog),confint(grouplog))))
        grouplog$call <- paste(deparse(logfxn[[i]]))
        grouplog$filename <- grouptitle
        ORresults <- rbind(ORresults, grouplogOR)
      }
      
      write.csv(logresults,file=paste0(startpath,groupstr,"/logistic_regression/",
                                       grouptitle,".log.csv"),row.names=TRUE)
      write.csv(ORresults,file=paste0(startpath,groupstr,"/logistic_regression/oddsratio/",
                                      grouptitle,".or.csv"),row.names=TRUE)
    }
    
    #Run loganalysis function
    loganalysis(logfxn, thickness)
    #===========================================================================
    
    #ROC Curves
    #===========================================================================
    rocanalysis <- function(throcfxn){
      #ROC model
      thROC <- roc(throcfxn,thickness,levels=c("negative","positive"),
                   auc=TRUE,smooth=FALSE)
      thROCparam <- data.frame(coords(thROC,"best",ret=c("threshold","specificity","sensitivity"),
                                      best.method="youden",transpose=FALSE))
      rownames(thROCparam)[1] <- paste0(grouptitle)
      thROCparam$auc <- thROC$auc
      write.csv(thROCparam,file=paste0(startpath,groupstr,"/roc/",grouptitle,".throc.csv"),
                row.names=TRUE)
      #ROC curve plot
      png(paste0(startpath,groupstr,"/roc/rocplot/",grouptitle,".throcplot.png"))
      thROCplot <- plot(thROC,print.thres="best",print.thres.best.method="youden",
                        print.auc=TRUE,main=paste0(grouptitle),cex.main=0.8)
      dev.off()
    }
    
    #Run rocanalysis function
    rocanalysis(throcfxn)
    #===========================================================================
    
    #Residual, ThicknessvAge, AmyloidvThickness Plots
    #===========================================================================
    statplots <- function(adjfxn1, adjfxn2, ageplotaes, amyloidplotaes,
                          centplotaes){
      #Linear model to plot adjusted thickness and adjusted amyloid (age+gender)
      adjlm1 <- lm(adjfxn1, thickness)
      adjlm2 <- lm(adjfxn2, thickness)
      
      #ggplot of LM models run above.
      adjplot <- ggplot()+
        geom_point(aes(adjlm1$residuals,adjlm2$residuals),size=2)+
        labs(color="Amyloid Status")+
        xlab(deparse(adjfxn1))+
        ylab(deparse(adjfxn2))+
        ggtitle("Residual Plot")+
        theme(plot.title=element_text(hjust=0.5))
      ggsave(paste0(startpath,groupstr,"/plots/adjustedplot/",grouptitle,".adjustedplot.png"),
             plot=adjplot,width=20,height=15,units="cm")
      
      #ThicknessvAge Plot
      ageplot <- ggplot(data=thickness,ageplotaes)+
        geom_point()+
        xlab("Age")+
        ylab("Thickness (mm)")+
        ggtitle(paste0(grouptitle))+
        geom_smooth(method="lm",se=FALSE)
      ggsave(paste0(startpath,groupstr,"/plots/thickness.age/",grouptitle,".ageplot.png"),
             plot=ageplot,width=20, height=15,units="cm")
      
      #DIAN PIB plot
      if(DIANcohort==1 & is.null(centplotaes)){
        pibplot <- ggplot(data=thickness, amyloidplotaes)+
          geom_point(size=2)+
          labs(color="Amyloid Status")+
          xlab("PIB - PET_fSUVR_rsf_TOT_CORTMEAN")+
          ylab("Thickness (mm)")+
          ggtitle(paste0(grouptitle))+
          geom_smooth(method="lm",se=FALSE)
        ggsave(paste0(startpath,groupstr,"/plots/thickness.pib/",grouptitle,".pibplot.png"),
               plot=pibplot,width=20,height=15,units="cm")
      }
      
      #ADRC Amyloid plots
      if(ADRCcohort==1 & !is.null(centplotaes)){
        pibonly <-thickness[thickness$Tracer == "PIB", ] 
        pibplot <- ggplot(data=pibonly, amyloidplotaes)+
          geom_point(size=2)+
          labs(color="Amyloid Status")+
          xlab("PIB - PET_fSUVR_rsf_TOT_CORTMEAN")+
          ylab("Thickness (mm)")+
          ggtitle(paste0(grouptitle))+
          geom_smooth(method="lm",se=FALSE)
        ggsave(paste0(startpath,groupstr,"/plots/thickness.pib/",grouptitle,".pibplot.png"),
               plot=pibplot,width=20,height=15,units="cm")
        
        av45only <-thickness[thickness$Tracer == "AV45", ] 
        av45plot <- ggplot(data=av45only, amyloidplotaes)+
          geom_point(size=2)+
          labs(color="Amyloid Status")+
          xlab("AV45 - PET_fSUVR_rsf_TOT_CORTMEAN")+
          ylab("Thickness (mm)")+
          ggtitle(paste0(grouptitle))+
          geom_smooth(method="lm",se=FALSE)
        ggsave(paste0(startpath,groupstr,"/plots/thickness.av45/",grouptitle,".av45plot.png"),
               plot=av45plot,width=20,height=15,units="cm")
        
        centplot <- ggplot(data=thickness, centplotaes)+
          geom_point(size=2)+
          labs(color="Amyloid Status")+
          xlab("Centiloid")+
          ylab("Thickness (mm)")+
          ggtitle(paste0(grouptitle))+
          geom_smooth(method="lm",se=FALSE)
        ggsave(paste0(startpath,groupstr,"/plots/thickness.centiloid/",grouptitle,".centplot.png"),
               plot=centplot,width=20,height=15,units="cm")
      }
    }
    
    #Run statplots function
    statplots(adjfxn1, adjfxn2, ageplotaes, amyloidplotaes,centplotaes)
    #===========================================================================
  }
  
  #ROC curve for total hippocampal volume (run once for each group)
  #=============================================================================
  #ROC model
  hcvROC <- roc(hcvrocfxn,thickness,levels=c("negative","positive"),
                auc=TRUE,smooth=FALSE)
  hcvROCparam <- data.frame(coords(hcvROC,"best",ret=c("threshold","specificity","sensitivity"),
                                   best.method="youden",transpose=FALSE))
  rownames(hcvROCparam)[1] <- paste0(groupstr,".tHCVVolume")
  hcvROCparam$auc <- hcvROC$auc
  write.csv(hcvROCparam,file=paste0(startpath,groupstr,"/roc/",groupstr,".hcvroc.csv"),
            row.names=TRUE)
  #ROC curve plot
  png(paste0(startpath,groupstr,"/roc/rocplot/",groupstr,".hcvrocplot.png"))
  hcvROCplot <- plot(hcvROC,print.thres="best",print.thres.best.method="youden",
                     print.auc=TRUE,main=paste0(groupstr,".tHCV"),
                     cex.main=0.8)
  dev.off()
  #=============================================================================
  
  #Concatenate data 
  #=============================================================================
  if(IndividROI==1){
    #Linear model
    mergelmpath <- paste0(startpath,groupstr,"/linear_model/")
    mergelmimport <- list.files(path=paste0(mergelmpath),full.names=TRUE,
                               pattern=".lm.csv")
    mergelm <- do.call("rbind",lapply(mergelmimport,function(i){
      read.csv(i,stringsAsFactors=FALSE)}))
    write.csv(mergelm,file=paste0(mergelmpath,groupstr,".mergedlm.csv"),
              row.names=FALSE)
    
    #Logisitic Regression
    mergelogpath <- paste0(startpath,groupstr,"/logistic_regression/")
    mergelogimport <- list.files(path=paste0(mergelogpath),full.names=TRUE,
                               pattern=".log.csv")
    mergelog <- do.call("rbind",lapply(mergelogimport,function(i){
      read.csv(i,stringsAsFactors=FALSE)}))
    write.csv(mergelog,file=paste0(mergelogpath,groupstr,".mergedlog.csv"),
              row.names=FALSE)
    
    #ROC
    mergerocpath <- paste0(startpath,groupstr,"/roc/")
    mergerocimport <- list.files(path=paste0(mergerocpath),full.names=TRUE,
                                pattern=".throc.csv")
    mergeroc <- do.call("rbind",lapply(mergerocimport,function(i){
      read.csv(i,stringsAsFactors=FALSE)}))
    write.csv(mergeroc,file=paste0(mergerocpath,groupstr,".mergedthroc.csv"),
              row.names=FALSE)
  }
  #=============================================================================
} #End of statanalysis function
#===============================================================================
#*******************************************************************************

#Run statanalysis function
#===============================================================================
if(DIANcohort==1){
  statanalysis(PreclinicalADADControls,startpath,ADRCcohort,DIANcohort,
               IndividROI,lmfxn,logfxn,throcfxn,hcvrocfxn,adjfxn1,adjfxn2,
               ageplotaes,amyloidplotaes,centplotaes=NULL)
  statanalysis(PreclinicalADADQ4Controls,startpath,ADRCcohort,DIANcohort,
               IndividROI,lmfxn,logfxn,throcfxn,hcvrocfxn,adjfxn1,adjfxn2,
               ageplotaes,amyloidplotaes,centplotaes=NULL)
}

if(ADRCcohort==1){
  statanalysis(PreclinicalADControls,startpath,ADRCcohort,DIANcohort,
               IndividROI,lmfxn,logfxn,throcfxn,hcvrocfxn,adjfxn1,adjfxn2,
               ageplotaes,amyloidplotaes,centplotaes)
  statanalysis(PreclinicalADQ4Controls,startpath,ADRCcohort,DIANcohort,
               IndividROI,lmfxn,logfxn,throcfxn,hcvrocfxn,adjfxn1,adjfxn2,
               ageplotaes,amyloidplotaes,centplotaes)
}
#===============================================================================