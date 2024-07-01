''# =====================================================================================================================
# Model Mangrove Vegetation Height from SAR VV and VH + LS reflectance and Indices
# =====================================================================================================================
# version:	v1 
# date: 	2022-03-15
# author: 	Daniel Gann, Boya Zhang
# 
# data used:
#	- gridded GliHT data (1m resolutin)
#	- SAR -- VH and VV time series 
#	- Lansat data: reflectance values + indices --list
# =====================================================================================================================
# PROCESSING SPECIFIC PARAMETERS AND DATA PATHS
# =================================================================================================
# libraries
#library(raster)
#library(rgdal)
library(tools)
library(parallel)
library(foreach)
library(sp)
library(caret)
library(sampling)
library(dplyr)
library(scales)
library(Metrics)
library(doParallel)
#library(sf)
#install.packages("pacman")
#library("pacman")
#pacman::p_load(pacman, dplyr, GGally, ggplot2, ggthemes, ggvis, httr, lubridate,plotly, rio, rmarkdown, shiny, stringr,tidyr)
# -------------------------------------------------------------------------------------------------
# set memory allocation for current session:  Free pysical memory - 1GB in MB
#memory.limit(size=round(as.numeric(gsub('\r','',gsub('FreePhysicalMemory=','',system('wmic OS get FreePhysicalMemory /Value',intern=TRUE)[3])))/1024 - 1000,0))
# set number of cores to use for parallel processing
#numCORE <- detectCores()
# -------------------------------------------------------------------------------------------------
# set raster processing options
#rasterOptions(tmpdir='')
# ---------------------------------------------------------------------------------------------

setwd('/beanstore/Personal/Paul/Third/')

#cl1 <- makeCluster(numCORE, outfile ='doParLog_simulation.txt')
#doParallel::registerDoParallel(cl1)
# =====================================================================================================================
meanst <- 'true'
prop <- 0.01 # This is the proportion for selecting training samples
nSAR <- 12
option <- 1 ## option 1 is sar, 2 is optical, 3 is both
iter <- 1 ## 50 iterations for differnt samplings
getwd()

#load('/beanstore/Personal/Paul/Third/Pre/smpdf.RData')
memory.limit(size=round(as.numeric(gsub('\r','',gsub('FreePhysicalMemory=','',system('wmic OS get FreePhysicalMemory /Value',intern=TRUE)[3])))/1024 - 1000,0))
# set number of cores to use for parallel processing
numCORE <- detectCores()
cl1 <- makeCluster(numCORE, outfile ='doParLog_simulation.txt')
doParallel::registerDoParallel(cl1)
SNWE <- c(2811613, 2812264, 485160, 491034)
print('option is')
print(option)
pathout <- '/beanstore/Personal/Paul/Third/Result/Temp/Mix/'
FID <- function(tm) {
  
  nSAR <- 12
  preProc <- 'no'
  setwd(paste('.',tm,sep='/'))
  getwd()
  
  # ----------------------------------------------------------------------------------------------------------------------
  #d <- read.csv('__dataraw_sw.csv') ## sw instead of no sw
  #d <- read.csv('__dataraw.csv') ## original data
  
  d <- read.csv('__data_allvar.csv')
  #d <- read.csv('__dataraw_sw.csv')
  colnames(d)[which(names(d) == 'ch_rfull')] <- 'LDR_HGT'
  colnames(d)[which(names(d) == 'class')] <- 'L2_Name'
  
  # process type: ldrMdl = model from pre-processed metrics files (config default); preProc = preprocess metric files for modeling; 
  procType <- 'mdl'    # model = mdl; predict = prd
  
  
  
  # ==========================================================================================================================
  # pre-process: join reference data and SAR metrics, and save as standardized data file
  
  varDpd <- 'LDR_HGT'
  mdlLst <- c('RF')
    
  table(d$VegCode);table(d$L2_Name);table(d$L4_Name)
  #table(d$L4_Name,d$L2_Name)	
  # ----------------------------------------------------------------------------------------------------------------------
  # set variable
  varLst00 <- 'L2_Name'
  varLst01 <- c('B2','B3','B4','B5','swir1', 'swir2', 'ndvi','ndmi','savi','wavi')
  varLst02 <- c(paste('vv',c(1:nSAR),sep=''))
  exclude <- c('vv1', 'vv8')
  varLst02 <- varLst02[!varLst02 %in% exclude]
  varLst03 <- c(paste('vh',c(1:nSAR),sep=''))
  exclude <- c('vh1', 'vh8')
  varLst03 <- varLst03[!varLst03 %in% exclude]
  varLst04 <-  c('vvmean','vvstd','vhmean','vhstd')
  varLst05 <- c('vv6', 'vh6') ## just single vv6 and vh6 value  
  varLst06 <- c('coh')
  ###### -----------------------------calculate sum and std for vv and vh---------
  copol <- d[varLst02]
  crosspol <- d[varLst03]
  sumv <- Reduce("+", copol) 
  sumh <- Reduce("+", crosspol) 
  meanv <- t(t(sumv) / 12)
  meanh <- t(t(sumh) / 12)
  d$vvmean <- meanv
  d$vhmean <- meanh
  sdv <- apply(copol, 1, sd) 
  sdh <- apply(crosspol, 1, sd)
  d$vvstd <- sdv
  d$vhstd <- sdh
  d <- d %>% mutate(across(is.numeric, round, digits=2))
  ################################# The end #########################
  if (meanst == 'false') {
    mdlVarLst0 <- c(varLst00,varLst01,varLst02)
    mdlVarLst1 <- c(varLst01,varLst02)
    mdlVarLst2 <- c(varLst00,varLst01,varLst03)
    mdlVarLst3 <- c(varLst01,varLst03)
    mdlVarLst4 <- c(varLst00,varLst01,varLst02,varLst03)
    mdlVarLst5 <- c(varLst01,varLst02,varLst03)
    varLstOLst <- list(mdlVarLst0,mdlVarLst1,mdlVarLst2,mdlVarLst3,mdlVarLst4,mdlVarLst5) # 
    varLstOLstNme <- c('vv_cls','vv','vh_cls','vh','vv_vh_cls','vv_vh') 
    
  } else {
    
    if (option == 1){
      mdlVarLst5 <- c(varLst00,varLst02,varLst03)
    }else if (option == 2){
      mdlVarLst5 <- c(varLst00,varLst01)
    } else {
      mdlVarLst5 <- c(varLst00,varLst01, varLst02,varLst03)
    }
    varLstOLst <- list(mdlVarLst5) # 
    varLstOLstNme <- c('s_vv_vh_coh') 
  }
  # ----------------------------------------------------------------------------------------------------------------------
  # subset predicton data, which is all of the data. Each obs includes 141 variables from all 6 combinations
  d.prd <- d[complete.cases(d[unlist(varLstOLst)]),]
  # ----------------------------------------------------------------------------------------------------------------------
  # subset training data # d.trn includes all the data with CH values.
  d.trn <- d[complete.cases(d), ]
  opt1 <- d.trn$B2
  diff(range(opt1))
  
  ### filtering data based on optical and SINGLE sar value #######
  d.filt <- subset(d.trn, d.trn$B2 > -1000 # the nan value is -9999
                    & d.trn$B3 > -1000
                    & d.trn$B4 > -1000
                    & d.trn$B5 > -1000
                    & d.trn$swir1 > -1000
                    & d.trn$swir2 > -1000
                    & d.trn$ndvi > -1000
                    & d.trn$ndmi > -1000
                    & d.trn$savi > -1000
                    & d.trn$wavi > -1000, ) ##& state %in% c('CA','NY'))
  d.filt2 <- subset(d.filt, d.filt$vv1 < 0 & d.filt$vh1 < 0 &
                            d.filt$vv2 < 0 & d.filt$vh2 < 0 &
                            d.filt$vv3 < 0 & d.filt$vh3 < 0 &
                            d.filt$vv4 < 0 & d.filt$vh4 < 0 &
                            d.filt$vv5 < 0 & d.filt$vh5 < 0 &
                            d.filt$vv6 < 0 & d.filt$vh6 < 0 &
                            d.filt$vv7 < 0 & d.filt$vh7 < 0 &
                            d.filt$vv8 < 0 & d.filt$vh8 < 0 &
                            d.filt$vv9 < 0 & d.filt$vh9 < 0 &
                            d.filt$vv10 < 0 & d.filt$vh10 < 0 &
                            d.filt$vv11 < 0 & d.filt$vh11 < 0 &
                            d.filt$vv12 < 0 & d.filt$vh12 < 0,  )
  
  #d.filt2 <- subset(d.filt, d.filt$vv6 < 0 & d.filt$vh6 < 0)
  
  ################ filter from the d.trn the rows with full features #######
  table(d.trn$L2_Name)
  # ----------------------------------------------------------------------------------------------------------------------
  
  # MODEL LiDAR ELEVATION
  
  # remove non-tree data (filter: height > 2 m and less than 35 m and NDVI != -9999)
  
  #hmin <- 1.99
  hmin <- 1.99
  hmax <- 25.01
  
  ##d.trn.tr <- d.trn[d.trn$LDR_HGT >= hmin & d.trn$LDR_HGT <= hmax, ] #& !d.trn$ndvi == -9999,] ## now the data is clean no need for ndvi filter
  d.filt2.tr <- d.filt2[d.filt2$LDR_HGT >= hmin & d.filt2$LDR_HGT <= hmax, ] #& !d.trn$ndvi == -9999,] ## now the data is clean no need for ndvi filter
  
  min(d.filt2.tr$ndvi)
  length(which(d.filt2.tr$ndvi < 0.4))
  #write.csv(d.trn.tr,"D:/Paul/re/CH/Mosaic/index.csv", row.names = FALSE)
  d.filt2.tr %>%
    group_by(d.filt2.tr$L2_Name) %>%
    summarise(mn=mean(LDR_HGT),sd=sd(LDR_HGT))
  
  table(d.filt2.tr$L2_Name)
  
  #return(d.filt2.tr)
  return(d.filt2.tr)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}

print('start')
tm <- 'Pre'   # 'Mosaic')
avail_pre = FID(tm)
print('Pre Finished')

setwd('/beanstore/Personal/Paul/Third/')
tm <- 'Post'   # 'Mosaic')
avail_post = FID(tm)
print('Post Finished')

selected_pre <- avail_pre[avail_pre$x > SNWE[3] & avail_pre$x < SNWE[4] & avail_pre$y > SNWE[1] & avail_pre$y < SNWE[2], ]
selected_post <- avail_post[avail_post$x > SNWE[3] & avail_post$x < SNWE[4] & avail_post$y > SNWE[1] & avail_post$y < SNWE[2], ]
print(nrow(avail_pre))
print(nrow(avail_post))
print('pre before the intersection')
print(nrow(selected_pre))
print('post before the intersection')
print(nrow(selected_post))

fid_pre = avail_pre$FID
fid_post = avail_post$FID
common = intersect(fid_pre, fid_post)

com_pre = avail_pre[avail_pre$FID %in% common, ] 
com_post = avail_post[avail_post$FID %in% common, ]

spre <- com_pre[com_pre$x > SNWE[3] & com_pre$x < SNWE[4] & com_pre$y > SNWE[1] & com_pre$y < SNWE[2], ]
spost <- com_post[com_post$x > SNWE[3] & com_post$x < SNWE[4] & com_post$y > SNWE[1] & com_post$y < SNWE[2], ]
print('pre after the intersection')
print(nrow(spre))
print('post after the intersection')
print(nrow(spost))

print('finished defining pre and post')
if(meanst == 'true'){
  vind = 1
}else {
  
  vind = 6
}

varLst00 <- 'L2_Name'
varLst01 <- c('B2','B3','B4','B5','swir1', 'swir2', 'ndvi','ndmi','savi','wavi')
varLst02 <- c(paste('vv',c(1:nSAR),sep=''))
exclude <- c('vv1', 'vv8')
varLst02 <- varLst02[!varLst02 %in% exclude]
varLst03 <- c(paste('vh',c(1:nSAR),sep=''))
exclude <- c('vh1', 'vh8')
varLst03 <- varLst03[!varLst03 %in% exclude]
varLst04 <-  c('vvmean','vvstd','vhmean','vhstd')
varLst05 <- c('vv6', 'vh6') ## just single vv6 and vh6 value  

if (meanst == 'false') {
  mdlVarLst0 <- c(varLst00,varLst01,varLst02)
  mdlVarLst1 <- c(varLst01,varLst02)
  mdlVarLst2 <- c(varLst00,varLst01,varLst03)
  mdlVarLst3 <- c(varLst01,varLst03)
  mdlVarLst4 <- c(varLst00,varLst01,varLst02,varLst03)
  mdlVarLst5 <- c(varLst01,varLst02,varLst03)
  varLstOLst <- list(mdlVarLst0,mdlVarLst1,mdlVarLst2,mdlVarLst3,mdlVarLst4,mdlVarLst5) # 
  varLstOLstNme <- c('vv_cls','vv','vh_cls','vh','vv_vh_cls','vv_vh') 
  
} else {
  
  if (option == 1){
		mdlVarLst5 <- c(varLst00,varLst02,varLst03)
	}else if (option == 2){
		mdlVarLst5 <- c(varLst00,varLst01)
	} else {
		mdlVarLst5 <- c(varLst00,varLst01, varLst02,varLst03)
    }
  
  varLstOLst <- list(mdlVarLst5) # 
  varLstOLstNme <- c('s_vv_vh') 
}

print('finished determining variables to use')
######### Random Forest training and testing #######
# Use Pre for training Post for testing
# number of resample iterations

# random samples per class
smpNum <- prop
strVar <- 'L2_Name'
#strVar <- c(3,1,4,2)
a <- table(com_pre$L2_Name) ## Pre for training
freqclass <- as.numeric(a)
trfreq <- round(freqclass * prop, digits = 0)
#sizeArg <- rep(smpNum,length(unique(d[strVar])[,1]))
sizeArg <- trfreq
com_pre=com_pre[order(com_pre$L2_Name),] ## reordered d.trn.t
com_post = com_post[order(com_post$L2_Name),]

cntr <- trainControl(method='repeatedcv', 
                     number=3, 
                     repeats=3)


# collector data frame
mdlRes <- data.frame(mdlNme=character(),
                     iter=numeric(), 
                     smp=numeric(),
                     rmse_cv=numeric(),
                     mae_cv=numeric(),
                     rae_cv=numeric(),
                     r2_cv = numeric(),
                     rmse_val=numeric(),
                     mae_val=numeric(),
                     rae_val=numeric(),
                     r2_val = numeric(),
                     varLstNme=character(),
                     stringsAsFactors=FALSE)
i <- 1
# # run iterations
print("finished before the smp set up")
print(dim(sizeArg))
print(dim(com_pre))
print(dim(com_post))
smp <- sampling:::strata(com_pre, strVar, sizeArg, method=c('srswor'))# used re_ordered here.
print("before FID")
smp <- smp$FID

load("/beanstore/Personal/Paul/Third/smpdf_all_filter.Rdata") 

iter <- 10
for (i in 1:iter){
  
	print(paste('iteration ',i,' of ',iter))
	# model name
	smpIter <- paste('smp',smpNum,'_i',i,sep='')
	if (i != 2) {
		# Continue to the next iteration when i is 3
		next
	}
	# stratified random sample for calibration
	#smp <- sampling:::strata(com_pre, strVar, sizeArg, method=c('srswor'))# used re_ordered here.
	#save(smp, file = "smp.Rdata")
	#load('/beanstore/Personal/Paul/Third/smp.Rdata')
	#smp <- smpdf[,i] ## use different sampling for training and testing for each time
	
	#calPnt <- sampling:::getdata(com_pre, smp)	
	#calPnt <- com_pre[smp,]
	#valPnt <- com_post[smp,]
	
	smp <- smpdf[,i]
	calPnt <- com_pre[smp,]
	valPnt <- com_post[smp,]
	print("the size matches")
	# remove calibration data from data set to create validation set
	#valPnt <- d.filt2.tr[!d.filt2.tr$FID %in% calPnt$FID ,]
	#valPnt <- sampling:::getdata(com_post, smp)


	varLst <- vind
	mdlVarLst <- varLstOLst[[varLst]]
	mdlVarNme <- varLstOLstNme[varLst]
	print(mdlVarNme)

	# subset calibration and validation data
	calDat <- calPnt[mdlVarLst]
	valDat <- valPnt[mdlVarLst]
	  # ----------------------------------------------------------------------------------------------------------------------	
	  
	  # RF
	print("before training")
	print(nrow(calPnt))
	print(ncol(calPnt))
	print(nrow(calPnt$LDR_HGT))
	print(nrow(com_pre))
	print(nrow(com_post))
	
	(start_time <- Sys.time())
	summary(RF <- caret::train(calDat,calPnt$LDR_HGT,method='rf', importance=T, trControl = cntr, ntree=500))
	(end_time <- Sys.time())
	write.csv(round(randomForest::importance(RF$finalModel, scale=FALSE), 2),paste(pathout, '_mdl_RF',as.character(option),smpIter,'_',mdlVarNme,'__varImp.csv',sep=''))

	# save the models to disk
	#saveRDS(RF, paste('_mdl_RF',as.character(option),smpIter,'_',mdlVarNme,'.rds',sep=''))

	# predict calibration data
	#cal.prdRF <- predict(RF)

	# apply to validation set and calculate rmse and rae
	val.prdRF <- predict(RF, newdata = valDat)

	# write results to data frame

	rsq <- function (x, y) cor(x, y) ^ 2
	#r2_cal <- rsq(calPnt$LDR_HGT, cal.prdRF)
	r2_val <- rsq(valPnt$LDR_HGT, val.prdRF)

	mdlRes <- rbind(mdlRes,c('RF',i,smpNum,
							 #round(rmse(calPnt$LDR_HGT,cal.prdRF),3),round(mae(calPnt$LDR_HGT,cal.prdRF),3),round(rae(calPnt$LDR_HGT,cal.prdRF),3), r2_cal,
							 round(rmse(valPnt$LDR_HGT,val.prdRF),3),round(mae(valPnt$LDR_HGT,val.prdRF),3),round(rae(valPnt$LDR_HGT,val.prdRF),3), r2_val,
							 mdlVarNme))
	#print(rmse(calPnt$LDR_HGT,cal.prdRF),3)
	print(rmse(valPnt$LDR_HGT,val.prdRF),3)
	#print(r2_cal)
	print(r2_val)
	print("This hits the end of current training/testing")
	#output <- data.frame(smp, com_pre[smp, 2],com_pre[smp, 3], val.prdRF, valPnt$LDR_HGT, val.prdRF - valPnt$LDR_HGT)
	output <- data.frame(seq(1, nrow(valPnt)), valPnt[, 1], valPnt[, 2],valPnt[, 3], val.prdRF, valPnt$LDR_HGT, val.prdRF - valPnt$LDR_HGT)
	colnames(output) <- c("ID", "FID", "x", "y", "pred", "true", "error")
	write.table(output, paste0(pathout, "predictions", as.character(option), smpIter,".csv"), row.names=FALSE, col.names=TRUE, sep=',')
	write.table(mdlRes, paste0(pathout, "uncertainty", as.character(option), smpIter, ".csv"), row.names=FALSE, col.names=TRUE, sep=',')
}

write.csv(mdlRes, paste0("Summary_uncertainty_", as.character(option), ".csv"), row.names=FALSE)
stopCluster(cl1);rm()
#names(mdlRes) <- c('mdlNme','iter','smp','rmse_cv','mae_cv','rae_cv','r2_cal','rmse_val','mae_val','rae_val','r2_val','varLstNme')

# write models results to file
#write.csv(mdlRes,paste('mdlComp_',tm,'_',smpIter,'.csv',sep=''))

# plt <- data.frame(valPnt$LDR_HGT,val.prdRF)

# b <- ggplot(plt, aes(x=valPnt$LDR_HGT, y=val.prdRF))
# b + geom_point(size = 0.1)

# ggplot(plot, aes(x=valPnt$LDR_HGT, y=valPnt$LDR_HGT)) +
  # geom_point(size=2, shape=23)
# err = valPnt$LDR_HGT - val.prdRF
# hist(err)
# loc <- which(err > 10)

# dif <- calDat - valDat
# b2dif <- calDat$B2 - valDat$B2
# chdif <- calPnt$LDR_HGT - valPnt$LDR_HGT

# difcheck <- data.frame(valPnt$LDR_HGT,val.prdRF)


    # apply to full data set
    #prdRF <- predict(RF, newdata = d.prd[mdlVarLst])
    
    # write to file
    #write.csv(prdRF,paste('_mdl_RF_',smpIter,'_',mdlVarNme,'_pred.csv',sep=''))			}
  

#save(smp, file = "smp.Rdata")
# add names to model result data frame
#names(mdlRes) <- c('mdlNme','iter','smp','rmse_cv','mae_cv','rae_cv','rmse_val','mae_val','rae_val','varLstNme')

############## output ########################

#write.csv(com_pre, "./common_pre.csv", row.names=FALSE)
#write.csv(com_post, "./common_post.csv", row.names=FALSE)
# write models results to file
#write.csv(mdlRes,paste('mdlComp_',tm,'_',smpIter,'.csv',sep=''))

# if (length(mdlLst) > 1){
#   # summarize stats
#   mdlResSum <- mdlRes %>% 
#     group_by(mdlNme) %>%
#     summarise(rmse_cv_mean=mean(as.numeric(rmse_cv)),rmse_cv_sd=sd(as.numeric(rmse_cv)),
#               rae_cv_mean=mean(as.numeric(rae_cv)),rae_cv_sd=sd(as.numeric(rae_cv)),
#               rmse_val_mean=mean(as.numeric(rmse_val)),rmse_val_sd=sd(as.numeric(rmse_val)),
#               rae_val_mean=mean(as.numeric(rae_val)),rae_val_sd=sd(as.numeric(rae_val)))
#   
#   # write model summary to file
#   write.csv(mdlResSum,paste('_mdlCompSum_',tm,'_',smpIter,'.csv',sep=''))			
# }
# model <- lm(valPnt$LDR_HGT ~ val.prdRF)
# summary(model)

#model <- lm(valPnt_t$LDR_HGT ~ val_t.prdRF)
#summary(model)

#prdRF <- predict(RF, newdata = d_t.prd[mdlVarLst])

# write =to file
# write.csv(mdlRes_t,paste('mdlComp_cross_',tm,'_',smpIter,'.csv',sep=''))
# op <- data.frame(valPnt_t$LDR_HGT, val_t.prdRF, valPnt_t$L2_Name)
# opself <- data.frame(valPnt$LDR_HGT, val.prdRF, valPnt$L2_Name)
# 
# write.csv(op,paste('mdlComp_cross_compare.csv',sep=''))
# write.csv(opself,paste('mdlComp_compare.csv',sep=''))
# mdlResSum <- mdlRes %>% 
#  group_by(mdlNme) %>%
#  summarise(rmse_val_mean=mean(as.numeric(rmse_val)),rmse_val_sd=sd(as.numeric(rmse_val)),
#            rae_val_mean=mean(as.numeric(rae_val)),rae_val_sd=sd(as.numeric(rae_val)))

# write model summary to file
#write.csv(mdlResSum,paste('_mdlCompSum_cross_',tm,'_',smpIter,'.csv',sep=''))	

# stop and remove parallel processing cluster
#stopCluster(cl1);rm()
# --------------------------------------------------------------------------------------------------------------------------
#EOF
