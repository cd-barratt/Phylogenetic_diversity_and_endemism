# This script is modified from Dan Rosauer's to calculate best combinations of predictor variables to explain PE (response variable)

#################################################################################################
# PE 35 (= all 35 intraspecific lineages)
################################f#################################################################
rm(list=ls())
library(SDMTools)
library(raster)
#library(relaimpo)
library(glmulti)

setwd("/Users/Chris/Desktop/glmulti_SARs/")

# first a function - main script follows below
glm_process = function( result_frame,
                        model_number,
                        predictor_text,
                        response_text,        
                        response_text_logistic = "div", # the response variable name, not values                        the_data,
                        do_spatial_LM = TRUE,
                        weight_list = NULL,
                        do_logistic = TRUE,
                        the_data
)
{
  cat("\n\n** ",model_number,result_frame[model_number,1],result_frame[model_number,2],result_frame[model_number,3]," **\n")
  
  the.formula <- paste(response_text,"~",predictor_text)
  glm_gauss <- glm(the.formula, data=the_data, family="gaussian")
  print(summary(glm_gauss))
  cat("\nGaussian model\naic:",glm_gauss$aic)
  dev_exp <- (glm_gauss$null.deviance-glm_gauss$deviance)/glm_gauss$null.deviance
  cat("\nDeviance explained:",dev_exp,"\n")
  result_frame[model_number,"glm_r2"] <- round(dev_exp,4)
  result_frame[model_number,"glm_aic"] <- round(glm_gauss$aic,1)
  
  if (do_logistic) {
    the.formula <- paste(response_text_logistic,"~",predictor_text)    
    glm_logistic <- glm(the.formula, data=the_data, family = "binomial")
    print(summary(glm_logistic))
    cat("\nLogistic model \naic:",glm_logistic$aic)
    dev_exp <- (glm_logistic$null.deviance-glm_logistic$deviance)/glm_logistic$null.deviance
    cat("\nDeviance explained:",dev_exp,"\n")
    result_frame[model_number,"glm_logistic_r2"] <- round(dev_exp,4)
    result_frame[model_number,"glm_logistic_aic"] <- round(glm_logistic$aic,1)
  }
  
  if (do_spatial_LM) {
    
    SAR_gauss <- errorsarlm(glm_gauss, data=the_data, listw=weight_list, quiet=FALSE, na.omit, zero.policy=TRUE, tol.solve=1e-11)    
    #calculate AIC
    SAR_gauss.aic <- (-2*SAR_gauss$LL)+(2*SAR_gauss$parameters)
    
    print(summary(SAR_gauss))
    cat("\nSAR gaussian \naic:",SAR_gauss.aic)
    result_frame[model_number,"sarlm_aic"] <- round(SAR_gauss.aic,1)
    
    if (do_logistic) {
      SAR_logistic <- errorsarlm(glm_logistic, data=the_data, listw=weight_list, quiet=FALSE, na.omit, zero.policy=TRUE, tol.solve=1e-11)
      #calculate AIC
      SAR_logistic.aic <- (-2*SAR_logistic$LL)+(2*SAR_logistic$parameters)
      
      print(summary(SARer_with_exp))
      cat("\nSARLM logistic \naic:",SAR_logistic.aic)
      result_frame[model_number,10] <- round(SAR_logistic.aic,1)
    }
    
  }
  
  return(result_frame)
}



base_path     <- "/Users/Chris/Desktop/glmulti_SARs/"
results.dir   <- paste(base_path,sep='')

regions <- data.frame()
do_logistic <- TRUE
logistic_threshold <- 0.95  # for a logistic model to predict membership of top diversity

#whether to do spatial autocorrelation, and at what radius
do_spatial_LM <- FALSE
spatial_LM_radius <- 3000

# whether to do delta AIC table and glmulti
do_glmulti <- TRUE
do_plots   <- FALSE

i <- 1
regions[i,"region"]        <- ''
regions[i,"current_forest"] <- paste(regions[i],"/current_forest_clipped.asc",sep='')

diversities <- data.frame(taxon="PE_35",
                          level="lineage", 
                          metric="endemism", 
                          grid="/phylogenetic_endemism35_clipped.asc",
                          transform="log",
                          stringsAsFactors = F)
j <- 2
diversities[j,1:4] <- c("PE_35","lineage","endemism","/phylogenetic_endemism35_clipped.asc")

#output <- data.frame()
k  <- 0

# create a data frame to hold the model results
result_frame <- result_frame <- data.frame(region="",response_description="",predictor_description="",n=0,glm_r2=0,glm_aic=0,glm_delta_aic=0,sarlm_aic=0,sarlm_delta_aic=0,glm_logistic_r2=0,glm_logistic_aic=0,glm_logistic_delta_aic=0,sarlm_logistic_aic=0,sarlm_logistic_delta_aic=0,stringsAsFactors = F)

if (do_spatial_LM) {
  library(spdep)
  
  
}

if (do_glmulti) {
  # create a result data frame
  glmulti_blank_results <- data.frame(region="",max_pred=0,AIC=-9999,deltaAIC=-9999,r2=0,formula="",stringsAsFactors = F)
 
  # add predictor columns
    glmulti_blank_results <- cbind(glmulti_blank_results,data.frame(current_forest=0,bio1=0,bio4=0,bio12=0,bio14=0,foreststability_120k=0,topographic_heterogeneity=0,anom_bio1=0,anom_bio12=0,stringsAsFactors = 0))
  glmulti_full_results <- glmulti_blank_results[-1,]
}

#Loop through each region, and within it, each diversity metric
for (i in 1:nrow(regions)) {
  
  # for the all region model, assign region sizes to cells
  
  #define some basic data
  rf.ras <- raster("current_forest_clipped.asc")# read in the vegetation grid
  rf.ras[which(is.finite(rf.ras[]) & rf.ras[] !=1)] <- 0  #set all veg != 1 (rainforests) to 0
  #rf_region.asc     <- read.asc(regions$region[i])
  current_forest.ras <- raster("current_forest_clipped.asc")
    
  #crop the rainforest raster to match the stability raster
  rf.ras = crop(x=rf.ras,y=current_forest.ras)
  rf.asc = asc.from.raster(rf.ras)
  current_forest.asc = asc.from.raster(current_forest.ras)
  
  j=1 # for now only do lineage endemism
  #for (j in 1:nrow(diversities)) {
  
  #load the diversity result at 0.01 degree resolution and resample to match stability
  div.ras     <-  raster(paste(results.dir,diversities$grid[j],sep=''))
  div_resample.ras <- resample(div.ras,current_forest.ras,method="bilinear")
  div_resample.asc <- asc.from.raster(div_resample.ras)
  
  #create pos_rf the table of grid cell values for rainforest
  pos_rf <- as.data.frame(which(is.finite(rf.asc),arr.ind=TRUE)) #get all points that have data
  
  xy_rf <- getXYcoords(rf.asc)
  pos_rf$long  <- xy_rf$x[pos_rf$row]
  pos_rf$lat  <- xy_rf$y[pos_rf$col]
  rm(xy_rf)
  
  pos_rf$rf = rf.asc[cbind(pos_rf$row,pos_rf$col)]            #append the current forest data
  #pos_rf = pos_rf[which(pos_rf$rf==1),]  # filter to rainforest areas
  pos_rf$div = div_resample.asc[cbind(pos_rf$row,pos_rf$col)] #append the PE data
  pos_rf = pos_rf[which(pos_rf$div > 0),]  # filter to areas with a PE score
  
  # add the predictors based on rainforest niche models
  pos_rf$current_forest = current_forest.asc[cbind(pos_rf$row,pos_rf$col)] #append the current forest scores
  
    # add other predictors
  predictors <- data.frame(name="bio1",description="bio1",path="bio1_clipped.asc",resample=FALSE,stringsAsFactors = F)
  predictors[2,] <- c("bio4","bio4","bio4_clipped.asc",FALSE)
  predictors[3,] <- c("bio12","bio12","bio12_clipped.asc",FALSE)
  predictors[4,] <- c("bio14","bio14","bio14_clipped.asc",FALSE)
  predictors[5,] <- c("topographic_heterogeneity","topographic_heterogeneity","topographic_heterogeneity_clipped.asc",FALSE)
  predictors[6,] <- c("foreststability_120k","foreststability_120k","120k_foreststability_clipped.asc",FALSE)
  predictors[7,] <- c("anom_bio1","anom_bio1","anom_bio1_clipped.asc",FALSE)
  predictors[8,] <- c("anom_bio12","anom_bio12","anom_bio12_clipped.asc",FALSE)

  
  # patch size
  
  for (p in 1:nrow(predictors) ) {
    path <- predictors$path[p]
    env.ras <- raster(path)
    if (predictors$resample[p]) {
      env.ras <- resample(env.ras,current_forest.ras,method="bilinear")
      cat(predictors[p,"name"], "done\n")
    }
    
    pos_rf[,predictors$name[p]] <- extract(env.ras,pos_rf[,c("long","lat")])
  }
  
  # transform the response variable if needed (as set in the diversities data frame)
  if (diversities[j,"transform"] == "log") {
    pos_rf$div <- log(pos_rf$div)
  }
  
  # rescale the predictors - KEEP THE COLUMN NUMBERS CORRECT
  #pos_rf[,7:15] <- scale(pos_rf[,7:15])
  pos_rf[,7:15] <- scale(pos_rf[,7:15])
  
    # add a binary variable for membership of top class
  thresh <- quantile(pos_rf$div,logistic_threshold)
  pos_rf$div_binary <- rep(0,nrow(pos_rf))
  pos_rf$div_binary[pos_rf$div >= thresh] <- 1
  
  # calculate the distance matrix for SARLM    
  if (do_spatial_LM) {
    cat("\n\nPreparing neighbourhood weights for SARLM\n\n")
    coords<-as.matrix(cbind(pos_rf$row, pos_rf$col))
    cont.nb <- dnearneigh(coords,0,spatial_LM_radius,longlat=TRUE)
    weight_list <- nb2listw(cont.nb, glist=NULL, style="W", zero.policy=TRUE)
  }
  
  n <- nrow(pos_rf)
  
  if (do_glmulti) {
    
    # automated predictor selection
    
        formula <- "div~foreststability_120k+current_forest+topographic_heterogeneity+bio1+bio4+bio12+bio14+anom_bio1+anom_bio12" 
        
  }
  fullGLM      <- glm(data=pos_rf, formula=formula)
  region_startrow <- nrow(glmulti_full_results) + 1  # use this to calculate deltaAIC for just the region
  
  for (model_size in 1:9) {
    rm(glmulti_res)
    
    model_name <- paste(diversities[j,1],diversities[j,2],diversities[j,3],"for",regions[i,"region"])
    
    try(glmulti_res <- glmulti(y=fullGLM,level=1, maxsize=model_size, crit="aicc", method="h",name=model_name, report=T, confsetsize=512, plotty=F))
    if (exists("glmulti_res")) {
      
      glmulti_sum <- summary(glmulti_res)
      
      #####
      # fit the best model
      bestGLM <- glm(data=pos_rf, formula=glmulti_sum$bestmodel)
      bestGLM_r2 <- (bestGLM$null.deviance - bestGLM$deviance) / bestGLM$null.deviance
      bestGLM_r2 <- round(bestGLM_r2,3)
      bestGLM_sum <- summary(bestGLM)
      
        
      # store results in the data frame
      glmulti_new_results <- glmulti_blank_results
      glmulti_new_results$region[1]      <- regions$region[i]
      glmulti_new_results$max_pred[1]    <- model_size
      glmulti_new_results$AIC[1]         <- glmulti_sum$bestic
      glmulti_new_results$formula[1]     <- glmulti_sum$bestmodel
      glmulti_new_results$r2[1]          <- bestGLM_r2
      
      # get model terms and co-efficents
      terms <- attributes(bestGLM_sum$terms)$term.labels
      coef  <- bestGLM_sum$coefficients
      for (pred in terms) {
        glmulti_new_results[1,pred] <- round(coef[pred,"Estimate"],4)
      }
      
      glmulti_full_results <- rbind(glmulti_full_results,glmulti_new_results)
      
      cat("\n**********************************\nBest",model_size,"predictor model of",diversities[j,1],diversities[j,2],diversities[j,3],"for",regions[i,"region"])
      cat("\nBest model formula:",glmulti_sum$bestmodel,
          "\nBest model AICC:   ",glmulti_sum$bestic,
          "\nBest model r^2:    ",bestGLM_r2,"\n**********************************\n")
      
      print(summary(bestGLM))
    } else {
      glmulti_new_results             <- glmulti_blank_results
      glmulti_new_results$formula[1]  <- "FAILED"
      
      cat("\n**********************************\nBest",model_size,"predictor model of",diversities[j,1],diversities[j,2],diversities[j,3],"for",regions[i,"region"],
          "\n **  FAILED  **",
          "\n**********************************\n")
      
    }
  }
  
  # calculate deltaAIC
  rows    <- region_startrow:nrow(glmulti_full_results)
  minAIC  <- min(glmulti_full_results$AIC,na.rm=T)
  glmulti_full_results$deltaAIC[rows] <- glmulti_full_results$AIC[rows] - minAIC
}
}

if (do_glmulti) {
  setwd(results.dir)
  write.csv(glmulti_full_results,"Result_table_PE_35.csv")
  rm(rf.ras, rf.asc, div.ras, div_resample.ras, current_forest.asc, current_forest.ras,glmulti_blank_results,glmulti_new_results, rows,minAIC)
}

if (do_plots) {
  # now plot current suitability v stability (past suitability), coloured by endemism
  library(maptools)
  library(classInt)
  #windows()
  class_count <- 12
  my.class <- classIntervals(pos_rf$div,n=class_count,style="quantile", digits=2)
  my.class_breaks <- round(my.class[[2]],4)
  my.pal <- c("darkblue","green2","yellow","red")
  my.col <-findColours(my.class,my.pal)
  legend_cols <- attr(my.col,"palette")
  plot(pos_rf$foreststability_120k,pos_rf$current_forest,xlab="Current suitability",ylab="Stability",col=my.col)
  legend(x="topleft",legend=my.class_breaks[1:class_count+1],fill=legend_cols)
}

#write the weights

tmp <- weightable (glmulti_res)
tmp
write.csv(tmp,"Weight_table_PE_35.csv")
coef.glmulti(glmulti_res) # FROM DAN

############################################################################
#####################    SPATIAL    ########################################
############################################################################
############################################################################

#Load packages
library(spdep)
library(akima)
library(raster)


#----------------------------- PE DATASET -----------------------------------#

#########################
#Examine the structure of the PE dataset
PE <- pos_rf
names(PE)
str(PE)
#PE[complete.cases(PE[,7:11]),]

#Make maps
par(mfrow=c(4,3))
image(interp(PE$long, PE$lat, PE$div), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="PE")
image(interp(PE$long, PE$lat, PE$current_forest), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="current_forest")
image(interp(PE$long, PE$lat, PE$foreststability_120k), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="foreststability_120k")
image(interp(PE$long, PE$lat, PE$topographic_heterogeneity), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="topographic_heterogeneity")
image(interp(PE$long, PE$lat, PE$bio1), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="annual mean temp")
image(interp(PE$long, PE$lat, PE$bio12), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="annnual precipitation")
image(interp(PE$long, PE$lat, PE$bio14), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="precipitation of driest quarter")
image(interp(PE$long, PE$lat, PE$bio4), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="tenperature seasonality")
image(interp(PE$long, PE$lat, PE$anom_bio1), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="Quaternary temp change")
image(interp(PE$long, PE$lat, PE$anom_bio12), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="Quaternary precipitation change")

#########################
#Make a regression model to explain PE
hist(PE$div)
hist(log(PE$div))
formula <- "div~foreststability_120k+topographic_heterogeneity+bio1+bio4+bio12+bio14+anom_bio1+anom_bio12" 
lm_PE <- glm(data=pos_rf, formula=formula)
summary(lm_PE)
hist(residuals(lm_PE))

#########################
#Spatial structure of residuals
library(ncf)

par(mfrow=c(2,2))
#Correlograms with latlon = FALSE
cor.OBL<-correlog(PE$long, PE$lat, z=PE$div, na.rm=T, increment=0.5, resamp=1, latlon = FALSE)    #makes 10 distance classes
cor.OBL
plot(cor.OBL$correlation, type="b", pch=16, cex=1.2, lwd=1.5, ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)
title(main="0.5 degree")

#With latlon = TRUE
cor.OBL_latlon<-correlog(PE$long, PE$lat, z=PE$div, na.rm=T, increment=0.5, resamp=1, latlon = TRUE)  #uses km because latlon = TRUE
plot(cor.OBL_latlon$correlation, type="b", pch=16, cex=1.2, lwd=1.5, ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)    
title(main="0.5km")

#With latlon = TRUE and increment=10
cor.OBL_10<-correlog(PE$long, PE$lat, z=PE$div, na.rm=T, increment=10, resamp=999, latlon = TRUE)  #uses km because latlon = TRUE
plot(cor.OBL_10$correlation, type="b", pch=16, cex=1.2, lwd=1.5, ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)
title(main="10km")

summary(cor.OBL_10)
cor.OBL_10$mean.of.class
cor.OBL_10$n
cor.OBL_10$correlation      
cor.OBL_10$p                

#Correlogram for residuals
cor.res_10<-correlog(PE$long, PE$lat, z=residuals(lm_PE), na.rm=T, increment=10, resamp=999, latlon = TRUE)  #uses km because latlon = TRUE

#Plot both residuals and raw data
plot(cor.OBL_10$correlation, type="b", pch=16, col="black", cex=1.2, lwd=1.5, xlim=c(0,20), ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)                  
points(cor.res_10$correlation, pch=1, cex=1.2, col="blue")
lines(cor.res_10$correlation, lwd=1.5, col="blue")
title(main="10km raw [black] vs residuals [blue]")

#########################
#Implementing a spatial model

#Make coordinate list
coords_PE<-as.matrix(cbind(PE$long,PE$lat))
plot(coords_PE)

#Minimum distance to connect to at least one neighbor
PE_knear <- knn2nb(knearneigh(coords_PE, k=1))
summary(PE_knear)
dsts_PE<-unlist(nbdists(PE_knear, coords_PE, longlat = TRUE))
summary(dsts_PE)
max(dsts_PE)

#Neighbour defined by dnearneigh()
par(mfrow=c(1,1))
PE_nb<-dnearneigh(coords_PE,0,max(dsts_PE), longlat=T)
par(mfrow=c(1,1))
plot(PE_nb, coords_PE, pch=20, lwd=2)
summary(PE_nb)

PE_nb_10<-dnearneigh(coords_PE,0,10, longlat=T) #max(dsts_PE)=12.99333 #use 8? (Moran's I = 0 there)
summary(PE_nb_10)
plot(PE_nb_10, coords_PE, pch=20, lwd=0.5)

#Defining the spatial weights matrix
nb1_w<-nb2listw(PE_nb, glist=NULL, style="W", zero.policy=TRUE)
summary(nb1_w)
plot(nb1_w, coords_PE)

#Spatial autoregressive error model
formula <- "div~foreststability_120k+bio1+bio4+bio12+bio14+topographic_heterogeneity+anom_bio1+anom_bio12" 
the_data <- PE
sem_error_nb1_w<-errorsarlm(formula, data=the_data, listw=nb1_w, zero.policy=TRUE)
summary(sem_error_nb1_w)
errorsalm_summary <- summary(sem_error_nb1_w)
capture.output(errorsalm_summary, file ="SAR_summary_35.txt")

#########################
#Testing for spatial structure in the residuals of the spatial model

#Correlogram for spatial model
cor.SEM_10<-correlog(PE$long, PE$lat, z=residuals(sem_error_nb1_w), increment=10, resamp=1, na.rm=T, latlon = TRUE)  #uses km because latlon = TRUE

#Plot both residuals and raw data
plot(cor.OBL_10$correlation, type="b", pch=16, cex=1.2, lwd=1.5, xlim=c(0,15), ylim=c(-0.5, 1), xlab="Distance class (x 10 km)", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)                  
points(cor.res_10$correlation, pch=1, cex=1.2, col="blue")
lines(cor.res_10$correlation, lwd=1.5, col="blue")
lines(cor.SEM_10$correlation, lwd=2.5, col="red")
title(main="10km raw [black] vs residuals [blue] vs SEM model [red]")

#Calculate Moran's I value individually
#PE_nb_moran10<-dnearneigh(coords_PE,0,10, longlat=T)
#nb_moran10<-nb2listw(PE_nb_moran10, glist=NULL, style="W", zero.policy=TRUE)

#moran.test(PE$div, listw=nb_moran10, randomisation=TRUE, zero.policy=TRUE, alternative="two.sided")
#moran.test(residuals(lm_PE2), listw=nb_moran10, randomisation=TRUE, zero.policy=TRUE, alternative="two.sided")
#moran.test(residuals(sem_error_nb1_w), listw=nb_moran10, randomisation=TRUE, zero.policy=NULL, alternative="two.sided")

#moran.test(PE$div, listw=nb_moran10, randomisation=TRUE, zero.policy=TRUE, alternative="greater")
#moran.test(residuals(lm_frugi), listw=nb_moran1000, randomisation=TRUE, zero.policy=NULL, alternative="greater")
#moran.test(residuals(sem_error_nb1_w), listw=nb_moran1000, randomisation=TRUE, zero.policy=NULL, alternative="greater")                                   

#Residual maps
#par(mfrow=c(2,2))
#image(interp(birds$X, birds$Y, birds$Frugivores), col = terrain.colors(12), main="Frugivore richness")
#image(interp(birds$X, birds$Y, residuals(lm_frugi)), col = terrain.colors(12), main="OLS model residuals")
#image(interp(birds$X, birds$Y, residuals(sem_error_nb1_w)), col = terrain.colors(12), main="SEM residuals")

#plot(coords_birds, col=c("blue", "red")[sign(residuals(lm_frugi))/2+1.5], pch=19,
#cex=abs(residuals(lm_frugi))/max(residuals(lm_frugi)), xlab="geographical x-coordinates", ylab="geographical y-coordinates")

#plot(coords_birds, col=c("blue", "red")[sign(residuals(sem_error_nb1_w))/2+1.5], pch=19,
#    cex=abs(residuals(sem_error_nb1_w))/max(residuals(sem_error_nb1_w)), xlab="geographical x-coordinates", ylab="geographical y-coordinates")

#########################
#Improving the spatial model

#Second spatial weights matrix
#birds_nb_200<-dnearneigh(coords_birds,0,200, longlat=T)
#summary(birds_nb_200)
#nb2_w<-nb2listw(birds_nb_200, glist=NULL, style="W", zero.policy=FALSE)
#summary(nb2_w)
#summary(nb1_w)

#second spatial model
#sem_error_nb2_w<-errorsarlm(lm_frugi, listw=nb2_w)
#moran.test(residuals(sem_error_nb2_w), listw=nb_moran1000, randomisation=TRUE, zero.policy=NULL, alternative="greater")                                   


#################################################################################################
# PE 28 (= sensitivity analysis with all bio 12 contributions >20 removed)
#################################################################################################
rm(list=ls())
library(SDMTools)
library(raster)
#library(relaimpo)
library(glmulti)

setwd("/Users/Chris/Desktop/glmulti_SARs/")

# first a function - main script follows below
glm_process = function( result_frame,
                        model_number,
                        predictor_text,
                        response_text,        
                        response_text_logistic = "div", # the response variable name, not values                        the_data,
                        do_spatial_LM = TRUE,
                        weight_list = NULL,
                        do_logistic = TRUE,
                        the_data
)
{
  cat("\n\n** ",model_number,result_frame[model_number,1],result_frame[model_number,2],result_frame[model_number,3]," **\n")
  
  the.formula <- paste(response_text,"~",predictor_text)
  glm_gauss <- glm(the.formula, data=the_data, family="gaussian")
  print(summary(glm_gauss))
  cat("\nGaussian model\naic:",glm_gauss$aic)
  dev_exp <- (glm_gauss$null.deviance-glm_gauss$deviance)/glm_gauss$null.deviance
  cat("\nDeviance explained:",dev_exp,"\n")
  result_frame[model_number,"glm_r2"] <- round(dev_exp,4)
  result_frame[model_number,"glm_aic"] <- round(glm_gauss$aic,1)
  
  if (do_logistic) {
    the.formula <- paste(response_text_logistic,"~",predictor_text)    
    glm_logistic <- glm(the.formula, data=the_data, family = "binomial")
    print(summary(glm_logistic))
    cat("\nLogistic model \naic:",glm_logistic$aic)
    dev_exp <- (glm_logistic$null.deviance-glm_logistic$deviance)/glm_logistic$null.deviance
    cat("\nDeviance explained:",dev_exp,"\n")
    result_frame[model_number,"glm_logistic_r2"] <- round(dev_exp,4)
    result_frame[model_number,"glm_logistic_aic"] <- round(glm_logistic$aic,1)
  }
  
  if (do_spatial_LM) {
    
    SAR_gauss <- errorsarlm(glm_gauss, data=the_data, listw=weight_list, quiet=FALSE, na.omit, zero.policy=TRUE, tol.solve=1e-11)    
    #calculate AIC
    SAR_gauss.aic <- (-2*SAR_gauss$LL)+(2*SAR_gauss$parameters)
    
    print(summary(SAR_gauss))
    cat("\nSAR gaussian \naic:",SAR_gauss.aic)
    result_frame[model_number,"sarlm_aic"] <- round(SAR_gauss.aic,1)
    
    if (do_logistic) {
      SAR_logistic <- errorsarlm(glm_logistic, data=the_data, listw=weight_list, quiet=FALSE, na.omit, zero.policy=TRUE, tol.solve=1e-11)
      #calculate AIC
      SAR_logistic.aic <- (-2*SAR_logistic$LL)+(2*SAR_logistic$parameters)
      
      print(summary(SARer_with_exp))
      cat("\nSARLM logistic \naic:",SAR_logistic.aic)
      result_frame[model_number,10] <- round(SAR_logistic.aic,1)
    }
    
  }
  
  return(result_frame)
}



base_path     <- "/Users/Chris/Desktop/glmulti_SARs/"
results.dir   <- paste(base_path,sep='')

regions <- data.frame()
do_logistic <- TRUE
logistic_threshold <- 0.95  # for a logistic model to predict membership of top diversity

#whether to do spatial autocorrelation, and at what radius
do_spatial_LM <- FALSE
spatial_LM_radius <- 3000

# whether to do delta AIC table and glmulti
do_glmulti <- TRUE
do_plots   <- FALSE

i <- 1
regions[i,"region"]        <- ''
regions[i,"current_forest"] <- paste(regions[i],"/current_forest_clipped.asc",sep='')

diversities <- data.frame(taxon="PE_28",
                          level="lineage", 
                          metric="endemism", 
                          grid="/phylogenetic_endemism28_clipped.asc",
                          transform="log",
                          stringsAsFactors = F)
j <- 2
diversities[j,1:4] <- c("PE_35","lineage","endemism","/phylogenetic_endemism28_clipped.asc")

#output <- data.frame()
k  <- 0

# create a data frame to hold the model results
result_frame <- result_frame <- data.frame(region="",response_description="",predictor_description="",n=0,glm_r2=0,glm_aic=0,glm_delta_aic=0,sarlm_aic=0,sarlm_delta_aic=0,glm_logistic_r2=0,glm_logistic_aic=0,glm_logistic_delta_aic=0,sarlm_logistic_aic=0,sarlm_logistic_delta_aic=0,stringsAsFactors = F)

if (do_spatial_LM) {
  library(spdep)
  
  
}

if (do_glmulti) {
  # create a result data frame
  glmulti_blank_results <- data.frame(region="",max_pred=0,AIC=-9999,deltaAIC=-9999,r2=0,formula="",stringsAsFactors = F)
  
  # add predictor columns
  glmulti_blank_results <- cbind(glmulti_blank_results,data.frame(current_forest=0,bio1=0,bio4=0,bio12=0,bio14=0,foreststability_120k=0,topographic_heterogeneity=0,anom_bio1=0,anom_bio12=0,stringsAsFactors = 0))
  glmulti_full_results <- glmulti_blank_results[-1,]
}

#Loop through each region, and within it, each diversity metric
for (i in 1:nrow(regions)) {
  
  # for the all region model, assign region sizes to cells
  
  #define some basic data
  rf.ras <- raster("current_forest_clipped.asc")# read in the vegetation grid
  rf.ras[which(is.finite(rf.ras[]) & rf.ras[] !=1)] <- 0  #set all veg != 1 (rainforests) to 0
  #rf_region.asc     <- read.asc(regions$region[i])
  current_forest.ras <- raster("current_forest_clipped.asc")
  
  #crop the rainforest raster to match the stability raster
  rf.ras = crop(x=rf.ras,y=current_forest.ras)
  rf.asc = asc.from.raster(rf.ras)
  current_forest.asc = asc.from.raster(current_forest.ras)
  
  j=1 # for now only do lineage endemism
  #for (j in 1:nrow(diversities)) {
  
  #load the diversity result at 0.01 degree resolution and resample to match stability
  div.ras     <-  raster(paste(results.dir,diversities$grid[j],sep=''))
  div_resample.ras <- resample(div.ras,current_forest.ras,method="bilinear")
  div_resample.asc <- asc.from.raster(div_resample.ras)
  
  #create pos_rf the table of grid cell values for rainforest
  pos_rf <- as.data.frame(which(is.finite(rf.asc),arr.ind=TRUE)) #get all points that have data
  
  xy_rf <- getXYcoords(rf.asc)
  pos_rf$long  <- xy_rf$x[pos_rf$row]
  pos_rf$lat  <- xy_rf$y[pos_rf$col]
  rm(xy_rf)
  
  pos_rf$rf = rf.asc[cbind(pos_rf$row,pos_rf$col)]            #append the current forest data
  #pos_rf = pos_rf[which(pos_rf$rf==1),]  # filter to rainforest areas
  pos_rf$div = div_resample.asc[cbind(pos_rf$row,pos_rf$col)] #append the PE data
  pos_rf = pos_rf[which(pos_rf$div > 0),]  # filter to areas with a PE score
  
  # add the predictors based on rainforest niche models
  pos_rf$current_forest = current_forest.asc[cbind(pos_rf$row,pos_rf$col)] #append the current forest scores
  
  # add other predictors
  predictors <- data.frame(name="bio1",description="bio1",path="bio1_clipped.asc",resample=FALSE,stringsAsFactors = F)
  predictors[2,] <- c("bio4","bio4","bio4_clipped.asc",FALSE)
  predictors[3,] <- c("bio12","bio12","bio12_clipped.asc",FALSE)
  predictors[4,] <- c("bio14","bio14","bio14_clipped.asc",FALSE)
  predictors[5,] <- c("topographic_heterogeneity","topographic_heterogeneity","topographic_heterogeneity_clipped.asc",FALSE)
  predictors[6,] <- c("foreststability_120k","foreststability_120k","120k_foreststability_clipped.asc",FALSE)
  predictors[7,] <- c("anom_bio1","anom_bio1","anom_bio1_clipped.asc",FALSE)
  predictors[8,] <- c("anom_bio12","anom_bio12","anom_bio12_clipped.asc",FALSE)
  
  
  # patch size
  
  for (p in 1:nrow(predictors) ) {
    path <- predictors$path[p]
    env.ras <- raster(path)
    if (predictors$resample[p]) {
      env.ras <- resample(env.ras,current_forest.ras,method="bilinear")
      cat(predictors[p,"name"], "done\n")
    }
    
    pos_rf[,predictors$name[p]] <- extract(env.ras,pos_rf[,c("long","lat")])
  }
  
  # transform the response variable if needed (as set in the diversities data frame)
  if (diversities[j,"transform"] == "log") {
    pos_rf$div <- log(pos_rf$div)
  }
  
  # rescale the predictors - KEEP THE COLUMN NUMBERS CORRECT
  #pos_rf[,7:15] <- scale(pos_rf[,7:15])
  pos_rf[,7:15] <- scale(pos_rf[,7:15])
  
  # add a binary variable for membership of top class
  thresh <- quantile(pos_rf$div,logistic_threshold)
  pos_rf$div_binary <- rep(0,nrow(pos_rf))
  pos_rf$div_binary[pos_rf$div >= thresh] <- 1
  
  # calculate the distance matrix for SARLM    
  if (do_spatial_LM) {
    cat("\n\nPreparing neighbourhood weights for SARLM\n\n")
    coords<-as.matrix(cbind(pos_rf$row, pos_rf$col))
    cont.nb <- dnearneigh(coords,0,spatial_LM_radius,longlat=TRUE)
    weight_list <- nb2listw(cont.nb, glist=NULL, style="W", zero.policy=TRUE)
  }
  
  n <- nrow(pos_rf)
  
  if (do_glmulti) {
    
    # automated predictor selection
    
    formula <- "div~foreststability_120k+current_forest+topographic_heterogeneity+bio1+bio4+bio12+bio14+anom_bio1+anom_bio12" 
    
  }
  fullGLM      <- glm(data=pos_rf, formula=formula)
  region_startrow <- nrow(glmulti_full_results) + 1  # use this to calculate deltaAIC for just the region
  
  for (model_size in 1:9) {
    rm(glmulti_res)
    
    model_name <- paste(diversities[j,1],diversities[j,2],diversities[j,3],"for",regions[i,"region"])
    
    try(glmulti_res <- glmulti(y=fullGLM,level=1, maxsize=model_size, crit="aicc", method="h",name=model_name, report=T, confsetsize=512, plotty=F))
    if (exists("glmulti_res")) {
      
      glmulti_sum <- summary(glmulti_res)
      
      #####
      # fit the best model
      bestGLM <- glm(data=pos_rf, formula=glmulti_sum$bestmodel)
      bestGLM_r2 <- (bestGLM$null.deviance - bestGLM$deviance) / bestGLM$null.deviance
      bestGLM_r2 <- round(bestGLM_r2,3)
      bestGLM_sum <- summary(bestGLM)
      
      
      # store results in the data frame
      glmulti_new_results <- glmulti_blank_results
      glmulti_new_results$region[1]      <- regions$region[i]
      glmulti_new_results$max_pred[1]    <- model_size
      glmulti_new_results$AIC[1]         <- glmulti_sum$bestic
      glmulti_new_results$formula[1]     <- glmulti_sum$bestmodel
      glmulti_new_results$r2[1]          <- bestGLM_r2
      
      # get model terms and co-efficents
      terms <- attributes(bestGLM_sum$terms)$term.labels
      coef  <- bestGLM_sum$coefficients
      for (pred in terms) {
        glmulti_new_results[1,pred] <- round(coef[pred,"Estimate"],4)
      }
      
      glmulti_full_results <- rbind(glmulti_full_results,glmulti_new_results)
      
      cat("\n**********************************\nBest",model_size,"predictor model of",diversities[j,1],diversities[j,2],diversities[j,3],"for",regions[i,"region"])
      cat("\nBest model formula:",glmulti_sum$bestmodel,
          "\nBest model AICC:   ",glmulti_sum$bestic,
          "\nBest model r^2:    ",bestGLM_r2,"\n**********************************\n")
      
      print(summary(bestGLM))
    } else {
      glmulti_new_results             <- glmulti_blank_results
      glmulti_new_results$formula[1]  <- "FAILED"
      
      cat("\n**********************************\nBest",model_size,"predictor model of",diversities[j,1],diversities[j,2],diversities[j,3],"for",regions[i,"region"],
          "\n **  FAILED  **",
          "\n**********************************\n")
      
    }
  }
  
  # calculate deltaAIC
  rows    <- region_startrow:nrow(glmulti_full_results)
  minAIC  <- min(glmulti_full_results$AIC,na.rm=T)
  glmulti_full_results$deltaAIC[rows] <- glmulti_full_results$AIC[rows] - minAIC
}
}

if (do_glmulti) {
  setwd(results.dir)
  write.csv(glmulti_full_results,"Result_table_PE_28.csv")
  rm(rf.ras, rf.asc, div.ras, div_resample.ras, current_forest.asc, current_forest.ras,glmulti_blank_results,glmulti_new_results, rows,minAIC)
}

if (do_plots) {
  # now plot current suitability v stability (past suitability), coloured by endemism
  library(maptools)
  library(classInt)
  #windows()
  class_count <- 12
  my.class <- classIntervals(pos_rf$div,n=class_count,style="quantile", digits=2)
  my.class_breaks <- round(my.class[[2]],4)
  my.pal <- c("darkblue","green2","yellow","red")
  my.col <-findColours(my.class,my.pal)
  legend_cols <- attr(my.col,"palette")
  plot(pos_rf$foreststability_120k,pos_rf$current_forest,xlab="Current suitability",ylab="Stability",col=my.col)
  legend(x="topleft",legend=my.class_breaks[1:class_count+1],fill=legend_cols)
}

#write the weights

tmp <- weightable (glmulti_res)
tmp
write.csv(tmp,"Weight_table_PE_28.csv")
coef.glmulti(glmulti_res) # FROM DAN

############################################################################
#####################    SPATIAL    ########################################
############################################################################
############################################################################

#Load packages
library(spdep)
library(akima)
library(raster)


#----------------------------- PE DATASET -----------------------------------#

#########################
#Examine the structure of the PE dataset
PE <- pos_rf
names(PE)
str(PE)
#PE[complete.cases(PE[,7:11]),]

#Make maps
par(mfrow=c(4,3))
image(interp(PE$long, PE$lat, PE$div), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="PE")
image(interp(PE$long, PE$lat, PE$current_forest), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="current_forest")
image(interp(PE$long, PE$lat, PE$foreststability_120k), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="foreststability_120k")
image(interp(PE$long, PE$lat, PE$topographic_heterogeneity), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="topographic_heterogeneity")
image(interp(PE$long, PE$lat, PE$bio1), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="annual mean temp")
image(interp(PE$long, PE$lat, PE$bio12), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="annnual precipitation")
image(interp(PE$long, PE$lat, PE$bio14), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="precipitation of driest quarter")
image(interp(PE$long, PE$lat, PE$bio4), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="tenperature seasonality")
image(interp(PE$long, PE$lat, PE$anom_bio1), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="Quaternary temp change")
image(interp(PE$long, PE$lat, PE$anom_bio12), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="Quaternary precipitation change")

#########################
#Make a regression model to explain PE
hist(PE$div)
hist(log(PE$div))
formula <- "div~foreststability_120k+topographic_heterogeneity+bio1+bio4+bio12+bio14+anom_bio1+anom_bio12" 
lm_PE <- glm(data=pos_rf, formula=formula)
summary(lm_PE)
hist(residuals(lm_PE))

#########################
#Spatial structure of residuals
library(ncf)

par(mfrow=c(2,2))
#Correlograms with latlon = FALSE
cor.OBL<-correlog(PE$long, PE$lat, z=PE$div, na.rm=T, increment=0.5, resamp=1, latlon = FALSE)    #makes 10 distance classes
cor.OBL
plot(cor.OBL$correlation, type="b", pch=16, cex=1.2, lwd=1.5, ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)
title(main="0.5 degree")

#With latlon = TRUE
cor.OBL_latlon<-correlog(PE$long, PE$lat, z=PE$div, na.rm=T, increment=0.5, resamp=1, latlon = TRUE)  #uses km because latlon = TRUE
plot(cor.OBL_latlon$correlation, type="b", pch=16, cex=1.2, lwd=1.5, ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)    
title(main="0.5km")

#With latlon = TRUE and increment=10
cor.OBL_10<-correlog(PE$long, PE$lat, z=PE$div, na.rm=T, increment=10, resamp=1, latlon = TRUE)  #uses km because latlon = TRUE
plot(cor.OBL_10$correlation, type="b", pch=16, cex=1.2, lwd=1.5, ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)
title(main="10km")

summary(cor.OBL_10)
cor.OBL_10$mean.of.class
cor.OBL_10$n
cor.OBL_10$correlation      
cor.OBL_10$p                

#Correlogram for residuals
cor.res_10<-correlog(PE$long, PE$lat, z=residuals(lm_PE), na.rm=T, increment=10, resamp=999, latlon = TRUE)  #uses km because latlon = TRUE

#Plot both residuals and raw data
plot(cor.OBL_10$correlation, type="b", pch=16, col="black", cex=1.2, lwd=1.5, xlim=c(0,20), ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)                  
points(cor.res_10$correlation, pch=1, cex=1.2, col="blue")
lines(cor.res_10$correlation, lwd=1.5, col="blue")
title(main="10km raw [black] vs residuals [blue]")

#########################
#Implementing a spatial model

#Make coordinate list
coords_PE<-as.matrix(cbind(PE$long,PE$lat))
plot(coords_PE)

#Minimum distance to connect to at least one neighbor
PE_knear <- knn2nb(knearneigh(coords_PE, k=1))
summary(PE_knear)
dsts_PE<-unlist(nbdists(PE_knear, coords_PE, longlat = TRUE))
summary(dsts_PE)
max(dsts_PE)

#Neighbour defined by dnearneigh()
PE_nb<-dnearneigh(coords_PE,0,max(dsts_PE), longlat=T)
par(mfrow=c(1,1))
plot(PE_nb, coords_PE, pch=20, lwd=2)
summary(PE_nb)

PE_nb_10<-dnearneigh(coords_PE,0,10, longlat=T) #max(dsts_PE)=12.99333 #use 8? (Moran's I = 0 there)
summary(PE_nb_10)
plot(PE_nb_10, coords_PE, pch=20, lwd=2)

#Defining the spatial weights matrix
nb1_w<-nb2listw(PE_nb, glist=NULL, style="W", zero.policy=TRUE)
summary(nb1_w)
plot(nb1_w, coords_PE)

#Spatial autoregressive error model
formula <- "div~foreststability_120k+bio1+bio4+bio12+bio14+topographic_heterogeneity+anom_bio1+anom_bio12" 
the_data <- PE
sem_error_nb1_w<-errorsarlm(formula, data=the_data, listw=nb1_w, zero.policy=TRUE)
summary(sem_error_nb1_w)
errorsalm_summary <- summary(sem_error_nb1_w)
capture.output(errorsalm_summary, file ="SAR_summary_28.txt")

#########################
#Testing for spatial structure in the residuals of the spatial model

#Correlogram for spatial model
cor.SEM_10<-correlog(PE$long, PE$lat, z=residuals(sem_error_nb1_w), increment=10, resamp=999, na.rm=T, latlon = TRUE)  #uses km because latlon = TRUE

#Plot both residuals and raw data
plot(cor.OBL_10$correlation, type="b", pch=16, cex=1.2, lwd=1.5, xlim=c(0,9), ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)                  
points(cor.res_10$correlation, pch=1, cex=1.2, col="blue")
lines(cor.res_10$correlation, lwd=1.5, col="blue")
lines(cor.SEM_10$correlation, lwd=2.5, col="red")
title(main="10km raw [black] vs residuals [blue] vs SEM model [red]")

#Calculate Moran's I value individually
#PE_nb_moran10<-dnearneigh(coords_PE,0,10, longlat=T)
#nb_moran10<-nb2listw(PE_nb_moran10, glist=NULL, style="W", zero.policy=TRUE)

#moran.test(PE$div, listw=nb_moran10, randomisation=TRUE, zero.policy=TRUE, alternative="two.sided")
#moran.test(residuals(lm_PE2), listw=nb_moran10, randomisation=TRUE, zero.policy=TRUE, alternative="two.sided")
#moran.test(residuals(sem_error_nb1_w), listw=nb_moran10, randomisation=TRUE, zero.policy=NULL, alternative="two.sided")

#moran.test(PE$div, listw=nb_moran10, randomisation=TRUE, zero.policy=TRUE, alternative="greater")
#moran.test(residuals(lm_frugi), listw=nb_moran1000, randomisation=TRUE, zero.policy=NULL, alternative="greater")
#moran.test(residuals(sem_error_nb1_w), listw=nb_moran1000, randomisation=TRUE, zero.policy=NULL, alternative="greater")                                   

#Residual maps
#par(mfrow=c(2,2))
#image(interp(birds$X, birds$Y, birds$Frugivores), col = terrain.colors(12), main="Frugivore richness")
#image(interp(birds$X, birds$Y, residuals(lm_frugi)), col = terrain.colors(12), main="OLS model residuals")
#image(interp(birds$X, birds$Y, residuals(sem_error_nb1_w)), col = terrain.colors(12), main="SEM residuals")

#plot(coords_birds, col=c("blue", "red")[sign(residuals(lm_frugi))/2+1.5], pch=19,
#cex=abs(residuals(lm_frugi))/max(residuals(lm_frugi)), xlab="geographical x-coordinates", ylab="geographical y-coordinates")

#plot(coords_birds, col=c("blue", "red")[sign(residuals(sem_error_nb1_w))/2+1.5], pch=19,
#    cex=abs(residuals(sem_error_nb1_w))/max(residuals(sem_error_nb1_w)), xlab="geographical x-coordinates", ylab="geographical y-coordinates")

#########################
#Improving the spatial model

#Second spatial weights matrix
#birds_nb_200<-dnearneigh(coords_birds,0,200, longlat=T)
#summary(birds_nb_200)
#nb2_w<-nb2listw(birds_nb_200, glist=NULL, style="W", zero.policy=FALSE)
#summary(nb2_w)
#summary(nb1_w)

#second spatial model
#sem_error_nb2_w<-errorsarlm(lm_frugi, listw=nb2_w)
#moran.test(residuals(sem_error_nb2_w), listw=nb_moran1000, randomisation=TRUE, zero.policy=NULL, alternative="greater")     


#################################################################################################
# PE 22 (= sensitivity analysis with all bio 14 contributions >20 removed)
#################################################################################################
rm(list=ls())
library(SDMTools)
library(raster)
#library(relaimpo)
library(glmulti)

setwd("/Users/Chris/Desktop/glmulti_SARs/")

# first a function - main script follows below
glm_process = function( result_frame,
                        model_number,
                        predictor_text,
                        response_text,        
                        response_text_logistic = "div", # the response variable name, not values                        the_data,
                        do_spatial_LM = TRUE,
                        weight_list = NULL,
                        do_logistic = TRUE,
                        the_data
)
{
  cat("\n\n** ",model_number,result_frame[model_number,1],result_frame[model_number,2],result_frame[model_number,3]," **\n")
  
  the.formula <- paste(response_text,"~",predictor_text)
  glm_gauss <- glm(the.formula, data=the_data, family="gaussian")
  print(summary(glm_gauss))
  cat("\nGaussian model\naic:",glm_gauss$aic)
  dev_exp <- (glm_gauss$null.deviance-glm_gauss$deviance)/glm_gauss$null.deviance
  cat("\nDeviance explained:",dev_exp,"\n")
  result_frame[model_number,"glm_r2"] <- round(dev_exp,4)
  result_frame[model_number,"glm_aic"] <- round(glm_gauss$aic,1)
  
  if (do_logistic) {
    the.formula <- paste(response_text_logistic,"~",predictor_text)    
    glm_logistic <- glm(the.formula, data=the_data, family = "binomial")
    print(summary(glm_logistic))
    cat("\nLogistic model \naic:",glm_logistic$aic)
    dev_exp <- (glm_logistic$null.deviance-glm_logistic$deviance)/glm_logistic$null.deviance
    cat("\nDeviance explained:",dev_exp,"\n")
    result_frame[model_number,"glm_logistic_r2"] <- round(dev_exp,4)
    result_frame[model_number,"glm_logistic_aic"] <- round(glm_logistic$aic,1)
  }
  
  if (do_spatial_LM) {
    
    SAR_gauss <- errorsarlm(glm_gauss, data=the_data, listw=weight_list, quiet=FALSE, na.omit, zero.policy=TRUE, tol.solve=1e-11)    
    #calculate AIC
    SAR_gauss.aic <- (-2*SAR_gauss$LL)+(2*SAR_gauss$parameters)
    
    print(summary(SAR_gauss))
    cat("\nSAR gaussian \naic:",SAR_gauss.aic)
    result_frame[model_number,"sarlm_aic"] <- round(SAR_gauss.aic,1)
    
    if (do_logistic) {
      SAR_logistic <- errorsarlm(glm_logistic, data=the_data, listw=weight_list, quiet=FALSE, na.omit, zero.policy=TRUE, tol.solve=1e-11)
      #calculate AIC
      SAR_logistic.aic <- (-2*SAR_logistic$LL)+(2*SAR_logistic$parameters)
      
      print(summary(SARer_with_exp))
      cat("\nSARLM logistic \naic:",SAR_logistic.aic)
      result_frame[model_number,10] <- round(SAR_logistic.aic,1)
    }
    
  }
  
  return(result_frame)
}



base_path     <- "/Users/Chris/Desktop/glmulti_SARs/"
results.dir   <- paste(base_path,sep='')

regions <- data.frame()
do_logistic <- TRUE
logistic_threshold <- 0.95  # for a logistic model to predict membership of top diversity

#whether to do spatial autocorrelation, and at what radius
do_spatial_LM <- FALSE
spatial_LM_radius <- 3000

# whether to do delta AIC table and glmulti
do_glmulti <- TRUE
do_plots   <- FALSE

i <- 1
regions[i,"region"]        <- ''
regions[i,"current_forest"] <- paste(regions[i],"/current_forest_clipped.asc",sep='')

diversities <- data.frame(taxon="PE_22",
                          level="lineage", 
                          metric="endemism", 
                          grid="/phylogenetic_endemism22_clipped.asc",
                          transform="log",
                          stringsAsFactors = F)
j <- 2
diversities[j,1:4] <- c("PE_35","lineage","endemism","/phylogenetic_endemism22_clipped.asc")

#output <- data.frame()
k  <- 0

# create a data frame to hold the model results
result_frame <- result_frame <- data.frame(region="",response_description="",predictor_description="",n=0,glm_r2=0,glm_aic=0,glm_delta_aic=0,sarlm_aic=0,sarlm_delta_aic=0,glm_logistic_r2=0,glm_logistic_aic=0,glm_logistic_delta_aic=0,sarlm_logistic_aic=0,sarlm_logistic_delta_aic=0,stringsAsFactors = F)

if (do_spatial_LM) {
  library(spdep)
  
  
}

if (do_glmulti) {
  # create a result data frame
  glmulti_blank_results <- data.frame(region="",max_pred=0,AIC=-9999,deltaAIC=-9999,r2=0,formula="",stringsAsFactors = F)
  
  # add predictor columns
  glmulti_blank_results <- cbind(glmulti_blank_results,data.frame(current_forest=0,bio1=0,bio4=0,bio12=0,bio14=0,foreststability_120k=0,topographic_heterogeneity=0,anom_bio1=0,anom_bio12=0,stringsAsFactors = 0))
  glmulti_full_results <- glmulti_blank_results[-1,]
}

#Loop through each region, and within it, each diversity metric
for (i in 1:nrow(regions)) {
  
  # for the all region model, assign region sizes to cells
  
  #define some basic data
  rf.ras <- raster("current_forest_clipped.asc")# read in the vegetation grid
  rf.ras[which(is.finite(rf.ras[]) & rf.ras[] !=1)] <- 0  #set all veg != 1 (rainforests) to 0
  #rf_region.asc     <- read.asc(regions$region[i])
  current_forest.ras <- raster("current_forest_clipped.asc")
  
  #crop the rainforest raster to match the stability raster
  rf.ras = crop(x=rf.ras,y=current_forest.ras)
  rf.asc = asc.from.raster(rf.ras)
  current_forest.asc = asc.from.raster(current_forest.ras)
  
  j=1 # for now only do lineage endemism
  #for (j in 1:nrow(diversities)) {
  
  #load the diversity result at 0.01 degree resolution and resample to match stability
  div.ras     <-  raster(paste(results.dir,diversities$grid[j],sep=''))
  div_resample.ras <- resample(div.ras,current_forest.ras,method="bilinear")
  div_resample.asc <- asc.from.raster(div_resample.ras)
  
  #create pos_rf the table of grid cell values for rainforest
  pos_rf <- as.data.frame(which(is.finite(rf.asc),arr.ind=TRUE)) #get all points that have data
  
  xy_rf <- getXYcoords(rf.asc)
  pos_rf$long  <- xy_rf$x[pos_rf$row]
  pos_rf$lat  <- xy_rf$y[pos_rf$col]
  rm(xy_rf)
  
  pos_rf$rf = rf.asc[cbind(pos_rf$row,pos_rf$col)]            #append the current forest data
  #pos_rf = pos_rf[which(pos_rf$rf==1),]  # filter to rainforest areas
  pos_rf$div = div_resample.asc[cbind(pos_rf$row,pos_rf$col)] #append the PE data
  pos_rf = pos_rf[which(pos_rf$div > 0),]  # filter to areas with a PE score
  
  # add the predictors based on rainforest niche models
  pos_rf$current_forest = current_forest.asc[cbind(pos_rf$row,pos_rf$col)] #append the current forest scores
  
  # add other predictors
  predictors <- data.frame(name="bio1",description="bio1",path="bio1_clipped.asc",resample=FALSE,stringsAsFactors = F)
  predictors[2,] <- c("bio4","bio4","bio4_clipped.asc",FALSE)
  predictors[3,] <- c("bio12","bio12","bio12_clipped.asc",FALSE)
  predictors[4,] <- c("bio14","bio14","bio14_clipped.asc",FALSE)
  predictors[5,] <- c("topographic_heterogeneity","topographic_heterogeneity","topographic_heterogeneity_clipped.asc",FALSE)
  predictors[6,] <- c("foreststability_120k","foreststability_120k","120k_foreststability_clipped.asc",FALSE)
  predictors[7,] <- c("anom_bio1","anom_bio1","anom_bio1_clipped.asc",FALSE)
  predictors[8,] <- c("anom_bio12","anom_bio12","anom_bio12_clipped.asc",FALSE)
  
  
  # patch size
  
  for (p in 1:nrow(predictors) ) {
    path <- predictors$path[p]
    env.ras <- raster(path)
    if (predictors$resample[p]) {
      env.ras <- resample(env.ras,current_forest.ras,method="bilinear")
      cat(predictors[p,"name"], "done\n")
    }
    
    pos_rf[,predictors$name[p]] <- extract(env.ras,pos_rf[,c("long","lat")])
  }
  
  # transform the response variable if needed (as set in the diversities data frame)
  if (diversities[j,"transform"] == "log") {
    pos_rf$div <- log(pos_rf$div)
  }
  
  # rescale the predictors - KEEP THE COLUMN NUMBERS CORRECT
  #pos_rf[,7:15] <- scale(pos_rf[,7:15])
  pos_rf[,7:15] <- scale(pos_rf[,7:15])
  
  # add a binary variable for membership of top class
  thresh <- quantile(pos_rf$div,logistic_threshold)
  pos_rf$div_binary <- rep(0,nrow(pos_rf))
  pos_rf$div_binary[pos_rf$div >= thresh] <- 1
  
  # calculate the distance matrix for SARLM    
  if (do_spatial_LM) {
    cat("\n\nPreparing neighbourhood weights for SARLM\n\n")
    coords<-as.matrix(cbind(pos_rf$row, pos_rf$col))
    cont.nb <- dnearneigh(coords,0,spatial_LM_radius,longlat=TRUE)
    weight_list <- nb2listw(cont.nb, glist=NULL, style="W", zero.policy=TRUE)
  }
  
  n <- nrow(pos_rf)
  
  if (do_glmulti) {
    
    # automated predictor selection
    
    formula <- "div~foreststability_120k+current_forest+topographic_heterogeneity+bio1+bio4+bio12+bio14+anom_bio1+anom_bio12" 
    
  }
  fullGLM      <- glm(data=pos_rf, formula=formula)
  region_startrow <- nrow(glmulti_full_results) + 1  # use this to calculate deltaAIC for just the region
  
  for (model_size in 1:9) {
    rm(glmulti_res)
    
    model_name <- paste(diversities[j,1],diversities[j,2],diversities[j,3],"for",regions[i,"region"])
    
    try(glmulti_res <- glmulti(y=fullGLM,level=1, maxsize=model_size, crit="aicc", method="h",name=model_name, report=T, confsetsize=512, plotty=F))
    if (exists("glmulti_res")) {
      
      glmulti_sum <- summary(glmulti_res)
      
      #####
      # fit the best model
      bestGLM <- glm(data=pos_rf, formula=glmulti_sum$bestmodel)
      bestGLM_r2 <- (bestGLM$null.deviance - bestGLM$deviance) / bestGLM$null.deviance
      bestGLM_r2 <- round(bestGLM_r2,3)
      bestGLM_sum <- summary(bestGLM)
      
      
      # store results in the data frame
      glmulti_new_results <- glmulti_blank_results
      glmulti_new_results$region[1]      <- regions$region[i]
      glmulti_new_results$max_pred[1]    <- model_size
      glmulti_new_results$AIC[1]         <- glmulti_sum$bestic
      glmulti_new_results$formula[1]     <- glmulti_sum$bestmodel
      glmulti_new_results$r2[1]          <- bestGLM_r2
      
      # get model terms and co-efficents
      terms <- attributes(bestGLM_sum$terms)$term.labels
      coef  <- bestGLM_sum$coefficients
      for (pred in terms) {
        glmulti_new_results[1,pred] <- round(coef[pred,"Estimate"],4)
      }
      
      glmulti_full_results <- rbind(glmulti_full_results,glmulti_new_results)
      
      cat("\n**********************************\nBest",model_size,"predictor model of",diversities[j,1],diversities[j,2],diversities[j,3],"for",regions[i,"region"])
      cat("\nBest model formula:",glmulti_sum$bestmodel,
          "\nBest model AICC:   ",glmulti_sum$bestic,
          "\nBest model r^2:    ",bestGLM_r2,"\n**********************************\n")
      
      print(summary(bestGLM))
    } else {
      glmulti_new_results             <- glmulti_blank_results
      glmulti_new_results$formula[1]  <- "FAILED"
      
      cat("\n**********************************\nBest",model_size,"predictor model of",diversities[j,1],diversities[j,2],diversities[j,3],"for",regions[i,"region"],
          "\n **  FAILED  **",
          "\n**********************************\n")
      
    }
  }
  
  # calculate deltaAIC
  rows    <- region_startrow:nrow(glmulti_full_results)
  minAIC  <- min(glmulti_full_results$AIC,na.rm=T)
  glmulti_full_results$deltaAIC[rows] <- glmulti_full_results$AIC[rows] - minAIC
}
}

if (do_glmulti) {
  setwd(results.dir)
  write.csv(glmulti_full_results,"Result_table_PE_22.csv")
  rm(rf.ras, rf.asc, div.ras, div_resample.ras, current_forest.asc, current_forest.ras,glmulti_blank_results,glmulti_new_results, rows,minAIC)
}

if (do_plots) {
  # now plot current suitability v stability (past suitability), coloured by endemism
  library(maptools)
  library(classInt)
  #windows()
  class_count <- 12
  my.class <- classIntervals(pos_rf$div,n=class_count,style="quantile", digits=2)
  my.class_breaks <- round(my.class[[2]],4)
  my.pal <- c("darkblue","green2","yellow","red")
  my.col <-findColours(my.class,my.pal)
  legend_cols <- attr(my.col,"palette")
  plot(pos_rf$foreststability_120k,pos_rf$current_forest,xlab="Current suitability",ylab="Stability",col=my.col)
  legend(x="topleft",legend=my.class_breaks[1:class_count+1],fill=legend_cols)
}

#write the weights

tmp <- weightable (glmulti_res)
tmp
write.csv(tmp,"Weight_table_PE_22.csv")
coef.glmulti(glmulti_res) # FROM DAN

############################################################################
#####################    SPATIAL    ########################################
############################################################################
############################################################################

#Load packages
library(spdep)
library(akima)
library(raster)


#----------------------------- PE DATASET -----------------------------------#

#########################
#Examine the structure of the PE dataset
PE <- pos_rf
names(PE)
str(PE)
#PE[complete.cases(PE[,7:11]),]

#Make maps
par(mfrow=c(4,3))
image(interp(PE$long, PE$lat, PE$div), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="PE")
image(interp(PE$long, PE$lat, PE$current_forest), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="current_forest")
image(interp(PE$long, PE$lat, PE$foreststability_120k), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="foreststability_120k")
image(interp(PE$long, PE$lat, PE$topographic_heterogeneity), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="topographic_heterogeneity")
image(interp(PE$long, PE$lat, PE$bio1), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="annual mean temp")
image(interp(PE$long, PE$lat, PE$bio12), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="annnual precipitation")
image(interp(PE$long, PE$lat, PE$bio14), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="precipitation of driest quarter")
image(interp(PE$long, PE$lat, PE$bio4), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="tenperature seasonality")
image(interp(PE$long, PE$lat, PE$anom_bio1), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="Quaternary temp change")
image(interp(PE$long, PE$lat, PE$anom_bio12), col = terrain.colors(20), xlab="X coordinate", ylab="Y coordinate", main="Quaternary precipitation change")

#########################
#Make a regression model to explain PE
hist(PE$div)
hist(log(PE$div))
formula <- "div~foreststability_120k+topographic_heterogeneity+bio1+bio4+bio12+bio14+anom_bio1+anom_bio12" 
lm_PE <- glm(data=pos_rf, formula=formula)
summary(lm_PE)
hist(residuals(lm_PE))

#########################
#Spatial structure of residuals
library(ncf)

par(mfrow=c(2,2))
#Correlograms with latlon = FALSE
cor.OBL<-correlog(PE$long, PE$lat, z=PE$div, na.rm=T, increment=0.5, resamp=1, latlon = FALSE)    #makes 10 distance classes
cor.OBL
plot(cor.OBL$correlation, type="b", pch=16, cex=1.2, lwd=1.5, ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)
title(main="0.5 degree")

#With latlon = TRUE
cor.OBL_latlon<-correlog(PE$long, PE$lat, z=PE$div, na.rm=T, increment=0.5, resamp=1, latlon = TRUE)  #uses km because latlon = TRUE
plot(cor.OBL_latlon$correlation, type="b", pch=16, cex=1.2, lwd=1.5, ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)    
title(main="0.5km")

#With latlon = TRUE and increment=10
cor.OBL_10<-correlog(PE$long, PE$lat, z=PE$div, na.rm=T, increment=10, resamp=1, latlon = TRUE)  #uses km because latlon = TRUE
plot(cor.OBL_10$correlation, type="b", pch=16, cex=1.2, lwd=1.5, ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)
title(main="10km")

summary(cor.OBL_10)
cor.OBL_10$mean.of.class
cor.OBL_10$n
cor.OBL_10$correlation      
cor.OBL_10$p                

#Correlogram for residuals
cor.res_10<-correlog(PE$long, PE$lat, z=residuals(lm_PE), na.rm=T, increment=10, resamp=999, latlon = TRUE)  #uses km because latlon = TRUE

#Plot both residuals and raw data
plot(cor.OBL_10$correlation, type="b", pch=16, col="black", cex=1.2, lwd=1.5, xlim=c(0,20), ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)                  
points(cor.res_10$correlation, pch=1, cex=1.2, col="blue")
lines(cor.res_10$correlation, lwd=1.5, col="blue")
title(main="10km raw [black] vs residuals [blue]")

#########################
#Implementing a spatial model

#Make coordinate list
coords_PE<-as.matrix(cbind(PE$long,PE$lat))
plot(coords_PE)

#Minimum distance to connect to at least one neighbor
PE_knear <- knn2nb(knearneigh(coords_PE, k=1))
summary(PE_knear)
dsts_PE<-unlist(nbdists(PE_knear, coords_PE, longlat = TRUE))
summary(dsts_PE)
max(dsts_PE)

#Neighbour defined by dnearneigh()
PE_nb<-dnearneigh(coords_PE,0,max(dsts_PE), longlat=T)
par(mfrow=c(1,1))
plot(PE_nb, coords_PE, pch=20, lwd=2)
summary(PE_nb)

PE_nb_10<-dnearneigh(coords_PE,0,10, longlat=T) #max(dsts_PE)=12.99333 #use 8? (Moran's I = 0 there)
summary(PE_nb_10)
plot(PE_nb_10, coords_PE, pch=20, lwd=2)

#Defining the spatial weights matrix
nb1_w<-nb2listw(PE_nb, glist=NULL, style="W", zero.policy=TRUE)
summary(nb1_w)
plot(nb1_w, coords_PE)

#Spatial autoregressive error model
formula <- "div~foreststability_120k+bio1+bio4+bio12+bio14+topographic_heterogeneity+anom_bio1+anom_bio12" 
the_data <- PE
sem_error_nb1_w<-errorsarlm(formula, data=the_data, listw=nb1_w, zero.policy=TRUE)
summary(sem_error_nb1_w)
errorsalm_summary <- summary(sem_error_nb1_w)
capture.output(errorsalm_summary, file ="SAR_summary_22.txt")

#########################
#Testing for spatial structure in the residuals of the spatial model

#Correlogram for spatial model
cor.SEM_10<-correlog(PE$long, PE$lat, z=residuals(sem_error_nb1_w), increment=10, resamp=999, na.rm=T, latlon = TRUE)  #uses km because latlon = TRUE

#Plot both residuals and raw data
plot(cor.OBL_10$correlation, type="b", pch=16, cex=1.2, lwd=1.5, xlim=c(0,9), ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)                  
points(cor.res_10$correlation, pch=1, cex=1.2, col="blue")
lines(cor.res_10$correlation, lwd=1.5, col="blue")
lines(cor.SEM_10$correlation, lwd=2.5, col="red")
title(main="10km raw [black] vs residuals [blue] vs SEM model [red]")

#Calculate Moran's I value individually
#PE_nb_moran10<-dnearneigh(coords_PE,0,10, longlat=T)
#nb_moran10<-nb2listw(PE_nb_moran10, glist=NULL, style="W", zero.policy=TRUE)

#moran.test(PE$div, listw=nb_moran10, randomisation=TRUE, zero.policy=TRUE, alternative="two.sided")
#moran.test(residuals(lm_PE2), listw=nb_moran10, randomisation=TRUE, zero.policy=TRUE, alternative="two.sided")
#moran.test(residuals(sem_error_nb1_w), listw=nb_moran10, randomisation=TRUE, zero.policy=NULL, alternative="two.sided")

#moran.test(PE$div, listw=nb_moran10, randomisation=TRUE, zero.policy=TRUE, alternative="greater")
#moran.test(residuals(lm_frugi), listw=nb_moran1000, randomisation=TRUE, zero.policy=NULL, alternative="greater")
#moran.test(residuals(sem_error_nb1_w), listw=nb_moran1000, randomisation=TRUE, zero.policy=NULL, alternative="greater")                                   

#Residual maps
#par(mfrow=c(2,2))
#image(interp(birds$X, birds$Y, birds$Frugivores), col = terrain.colors(12), main="Frugivore richness")
#image(interp(birds$X, birds$Y, residuals(lm_frugi)), col = terrain.colors(12), main="OLS model residuals")
#image(interp(birds$X, birds$Y, residuals(sem_error_nb1_w)), col = terrain.colors(12), main="SEM residuals")

#plot(coords_birds, col=c("blue", "red")[sign(residuals(lm_frugi))/2+1.5], pch=19,
#cex=abs(residuals(lm_frugi))/max(residuals(lm_frugi)), xlab="geographical x-coordinates", ylab="geographical y-coordinates")

#plot(coords_birds, col=c("blue", "red")[sign(residuals(sem_error_nb1_w))/2+1.5], pch=19,
#    cex=abs(residuals(sem_error_nb1_w))/max(residuals(sem_error_nb1_w)), xlab="geographical x-coordinates", ylab="geographical y-coordinates")

#########################
#Improving the spatial model

#Second spatial weights matrix
#birds_nb_200<-dnearneigh(coords_birds,0,200, longlat=T)
#summary(birds_nb_200)
#nb2_w<-nb2listw(birds_nb_200, glist=NULL, style="W", zero.policy=FALSE)
#summary(nb2_w)
#summary(nb1_w)

#second spatial model
#sem_error_nb2_w<-errorsarlm(lm_frugi, listw=nb2_w)
#moran.test(residuals(sem_error_nb2_w), listw=nb_moran1000, randomisation=TRUE, zero.policy=NULL, alternative="greater")     
