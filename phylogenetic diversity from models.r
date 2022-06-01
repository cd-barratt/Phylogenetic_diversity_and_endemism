## this script calculates richness and endemism from modelled suitability surfaces
## it requires all the model grids to have the same extent
rm(list=ls())

library(SDMTools)
library(raster)
library(ape)
library(phylobase)
library(foreach)
library(doParallel)
library(ggplot2)

  source('C:/Users/Chris/Desktop/PE/PE_paper_GLMs/generate_data_for_GLMs_NEW_NEW/phylogenetic endemism.r')
################################################################################
#first define some functions

map_raster = function(raster, output_file, title) {
  
  p         <- rasterToPoints(raster)
  p         <- data.frame(p)
  names(p) <- c("x", "y", "Model")
  colour_gradient <- scale_fill_gradientn(colours = rainbow(15), values=p$model)
  colour_gradient <- scale_fill_gradient2(low="white", mid="yellow", high="red",
                                          limits=c(min(p$Model),max(p$Model)), midpoint=quantile(p$Model, 0.75), space='Lab')
  m <- ggplot(data=p) + geom_tile(aes(x, y, fill=Model)) + coord_equal() + labs(x=NULL, y=NULL) + colour_gradient
  
  # delete a previous file if needed
  if (file.exists(output_file)) {
    file.remove(output_file)
    cat("Previous", output_file, "removed\n")
  }
  
  m <- m + ggtitle(title)
  m <- m + theme(axis.title=element_text(face="bold", size="18"))
  m <- m + theme(axis.text=element_text(face="bold", size="14"))
  m <- m + theme(plot.title=element_text(face="bold", size="24"))
  m <- m + xlab("longitude") + ylab("latitude")
  
  png(output_file, width=image.width, height=image.height)
  print(m)
  dev.off()
  m <- NULL
}

################################################################################
################################################################################

max.rows        <- 10000000
core_count      <- 4 # number of cores to use for parallel steps
write_matrices  <- TRUE

# size in pixels for maps
image.width=1400
image.height=1400

#define directories
base.dir        <- 'C:/Users/Chris/Desktop/PE/PE_paper_GLMs/generate_data_for_GLMs_NEW_NEW/'      # modify to the base directory for your lineage modelling
input.dir       <- 'C:/Users/Chris/Desktop/PE/PE_paper_GLMs/generate_data_for_GLMs_NEW_NEW/asc_aligned/'
output.dir      <- 'C:/Users/Chris/Desktop/PE/PE_paper_GLMs/generate_data_for_GLMs_NEW_NEW/PE_output/'  # output location where diversity results and maps will be saved
file_pattern    <- 'lin_model_'                     # modify this to a match the start of the name of all lineage model asc files

template_grid   <- 'C:/Users/Chris/Desktop/PE/PE_paper_GLMs/generate_data_for_GLMs_NEW_NEW/template.asc.gz'
group_lin_file  <- 'C:/Users/Chris/Desktop/PE/PE_paper_GLMs/generate_data_for_GLMs_NEW_NEW/asc_aligned/All_lineage_list.csv'

preface         <- ""

genus           <- ''
output_prefix   <- 'All_17_species'
threshold       <- 0.000000000000001  # this is not a species level threshold, but one used for each lineage model

####  end of parameters  ####

setwd(base.dir)
files <- list.files(path = input.dir, pattern = file_pattern, recursive = FALSE,ignore.case = TRUE, include.dirs = FALSE)

setwd(input.dir)
#template.asc = read.asc.gz('files[1])
template.asc = read.asc('template.asc')
model_rows=nrow(template.asc)
model_cols=ncol(template.asc)

# the original version, excluding NA cells
pos <- as.data.frame(which(is.finite(template.asc),arr.ind=TRUE)) #get all points that have data

###
xy_template <-getXYcoords(template.asc)
pos$long  <- xy_template$x[pos$row]
pos$lat  <- xy_template$y[pos$col]
rm(xy_template)

####
cat("\nLoading model rasters in parallel\n")

cl <- makeCluster(core_count)
registerDoParallel(cl)

pos_par <- foreach (j=1:length(files), .combine=cbind, .packages='SDMTools') %dopar% {
  tfile <- files[j]
  pos_temp <- pos
  checkname = unlist(strsplit(tfile,".",fixed=T))
  if (checkname[length(checkname)]=="asc") {   # only accept filenames ending in .asc
    
    cat(j)
    
    tasc = read.asc(tfile)                                #read in the data
    dataname=gsub(".asc",'',tfile)
    newname <- tolower(gsub(preface, "", dataname))
    
    cat("About to do pos\n")
    
    pos_temp[newname] <- tasc[cbind(pos_temp$row, pos_temp$col)]
    pos_temp[(which(pos_temp[newname]< threshold)), newname] <- 0    # set values below the threshold to 0
    pos_temp[(which(is.na(pos_temp[newname]))), newname]     <- 0    # set the nulls to 0
    cat("\n", j, newname, "loaded")
    newcol <- data.frame(pos_temp[, newname])
    newcol <- round(newcol,3)
    names(newcol) <- newname
    newcol
  }
}

pos <- cbind(pos,pos_par)
rm(pos_par)

setwd(base.dir)
group_lin_list <-read.csv(group_lin_file, stringsAsFactors=F)
group_lin_list$lineage <- tolower(group_lin_list$lineage)
group_lin_list <- group_lin_list[group_lin_list$genus == genus,]



if (write_matrices) {
  cat("\nWriting the GLM data to file\n")
  write.csv(pos,paste(output.dir,"GLM_data.csv",sep=''))
}
