####################################################
## Fluorescent microscopy image analysis ###########
## with EBImage - optimized for 4x or 10x magn. ###
## images GFP-channel for fluorescence  ###########
## DAPI-channel for hoechst-stained nuclei ########
##### Chris Nowak March 2016 ######################
library(EBImage)
library(metafor)
dir <- "/Users/nowakchr/Desktop/Microscopy analysis in R /xxxxx"
setwd(dir)

# load image pairs: gfp/dapi 
## example: 4x magnification, one central picture/well
## wells B2-B4 basal, B5-B7 insulin
b2_1<-normalize(channel(readImage("B02_1.png"),"grey")) # normalize pixel values &
b3_1<-normalize(channel(readImage("B03_1.png"),"grey")) # convert to monocolour
b4_1<-normalize(channel(readImage("B04_1.png"),"grey"))
b5_1<-normalize(channel(readImage("B05_1.png"),"grey"))
b6_1<-normalize(channel(readImage("B06_1.png"),"grey"))
b7_1<-normalize(channel(readImage("B07_1.png"),"grey"))
b2_1h<-normalize(channel(readImage("B02_1h.png"),"grey")) # paired dapi-pictures
b3_1h<-normalize(channel(readImage("B03_1h.png"),"grey"))
b4_1h<-normalize(channel(readImage("B04_1h.png"),"grey"))
b5_1h<-normalize(channel(readImage("B05_1h.png"),"grey"))
b6_1h<-normalize(channel(readImage("B06_1h.png"),"grey"))
b7_1h<-normalize(channel(readImage("B07_1h.png"),"grey"))

l1<-list(b2_1,b3_1,b4_1,b5_1,b6_1,b7_1)
l2<-list(b2_1h,b3_1h,b4_1h,b5_1h,b6_1h,b7_1h)

# to inspect raw images
par(mfrow=c(2,1))
display(combine(b2_1*10,b2_1h*5),method="raster",all=T)
display(rgbImage(green=b2_1*10,red=b2_1h*5),method="raster",all=T)  # overlay image to check cytoplasm/nucleus registration

# create empty lists
nucleusThresh<-list()
nucleusOpening<-list()
nucleusSeed <-list()
nucleusFill<-list()
nucleusRegions<-list()
cytoplasmThresh<-list()
cytoplasmOpening<-list()
cytoplasmRegions<-list()
cytoplasmCombined<-list()
Features<-list()

# not needed for image analysis, but useful for comparisons:
# a function for counting nuclei in the dapi-images
count <- function(im){
  w=makeBrush(size=5,shape="gaussian",sigma=0.3)
  img_flo = filter2(im,w)
  display(img_flo*4,method="raster")
  nmaskt = thresh(img_flo,w=30,h=30,offset=quantile(im,0.5)) # offset needs adjusting for each experiment, see below
  display(nmaskt,method="raster")             
  nucNo<-max(bwlabel(nmaskt))
  cat("Number of nuclei =",max(bwlabel(nmaskt)),'\n')
  return(nucNo)
}

# background intensity values vary between experiments
# so for each new experiment, check a few images to make sure the algorithm
# correctly identifies cell bodies and nuclei
nucleusThresh[[2]] <- thresh(l2[[2]], w = 15, h = 15, offset = quantile(l2[[2]],0.5)) 
nucleusOpening[[2]] = opening(nucleusThresh[[2]], kern=makeBrush(3, shape="disc")) 
nucleusSeed[[2]] = bwlabel(nucleusOpening[[2]])  
nucleusFill[[2]] = fillHull(thresh(l1[[2]], w = 20, h = 20, offset = quantile(l1[[2]],0.05)))
nucleusRegions[[2]] = propagate(l1[[2]], nucleusSeed[[2]], mask=nucleusFill[[2]]) 
display(combine((nucleusRegions[[2]]),l2[[2]]*5),method="raster",all=T) # cell nuclei identification
  # usually works wo. adjustment as Hoechst-stained nuclei offset strongly from the background
  # otherwise adjust offset-threshold & makeBrush(x, ....)

cytoplasmThresh[[2]] = thresh(l1[[2]], w = 150, h = 150, quantile(l1[[2]],0.0000001)) # adjust offset...
  # rather than absolute pixel values, I use a quantile-thresholds per image, set very low
  # adjust the w/h-window according to magnification (i.e. should be large enough to ignore meaningless
  # bright dots, but small enough to capture cell bodies), over-doing it would be to size exactly according
  # to magnification to, e.g. 7um, if that is a meaningful min-cell size...
cytoplasmOpening[[2]] = opening(cytoplasmThresh[[2]],kern=makeBrush(3,shape="disc")) # ajust makeBrush(x, ) to 3/5/9/11
cytoplasmRegions[[2]] = propagate(l1[[2]], nucleusRegions[[2]], lambda=1.0e-04, mask=cytoplasmOpening[[2]]) 
par(mfrow=c(2,1))
display(colorLabels(cytoplasmOpening[[2]]),method="raster")
display((l1[[2]]*5),method="raster",all=T) # compare extracted cell bodies with original image & adjust
 
cytoplasmCombined[[2]] = cytoplasmOpening[[2]] # combine nucleus and cell body location images so
  # nuclei are used to propagate the recognition of a cell body area
cytoplasmCombined[[2]][nucleusFill[[2]] > cytoplasmCombined[[2]]] =nucleusFill[[2]][nucleusFill[[2]] > cytoplasmCombined[[2]]]
cytoplasmRegions[[2]] = propagate(l1[[2]], nucleusRegions[[2]], lambda=1.0e-4, mask=cytoplasmCombined[[2]])
display(colorLabels(cytoplasmRegions[[2]]),method="raster") # check that most cells are identified by different colours
display((l1[[2]]*5),method="raster",all=T)

# carry the optimized offset/kern/window-size values forward to loop over all images
for (i in 1:length(l1)){
  nucleusThresh[[i]] <- thresh(l2[[i]], w = 15, h = 15, offset = quantile(l2[[i]],0.5)) # ncl. segmentation by adaptive thresholding
  nucleusOpening[[i]] = opening(nucleusThresh[[i]], kern=makeBrush(3, shape="disc")) # morphological opening
  nucleusSeed[[i]] = bwlabel(nucleusOpening[[i]])  # labeling
  nucleusFill[[i]] = fillHull(thresh(l1[[i]], w = 20, h = 20, offset = quantile(l1[[i]],0.05))) # 2nd segmentation by mask to fill
  # foreground-pixel-surrounded "holes" with less stringent offset, usually not needed for the strong dapi-stain
  nucleusRegions[[i]] = propagate(l1[[i]], nucleusSeed[[i]], mask=nucleusFill[[i]]) # advanced segm. by Voronoi tessellation
  cytoplasmThresh[[i]] = thresh(l1[[i]], w = 150, h = 150, quantile(l1[[i]],0.0000001)) # dto. for cell bodies 
  cytoplasmOpening[[i]] = opening(cytoplasmThresh[[i]],kern=makeBrush(9,shape="disc"))
  cytoplasmRegions[[i]] = propagate(l1[[i]], nucleusRegions[[i]], lambda=1.0e-04, mask=cytoplasmOpening[[i]]) 
  cytoplasmCombined[[i]] = cytoplasmOpening[[i]] # define image area mask covered by cell bodies, combine cytoplasm
  # and nucleus masks by their union
  cytoplasmCombined[[i]][nucleusFill[[i]] > cytoplasmCombined[[i]]] =nucleusFill[[i]][nucleusFill[[i]] > cytoplasmCombined[[i]]]
  cytoplasmRegions[[i]] = propagate(l1[[i]], nucleusRegions[[i]], lambda=1.0e-4, mask=cytoplasmCombined[[i]]) # Voronoi tess.
}

# to produce a pretty picture that combines segm. results from nuclei/cell bodies
Im<-paintObjects(cytoplasmRegions[[1]],paintObjects(nucleusRegions[[1]],toRGB(b2_1h),col="#ff00ff"),col="#ff0000")
par(mfrow=c(1,1))
display(Im,method="raster")

# extract Features: here mean-intensity/cell + mean_variance/cell for each picture
means<-list()
var<-list()
for (i in 1:length(l1)){
  Features[[i]] <-data.frame(computeFeatures(cytoplasmRegions[[i]], l1[[i]], xname="cell", refnames="cytoplasm",methods.ref="computeFeatures.basic"))
  means[[i]]<-mean(Features[[i]]$cell.cytoplasm.b.mean)
  var[i]<-mean(Features[[i]]$cell.cytoplasm.b.sd^2)
}
means<-unlist(means)
var<-unlist(var)

# depending on aim, average means-per-well per condition
# dto. with variance, and convert to SE 
means_well = c(mean(means[1:3]),mean(means[4:6])) # example: combine B2-B4 and B5-B7
se<-c(sqrt(mean(c(var[1],var[2],var[3])))/sqrt(3),
           sqrt(mean(c(var[4],var[5],var[6])))/sqrt(3))
# plots the results, e.g., forest plot
means_perc<-(means_well/means_well[1])*100 # express as basal = 100%
se_perc<-(se/means_well[1])*100
pdf(...)
forest.default(means_perc,sei=se_perc,psize=1,
               slab=c("Basal","Insulin"),
               main="Glucose uptake - cell type, treatment...",
               xlab="percentage basal +- SEM",cex.lab=1,refline=100)
dev.off()


#############################################
############################################
