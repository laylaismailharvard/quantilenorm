##### Load Libraries #####
library(NanoStringNCTools)
library(GeomxTools)
library(EnvStats)
library(ggiraph)
library(openxlsx)
library(edgeR)
library(DT)
library(lme4)
library(ggplot2)
library(ggridges)
library(matrixStats)
library(readxl)
library(sqldf)
library(BiocManager)
library(limma)
library(pheatmap)

##### Set Directory & Annotation Sheet #####
pathdir <- "/Users/laylaismail/Desktop/summer 25/hearttransplant"
annotationSheet <- "LW"  # Worksheet in annotation file

##### Data Loading Function (based on working BM code) #####
DataLoading <- function(pathdir, sheet) {
  # Step 1: List DCC and PKC Files
  DCCFiles <- dir(file.path(pathdir, "dcc"), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
  PKCFiles <- dir(file.path(pathdir, "pkc"), pattern = ".pkc$", full.names = TRUE, recursive = TRUE)
  
  # Step 2: Locate Sample Annotation File (more flexible search)
  SampleAnnotationFile <- list.files(file.path(pathdir, "annotationkidney"), 
                                     full.names = TRUE, recursive = TRUE, 
                                     pattern = "^[^~].*\\.xlsx$")
  
  # Take first xlsx file found
  if (length(SampleAnnotationFile) == 0) {
    stop("No annotation Excel file found in annotation directory")
  }
  SampleAnnotationFile <- SampleAnnotationFile[1]
  
  cat("Using annotation file:", SampleAnnotationFile, "\n")
  
  # Step 3: Read GeoMx Data
  GeoMxData <- readNanoStringGeoMxSet(
    dccFiles = DCCFiles,
    pkcFiles = PKCFiles,
    phenoDataFile = SampleAnnotationFile,
    phenoDataSheet = sheet,
    phenoDataDccColName = "Sample_ID",
    protocolDataColNames = c("Aoi", "Roi"),
    experimentDataColNames = c("Panel")
  )
  
  # Step 4: Return GeoMxData
  return(GeoMxData)
}

dat<- DataLoading(pathdir, annotationSheet)

dat <- shiftCountsOne(dat, useDALogic = TRUE)
dat <- setSegmentQCFlags(dat)
dat <- setBioProbeQCFlags(dat)
badprobes <- apply(fData(dat)$QCFlags[,1:2], 1, any)
dat <- dat[!badprobes,]

## function to subtract out the geometric mean of the negative controls first
bgFun <- function(obj){
  nc <- negativeControlSubset(obj)
  bg <- colMeans(log2(exprs(nc)))
  vals <- log2(exprs(obj))
  vals <- sweep(vals, 2, bg)
  assayDataElement(obj, "exprs", validate = TRUE) <- vals
  obj
}

genedat <- aggregateCounts(dat)[,-1]

## LOQ based on recommendations from nanostring
LOQ <- data.frame(row.names = colnames(genedat))
LOQ[,"HS_R_NGS_WTA_v1.0"] <- pmax(2, pData(genedat)$NegGeoMean_Hs_R_NGS_WTA_v1.0 * pData(genedat)$NegGeoSD_Hs_R_NGS_WTA_v1.0^2)
pData(genedat)$LOQ <- LOQ

##  We can now make a matrix that is the same dimension as our data, 
## but Boolean, to indicate if a given gene for a given sample is above or below the LOQ 
## (this is the first part of the addLOQ function).

loqmat <- t(t(exprs(genedat)) > pData(genedat)$LOQ[,1])

## TARGET FILTER option from GeoMx: filter on genes where at least 50/191 = 25% of samples are above LOQ
genedat <- genedat[rowSums(loqmat) > 50,]

samps <- cbind(pData(protocolData(dat))[-1,"Roi"],
               pData(dat)[-1, c("Patient", "Structure", "TCMR_Grade", "ABMR", "Quilty"), drop = FALSE])

colnames(samps)[1] <- "ROI_ID"

## ADD LOQ stuff
addLOQ <- function(geomx){
  loqmat <- t(t(exprs(geomx)) > pData(geomx)$LOQ[,1])
  pData(geomx)$GenesDetected <- colSums(loqmat)
  pData(geomx)$PctGenesDetected <- colSums(loqmat)/nrow(geomx)
  fData(geomx)$SamplesDetected <- rowSums(loqmat)
  fData(geomx)$PctSamplesDetected <- rowSums(loqmat)/ncol(geomx)
  geomx
}

genedat <- addLOQ(genedat)

########################### Quantile Normalization function##############################
myNorm4 <- function(obj, fromElt = "exprs", ...) {
  vals <- assayDataElement(obj, fromElt)
  vals <- normalizeBetweenArrays(log2(vals), "quantile")
  assayDataElement(obj, "exprs_norm", validate = TRUE) <- 2^vals
  obj
}

# Write Quantile normalized values to file
normed_data <- myNorm4(genedat)
ramdat <- (assayDataElement(normed_data, "exprs_norm"))

colnames(ramdat) <- with(samps, paste(Patient, Structure, TCMR_Grade, ABMR, Quilty, ROI_ID, sep = "_"))

#save quantile normalized data for kidney
write.csv(ramdat, "kidney.quantilenorm.data.csv", row.names = TRUE)