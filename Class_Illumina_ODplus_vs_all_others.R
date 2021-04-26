# Creating a model to distinguish patients with sepsis3 from all other patients
#
# Input data are expression values from Illumina microarray and 
#  metadata from patients (which day before diagnosis of sepsis, infection yes/no, organ dysfunction yes/no, ...
#
# Data are available upon request and from the GEO dataset GSE159750.
#
# author: Wolfgang Schmidt-Heck

# Pre-selection of relevant genes with BORUTA 
#
# With the relevant features found, a model is learned with 5-Fold cross validation
#   which is repeated 25 times.
# The importance of the features for the model are calculated. The feature 
#  with the smallest importance is removed and the model is learned with the reduced number of features.
#
# Output is a matrix of performance parameters
#colnames(oMetSum)
# [1] "AUC"                  "ACC"                  "C11"                  "C21"                 
# [5] "C12"                  "C22"                  "AnzFeat"              "Sensitivity"         
# [9] "Specificity"          "Pos Pred Value"       "Neg Pred Value"       "Precision"           
# [13] "Recall"               "F1"                   "Prevalence"           "Detection Rate"      
# [17] "Detection Prevalence" "Balanced Accuracy"   

library(caret)
library(Boruta)
require(ROCR)
library(randomForest)

# Reading in the expression matrix. 
# The rows correspond to the genes (Illumina IDs) and the columns correspond to the patients (sample_name).
yaveE <- read.table("Expression_illumina_Chips.txt", sep="\t")


# Reading in the metadata for patients
pheno_Illu <- read.table("pheno_illumina.txt", sep="\t", header = TRUE)

#    lfdNr patient Typ sample_name            group              Day    MetaYes kzID  classifier DysFunction1P DysFunction2P DysFunctionAllP
##     1   U1151   S    U1151A00           Pre-op          PreSurgery      1    SPre     S001             N             N          N
##     2   U1151   S    U1151A01           Day -1 Sepsis          1        1    S1       S001             N             N          N
##     3   U1390   C    U1390A00           Pre-op          PreSurgery      1    CPre     C001          <NA>          <NA>          N    
##     4   U1390   C    U1390A01           Day -1 Baseliner       1        1    C1       C001          <NA>          <NA>          N
##     5   J1018   S    J1018A00           Pre-op          PreSurgery      1    SPre     S003             Y             Y          N
##     6   J1018   S    J1018A01           Day -1 Sepsis          1        1    S1       S003             Y             Y          N  



# Reading in the assignment of official gene symbols to the Illumina IDs
illu_sym <- read.table("Illu_ID_Symbol_anno.txt.csv", sep="\t", header = TRUE, quote = "\"'")

#    IL.ID         Symbol
## 1 ILMN_1766896  ZDHHC19
## 2 ILMN_2116877  OLFM4
## 3 ILMN_1776519  RAP1GAP
## 4 ILMN_1674574  VNN1
## 5 ILMN_1711030  OPLAH
## 6 ILMN_3238814  None

# Patients are selected for classification on the basis of this attribute.
tc <- pheno_Illu$kzID

mat <- yaveE
matN <- 0*mat

# Preselection of genes; select Genes with fold-change greater 1.3
comL <- read.table("PreSel_IlluIDs_F13_uni_anno_and_pl-allOther.txt.csv", header=TRUE, sep="\t")  # Fold 1.3
#   IL.ID         Symbol
##1 ILMN_1805104     ABAT
##2 ILMN_1704579   ABCA13
##3 ILMN_1747627    ABCA2
##4 ILMN_3291778   ABCB10
##5 ILMN_3247608 ABCB10P4
##6 ILMN_2343048    ABCB9

glist <- row.names(yaveE)

cl <-as.character(comL$IL.ID)

ig <- intersect(glist, cl)

ix <-match(glist, cl)

idx <- which(!is.na(ix))


mat.s <- mat[idx, ]

aucA <- list()
selFA <- list()

# Organdysfunction Yes/No
tCDys <- pheno_Illu$DysFunction1P

# Was the patient included in the study?
meta <- pheno_Illu$MetaYes

X <- t(mat.s)


# Renaming of the genes in the expression matrix. The new name is the gene symbol.
cN <-colnames(X)

cNe <- cN
for (i in 1:length(cN)) {
  id <- which(!is.na(match(illu_sym$IL.ID,cN[i])))
  cNe[i]<-as.character( illu_sym$Symbol[id])
}
colnames(X) <- cNe   

# Patient classifier used. Classifiers starting with "S" are "sepsis patients" and those starting with "C" are the control group.
tc <- pheno_Illu$kzID

tcDay <- pheno_Illu$Day

# Creating a classifier for the discrimination "organdysfunction plus vs all other patients"  for the respective days before diagnosis

#Day1 (One day before diagnosis.)
id1 <- which(tcDay == 1 & meta == 1) 
cb1 <- tCDys[id1]
id <- which(is.na(cb1))

cb1[id] <- "N"
oMati1 <- X[id1, ]

# Day2 (Two day before diagnosis.)
id2 <- which(tcDay == 2 & meta == 1)
cb2 <- tCDys[id2]
id <- which(is.na(cb2))

cb2[id] <- "N"
oMati2 <- X[id2, ]


# Day3
id3 <- which(tcDay == 3 & meta == 1) 
cb3 <- tCDys[id3]
id <- which(is.na(cb3))

cb3[id] <- "N"
oMati3 <- X[id3, ]

# Labels for the naming of result tables
kz <- c("illumina_Day-1_DysP1_vs_all-others", "illumina_Day-2_DysP1_vs_all-others", "illumina_Day-3_DysP1_vs_others")

# Loop for calculating the classification models for the respective days (day -1, -2 and -3 before diagnosis)
require(e1071)
# Loop for calculating the classification models for the respective days (day -1, -2 and -3 before diagnosis)
for (ikk in 1:3) {
  
  # matrix of predictors -> xM
  # response vector - yM
  if (ikk == 1) {
    xM <- oMati1 
    yM <- as.factor(cb1)
  } else if (ikk == 2) {
    xM <- oMati2 
    yM <- as.factor(cb2)
  } else if (ikk == 3) {
    xM <- oMati3 
    yM <- as.factor(cb3)
    
  }
  clM <- matrix(0, nrow = length(yM), ncol = length(row.names(xM)))
  row.names(clM) <- row.names(xM)
  clMN <- row.names(clM)
  
  
  #Search for relevant features using BORUTA with 5-fold cross-validation repeated 25 times.
  fac <- yM
  nFold = 5
  nTimes = 25
  
  indexFG <- createMultiFolds(as.factor(fac), k = nFold, times = nTimes)
  
  # Selection of important features
  fName <- colnames(xM)
  anzR <- 1
  sMat <- matrix(0, nrow = length(fName), ncol = nFold*nTimes*anzR)
  
  
    for (j in 1:anzR) {
      for (i in 1:(nFold*nTimes)) {
        
        idx <- indexFG[[i]]
        
        b_out <- Boruta(xM[indexFG[[i]], ], yM[indexFG[[i]]] ,maxRuns = 1000, doTrace=0) 
        print(b_out)
        #  plot(b_out)
        
        idB1 <- which(b_out$finalDecision == "Confirmed" | b_out$finalDecision == "Tentative" )
        sMat[idB1, i] <- 1
      }
    }
    
  row.names(sMat) <- fName
  
  
  
  rs <- rowSums(sMat) /(nFold*nTimes/100) # Calculate in per cent
  
  rss <- sort(rs, decreasing = T)
  head(rss)
  rssi <- order(rs, decreasing = T)
  
  rss <- sort(rs)
  
  # Use of all features found to be relevant in at least one model.
  xMF <- xM[, rs>0] 
 # xMF <- xM ## test all genes
  #xMF <- xM[, 1:20] ## test all genes
  
  
  xMF_imp <- xM[, rssi[1:20]]
  
  write.table(rbind(rs[rs>0], xM[,rs>0]), file =paste0("ImpFeat_",kz[ikk], ".xls"), sep="\t")
  
  
  
  # Parameter Tuning for RandomForest (https://rpubs.com/phamdinhkhanh/389752)
  bT <- 0
  
  if (bT) {
    
    # Tuning RF
    metric<-'RMSE'
    
    customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
    
    customRF$parameters <- data.frame(parameter = c("maxnodes", "ntree"), class = rep("numeric", 2), label = c("maxnodes", "ntree"))
    
    customRF$grid <- function(x, y, len = NULL, search = "grid") {}
    
    customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
      randomForest(x, y, maxnodes = param$maxnodes, ntree=param$ntree, ...)
    }
    
    customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
      predict(modelFit, newdata)
    customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
      predict(modelFit, newdata, type = "prob")
    customRF$sort <- function(x) x[order(x[,1]),]
    customRF$levels <- function(x) x$classes
    
    # Set grid search parameters
    control <- trainControl(method="repeatedcv", number=10, repeats=3, search='grid')
    
    # Outline the grid of parameters
    tunegrid <- expand.grid(.maxnodes=seq(from=5, to=50, by=5), .ntree=seq(from=50, to=500, by = 50))
    #set.seed(seed)
    
    # Train the model
    rf_gridsearch <- train(x=xMF, y=yM, method=customRF, metric=metric, tuneGrid=tunegrid, trControl=control)
    
    plot(rf_gridsearch)
    
    # Best parameters:
    rf_gridsearch$bestTune
    #maxnodes ntree
    #51       10    150
  } # bT
  mxNodes <- 10
  nuTree <- 150  
  
  
  # for optimized parameter
  # Train the model 
  yM <- as.factor(yM)
  require(numbers)
  require(Metrics)
  require(ROCR)
  rowN <-  colnames(xMF)
  
  lVM <- length(rowN)
  impM <- matrix(0, nrow = lVM, ncol = length(indexFG))
  
  #  rowN <- row.names(iVar[[1]])
  row.names(impM) <- rowN
  oMet <- matrix(0, nrow=(length(indexFG)), ncol=18)
  colnames(oMet) <- c("AUC", "ACC", "C11", "C21", "C12", "C22", "AnzFeat", "Sensitivity" ,         "Specificity",          "Pos Pred Value"    ,   "Neg Pred Value"  ,     "Precision"   ,        
                      "Recall"   ,            "F1"     ,              "Prevalence"     ,      "Detection Rate" ,      "Detection Prevalence",
                      "Balanced Accuracy"  )
  
  oMetSum <- matrix(0, nrow=(lVM), ncol=18)
  
  colnames(oMetSum) <- c("AUC", "ACC", "C11", "C21", "C12", "C22", "AnzFeat", "Sensitivity" ,         "Specificity",          "Pos Pred Value"    ,   "Neg Pred Value"  ,     "Precision"   ,        
                         "Recall"   ,            "F1"     ,              "Prevalence"     ,      "Detection Rate" ,      "Detection Prevalence",
                         "Balanced Accuracy"  )
  
  
  for (icc in 1:(lVM-2)) {
    
    
    idI <- which(!is.na(match(colnames(xMF), rowN)))
    oMet <- matrix(0, nrow=(length(indexFG)), ncol=18)
    colnames(oMet) <- c("AUC", "ACC", "C11", "C21", "C12", "C22", "AnzFeat", "Sensitivity" ,         "Specificity",          "Pos Pred Value"    ,   "Neg Pred Value"  ,     "Precision"   ,        
                        "Recall"   ,            "F1"     ,              "Prevalence"     ,      "Detection Rate" ,      "Detection Prevalence",
                        "Balanced Accuracy"  )
    
    
    
    iz <- 0
    iVar <- list()
    
    for (i in 1:(length(indexFG))) {
      if (mod(i, nFold) == 1) { iz <- iz + 1}
      
      idx <- indexFG[[i]] 
      fac <- yM
      
      lxMF <- xMF[idx,idI]  # learn data
      txMF <- xMF[-idx,idI] # test data
      lyM  <- yM[idx]
      tyM  <-  as.factor(yM[-idx])
      
      tidx <- seq(from=1, to=length(yM), by=1)[-idx]
      regr <- randomForest(x=lxMF, y=lyM , replace= T, maxnodes = mxNodes, ntree = nuTree, importance=T)
      tmp <- predict(regr, newdata = txMF , type = "prob")
      red <- prediction((tmp[,2]), (tyM))  # only for two class problems
      perf <- performance(red, measure = "tpr", x.measure = "fpr")
      #  plot(perf)
      
      perf1 <- performance(red, measure = "auc")
      perf2 <- ROCR::performance(red, measure = "tpr")
      
      
      tm <- predict(regr, newdata = txMF)
      
      
      cM<- confusionMatrix(tm, tyM)
      
      oMet [(i), 1] <- as.numeric(perf1@y.values)
      oMet [(i), 2] <- cM$overall[1]
      oMet [(i), 3:6] <- cM$table
      if (as.numeric(perf1@y.values) != 0.5) {
        iVar[[i]]<- randomForest::importance(regr, type=2)
      }   
      oMet[i, 7] <- length(idI)
      oMet[i, 8:18] <- cM$byClass
      
      
      
    } # end:  for (i in 1:(length(indexFG)))
    
    
    
    impM <- matrix(0, nrow = length(rowN), ncol = ((length(indexFG))))
    #rowN <- feName
    row.names(impM) <- rowN
    
    #  impM[,1]<-unlist(iVar[[1]])
    
    for (i in 1:length(iVar)){
      if (!is.null(iVar[[i]])) {
        aRN <- rownames(iVar[[i]])
        idU <- which(!is.na(match(rowN, aRN)))
        
        impM[idU,i] <- unlist(iVar[[i]])#[idU]
      }   
    }
    
    rSum <- rowSums(impM)
    #    write.table(cbind(rSum, impM), file=paste0("Var_Use_", kz[ikk],"_IterStep ", as.character(icc), ".xls"), sep="\t")
    
    oMetSum[(icc),] <- colMeans(oMet, na.rm=T)
    
    rS1 <- rowSums(impM)
    
    idM <- which.min(rSum)
    write.table(rSum, file=paste0("Var_Use_it_",  kz[ikk], "_IterStep ", as.character(icc),  ".xls"), sep="\t")
    rowN<- rowN[-idM]
    
  } # icc
  
  
  
  write.table(oMetSum, file=paste0("Summary_Metrix__", kz[ikk], ".xls"), sep ="\t")
} #ikk

