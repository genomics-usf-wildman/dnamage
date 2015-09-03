# Steve Horvath: Estimating DNAm age.
# This file assumes a data frame exists called data whose rows correspond to CpGs
# and whose first column reports the CpG identifier
# and whose remaining columns corresponds to samples (e.g. Illumina arrays).

library("impute")

##' Estimation of DNAm Age
##'
##' .. content for \details{} ..
##' @title 
##' @param data -- data frame whose rows correspond to CpGs and first
##' column is the CpG identifier
##' @param goldstandard -- name of the gold standard probe 21k data file
##' @param dat.clock -- AdditionalFile3.csv from horvath
##' @param trafo -- age transformation function.
##' @param anti.trafo -- inverse age transformation function
##' @param fastImputation 
##' @param meanXChromosome 
##' @return 
##' @author Steve Horvath
estimate.dnam.age <- function(data,goldstandard="probeAnnotation21kdatMethUsed.csv",dat.clock="AdditionalFile3.csv",trafo=function(x,adult.age=20){x <- (x+1)/(1+adult.age); ifelse(x<=1,log(x),x-1)},anti.trafo=function(x,adult.age=20){ifelse(x<0,(1+adult.age)*exp(x)-1,(1+adult.age)*x+adult.age)},fastImputation=FALSE,normalizeData=TRUE,alwaysImpute=FALSE,meanXchromosome=NULL) {
    ## if goldstandard is a character or a connection it must be a
    ## filename; read it.
    if (inherits(goldstandard,"connection") ||
        is.character(goldstandard)) {
        goldstandard <- read.csv(goldstandard)
    }
    if (inherits(dat.clock,"connection") ||
        is.character(dat.clock)) {
        dat.clock <- read.csv(dat.clock)
    }
    
    
    ##STEP 1: DEFINE QUALITY METRICS

    meanMethBySample <- as.numeric(apply(as.matrix(data[,-1]),2,mean,na.rm=TRUE))
    minMethBySample <- as.numeric(apply(as.matrix(data[,-1]),2,min,na.rm=TRUE))
    maxMethBySample <- as.numeric(apply(as.matrix(data[,-1]),2,max,na.rm=TRUE))
    
    datMethUsed <- t(data[,-1])
    colnames(datMethUsed) <- as.character(data[,1])

    noMissingPerSample <- apply(as.matrix(is.na(datMethUsed)),1,sum)
    table(noMissingPerSample)

    nSamples <- nrow(datMethUsed)

    ##STEP 2: Imputing
    if (alwaysImpute) {
        dimnames1=dimnames(datMethUsed)
        capture.output(datMethUsed <- data.frame(t(impute.knn(t(datMethUsed))$data)),
                       file="/dev/null")
        dimnames(datMethUsed)=dimnames1
    } else if (! fastImputation & nSamples>1 &
        max(noMissingPerSample,na.rm=TRUE)<3000 ){

        ## run the following code if there is at least one missing
        if ( max(noMissingPerSample,na.rm=TRUE)>0 ){
            dimnames1=dimnames(datMethUsed)
            datMethUsed= data.frame(t(impute.knn(t(datMethUsed))$data))
            dimnames(datMethUsed)=dimnames1
        } # end of if
    } # end of if (! fastImputation )

    if ( max(noMissingPerSample,na.rm=TRUE)>=3000 ) {
        fastImputation <- TRUE
    }
    if ( fastImputation | nSamples==1 ){
        noMissingPerSample <- apply(as.matrix(is.na(datMethUsed)),1,sum)
        table(noMissingPerSample)
        if ( max(noMissingPerSample,na.rm=TRUE)>0 &
            max(noMissingPerSample,na.rm=TRUE) >= 3000 ) {
            normalizeData=FALSE
        }
        
        ## run the following code if there is at least one missing
        if ( max(noMissingPerSample,na.rm=TRUE)>0 &
            max(noMissingPerSample,na.rm=TRUE) < 3000 ){
            dimnames1=dimnames(datMethUsed)
            for (i in which(noMissingPerSample>0) ){
                selectMissing1=is.na(datMethUsed[i,])
                datMethUsed[i,selectMissing1] <-
                    as.numeric(goldstandard$goldstandard2[selectMissing1])
            } # end of for loop
            dimnames(datMethUsed)=dimnames1
        } # end of if
    } # end of if (! fastImputation )

    ## STEP 3: Data normalization (each sample requires about 8
    ## seconds). It would be straightforward to parallelize this
    ## operation.

    if (normalizeData ){
        rownames(goldstandard) <- goldstandard$Name
        datMethUsedNormalized <-
            BMIQcalibration(datM=datMethUsed[,as.character(goldstandard$Name)],
                            goldstandard.beta=goldstandard$goldstandard2,
                            plots=FALSE)
    } else { 
       datMethUsedNormalized <- datMethUsed
    }
    rm(datMethUsed); gc()

    ## STEP 4: Predict age and create a data frame for the output (referred to as datout)
    selectCpGsClock=is.element(dimnames(datMethUsedNormalized)[[2]],
        as.character(dat.clock$CpGmarker[-1]))
    if ( sum( selectCpGsClock) < dim(dat.clock)[[1]]-1 ) {
        stop(paste0("The CpGs listed in column 1 of the input data did ",
                    "not contain the CpGs needed for calculating DNAm age. ",
                    "Make sure to input cg numbers such as cg00075967."))
    }
    if ( sum( selectCpGsClock) > dim(dat.clock)[[1]]-1 ) {
        stop(paste0("ERROR: The CpGs listed in column 1 of the input data ",
                    "contain duplicate CpGs. Each row should report only ",
                    "one unique CpG marker (cg number)."))}
    if (nSamples>1 ) {
        datMethClock0=data.frame(datMethUsedNormalized[,selectCpGsClock])
        datMethClock= data.frame(datMethClock0[ as.character(dat.clock$CpGmarker[-1])])
        dim(datMethClock)
        predictedAge=as.numeric(anti.trafo(dat.clock$CoefficientTraining[1]+as.matrix(datMethClock)%*% as.numeric(dat.clock$CoefficientTraining[-1])))
    } # end of if


    if (nSamples==1 ) {
        datMethUsedNormalized2=data.frame(rbind(datMethUsedNormalized,datMethUsedNormalized))
        datMethClock0=data.frame(datMethUsedNormalized2[,selectCpGsClock])
        datMethClock= data.frame(datMethClock0[ as.character(dat.clock$CpGmarker[-1])])
        dim(datMethClock)
        predictedAge=as.numeric(anti.trafo(dat.clock$CoefficientTraining[1]+as.matrix(datMethClock)%*% as.numeric(dat.clock$CoefficientTraining[-1])))
        predictedAge=predictedAge[1]
    } # end of if
    


    ## Let's add comments to the age prediction
    Comment=ifelse ( predictedAge <0, "Negative DNAm age.",
        ifelse ( predictedAge >100, "Old DNAm age.",
                rep("",length(predictedAge))))

    Comment[is.na(predictedAge)]="Age prediction was not possible. "


    if ( sum( selectCpGsClock) < dim(dat.clock)[[1]]-1 ) {
        Comment=rep(paste0("ERROR: The CpGs listed in column 1 of the input data ",
            "did not contain the CpGs needed for calculating DNAm age. ",
            "Make sure to input cg numbers such as cg00075967."),
            length(predictedAge) )
    }
    if ( sum( selectCpGsClock) > dim(dat.clock)[[1]]-1 ) {
        Comment=rep(paste0("ERROR: The CpGs listed in column 1 of the input data ",
            "contain duplicate CpGs. Each row should report only one ",
            "unique CpG marker (cg number).",
            length(predictedAge) ))
    }


    restSamples=-minMethBySample>0.05 | maxMethBySample>1.05;
    restSamples[is.na(restSamples)]=FALSE
    lab1="MAJOR WARNING: Probably you did not input beta values since either minMethBySample<-0.05 or maxMethBySample>1.05.";
    Comment[restSamples]= paste(Comment[restSamples],lab1)

    restSamples= noMissingPerSample >0 & noMissingPerSample <=100;lab1="WARNING: Some beta values were missing, see noMissingPerSample.";
    Comment[restSamples]= paste(Comment[restSamples],lab1)
    restSamples= noMissingPerSample >3000;
    lab1="MAJOR WARNING: More than 3k missing values!!";
    Comment[restSamples]= paste(Comment[restSamples],lab1)

    restSamples= noMissingPerSample >100 & noMissingPerSample <=3000 ;
    lab1="MAJOR WARNING: noMissingPerSample>100"
    Comment[restSamples]= paste(Comment[restSamples],lab1)
    restSamples=meanMethBySample>.35;
    restSamples[is.na(restSamples)]=FALSE
    lab1="Warning: meanMethBySample is >0.35";
    Comment[restSamples]= paste(Comment[restSamples],lab1)
    restSamples=meanMethBySample<.25;
    restSamples[is.na(restSamples)]=FALSE;
    lab1="Warning: meanMethBySample is <0.25"
    Comment[restSamples]= paste(Comment[restSamples],lab1)
    datout=data.frame(SampleID=colnames(data)[-1],
        DNAmAge=predictedAge,
        Comment,
        noMissingPerSample,
        meanMethBySample,
        minMethBySample,
        maxMethBySample)
    
    if ( !is.null( meanXchromosome) ){  

        if ( length( meanXchromosome)==dim(datout)[[1]] ){
            predictedGender=ifelse(meanXchromosome>.4,"female",
                ifelse(meanXchromosome<.38,"male","Unsure"))
            datout=data.frame(datout,predictedGender=predictedGender,meanXchromosome=meanXchromosome)

        } # end of if 

    } # end of if
    return(datout)
}



