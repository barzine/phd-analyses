report4chap4<-function(set="1",
                       prefix='transcriptComp',
                       quantMethod='htseq',
                       normMethod='fpkm',
                       proteinCodingOnly=TRUE,
                       removeMitoch=TRUE,
                       cutoff=0,
                       corrMethod='pearson',
                       log2scaling=TRUE,
                       pseudocount=1,
                       TissueCol,
                       datasetCol,
                       #remove0,
                       newDir=FALSE,
                       setWorkingDir=TRUE,
                       workingDir='phd-analyses/chapter4'
                       ){

  #Comparison between 5 or 2 studies: 5DF for castle, brawand, ibm, uhlen and gtex
  #                                   2DF for uhlen and gtex
  switch(as.character(set),
         "1"={compDF="5DF"
              labelDF1='castle'
              labelDF2='brawand'
              labelDF3='ibm'
              labelDF4='uhlen'
              labelDF5='gtex'
              nameDFs=c('Castle','Brawand','IBM','Uhlén','GTEx')
         },
         "2"={compDF="2DF"
              labelDF1='uhlen'
              labelDF2='gtex'
              nameDFs=c('Uhlén','GTEx')
              },
         "5"={compDF="5DF"
         labelDF1='castle'
         labelDF2='brawand'
         labelDF3='ibm'
         labelDF4='uhlen'
         labelDF5='gtex'
         nameDFs=c('Castle','Brawand','IBM','Uhlén','GTEx')
         })

  if(compDF=="5DF"){
    studyTag<-'five transcriptomic studies'
  }else{
    studyTag<-'two transcriptomic studies'
  }

  #Quantification method
  if(missing("quantMethod")) quantMethod<-'htseq'
  #if(!quantMethod %in% c("htseq","cufflinks")) quantMethod<-'htseq'

  #Normalisation method
  if(missing("normMethod")) normMethod<-'fpkm'

  #Keeping only protein coding genes?
  if(missing("proteinCodingOnly"))  proteinCodingOnly<-TRUE
  if(proteinCodingOnly){
    proteinTag<-"proteinCodingOnly"
  }else{
    proteinTag<-'anyGene'
  }

  #Keeping mitochondria?
  if(missing("removeMitoch")) removeMitoch<-TRUE
  if(removeMitoch){
    mitochTag<-'noMitoch'
  }else{
    mitochTag<-'withMitoch'
  }

  #threshold of expression
  if(missing(cutoff)) cutoff<-0

  #Correlation method
  if(missing("corrMethod")) corrMethod<-'pearson'

  #Scaling method
  if(missing("log2scaling")) log2scaling<-TRUE
  if(log2scaling){
    tagScaling<-'log2'
  }else{
    tagScaling<-'noScaling'
  }

  if(log2scaling){
    if(pseudocount<=0){print('Pseudocount need to be greater than 0, thus reverted to default: 1')
      pseudocount<-1}
  }

  #Keeping null values?
  #if(missing("remove0")) remove0<-FALSE

  #constants
  ENSEMBL_VER=76 # put as constant since only data mapped with Ensembl 76 is included in barzinePhdData package
  analysis_level="gene level"

  if(missing(datasetCol)) datasetCol<-rlang::eval_tidy(rlang::parse_expr("barzinePhdData::datasetCol"))
  if(missing(TissueCol)){
    TissueCol<-rlang::eval_tidy(rlang::parse_expr("barzinePhdData::TissueCol_62"))
    TissueCol["Fallopian tube"]<-TissueCol["Fallopiantube"]
    TissueCol["Frontal cortex"]<-TissueCol["Frontalcortex"]
    TissueCol["Skeletalmuscle"]<-TissueCol["Skeletal muscle"]
    TissueCol["Lymph node"]<-TissueCol["Lymph.node"]
    TissueCol["Spinal cord"]<-TissueCol["Spinal.cord"]
    TissueCol["Urinary.bladder"]<-TissueCol["Urinary bladder"]
    TissueCol["Smooth muscle"]<-TissueCol["Smooth.muscle"]
    TissueCol["Bone marrow"]<-TissueCol["Bone.marrow"]
  }

  if(!missing('newDir')) newDir<-FALSE

  NameReport<-gsub(' ','_',prefix)
  NameReport<-paste0(NameReport,compDF,'_ens',ENSEMBL_VER,'_',analysis_level,'_',quantMethod,'_',normMethod)
  NameReport<-paste0(NameReport,'_',proteinTag,'_',mitochTag,'_expCut',cutoff,'_',corrMethod,'_',tagScaling)
  NameReport<-gsub(' ','_',NameReport)

  dynTitle<-paste("Analysis for",studyTag)

  #more definition:
  loadDF<-function(labelDFdata,quantM=quantMethod,normM=normMethod){
    suppressPackageStartupMessages(library(barzinePhdData))
    if(exists(paste0(labelDFdata,'.',quantM,'.',normM,'.pooled'),where='package:barzinePhdData')){
      DFdata<-rlang::eval_tidy(rlang::parse_expr(paste0("barzinePhdData::",labelDFdata,'.',quantM,'.',normM,'.pooled')))
    }else{
      DFdata<-rlang::eval_tidy(rlang::parse_expr(paste0("barzinePhdData::",labelDFdata,'.',quantM,'.',normM)))
    }
    return(DFdata)
  }

  protein.coding<-rlang::eval_tidy(rlang::parse_expr("barzinePhdData::g.Pcoding"))
  mitoch<-rlang::eval_tidy(rlang::parse_expr("barzinePhdData::g.mitoch[,1]"))
  if(compDF=="5DF"){
      TigerTSselect<-rlang::eval_tidy(rlang::parse_expr("barzinePhdData::Tiger4Tissues"))
  }else{
     TigerTSselect<-rlang::eval_tidy(rlang::parse_expr("barzinePhdData::Tiger23Tissues"))
  }

  if(setWorkingDir) setwd(workingDir)
  if(newDir){
    rmarkdown::render('templateChap4.Rmd',output_file=NameReport,output_dir=NameReport)
  }else{
    rmarkdown::render('templateChap4.Rmd',output_file=NameReport)
  }
}
