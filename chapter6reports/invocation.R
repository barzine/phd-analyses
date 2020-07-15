datasetCol<-rlang::eval_tidy(rlang::parse_expr("barzinePhdData::datasetCol"))
TissueCol<-rlang::eval_tidy(rlang::parse_expr("barzinePhdData::TissueCol_62"))
TissueCol["Fallopian tube"]<-TissueCol["Fallopiantube"]
TissueCol["Frontal cortex"]<-TissueCol["Frontalcortex"]
TissueCol["Skeletalmuscle"]<-TissueCol["Skeletal muscle"]
TissueCol["Lymph node"]<-TissueCol["Lymph.node"]
TissueCol["Spinal cord"]<-TissueCol["Spinal.cord"]
TissueCol["Urinary.bladder"]<-TissueCol["Urinary bladder"]
TissueCol["Smooth muscle"]<-TissueCol["Smooth.muscle"]
TissueCol["Bone marrow"]<-TissueCol["Bone.marrow"]

report4chap6_tissue<-function(set="1",
                       prefix='protRNAcomp_',
                       quantMethod='htseq',
                       normMethod='fpkm',
                       protQuant='ppkm',
                       proteinCodingOnly=TRUE,
                       removeMitoch=TRUE,
                       cutoffprot=0,
                       cutoff=0,
                       corrMethod='pearson',
                       log2scaling=TRUE,
                       pseudocount=1,
                       TissueCol=TissueCol,
                       datasetCol=datasetCol,
                       remove0=FALSE,
                       newDir=FALSE,
                       setWorkingDir=TRUE,
                       workingDir='/hps/nobackup2/ma/mitra/phd-analysis/chapter6'
){

  #Comparison between 3 or 2 studies: 3DF for pandey, uhlen and gtex,
  #                                   2DF for pandey and gtex,
  switch(as.character(set),
         "1"={compDF="3DF"
         labelDF1='pandey'
         labelDF2='uhlen'
         labelDF3='gtex'
         nameDFs=c('Pandey','Uhlén','GTEx')
         },
         "2"={compDF="2DF"
         labelDF1='pandey'
         labelDF2='uhlen'
         nameDFs=c('Pandey','Uhlén')
         },
         "3"={compDF="3DF"
         labelDF1='pandey'
         labelDF2='uhlen'
         labelDF3='gtex'
         nameDFs=c('Pandey','Uhlén','GTEx')
         })


  if(compDF=="3DF"){
    studyTag<-'three datasets, twelve tissues studies'
  }else{
    studyTag<-'two datasets, fifteen tissues studies'
  }

  #Quantification method
  if(missing("quantMethod")) quantMethod<-'htseq'
  #if(!quantMethod %in% c("htseq","cufflinks")) quantMethod<-'htseq'

  if(missing(protQuant)) protQuant='ppkm'
  protQuant<-tolower(protQuant)
  if(!protQuant%in%c('top3','ppkm')) protQuant<-'ppkm'

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
  if(missing(cutoffprot)) cutoffprot<-0

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
  if(missing("remove0")) remove0<-FALSE

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

  preselected12T3DF<-c("Adrenal", "Colon", "Oesophagus", "Heart", "Kidney", "Liver",
                       "Lung", "Ovary", "Pancreas", "Prostate", "Testis", "Urinarybladder")
  preselected15T2DF<-c("Adrenal", "Colon", "Oesophagus", "Gallbladder", "Heart", "Kidney",
                       "Liver", "Lung", "Ovary", "Pancreas", "Placenta", "Prostate",
                       "Rectum", "Testis", "Urinarybladder")


  NameReport<-gsub(' ','_',prefix)
  NameReport<-paste0(NameReport,compDF,'_ens',ENSEMBL_VER,'_',analysis_level,'_rna_',quantMethod,'_',normMethod)
  NameReport<-paste0(NameReport,'_expRnaCut',cutoff,"_")
  NameReport<-paste0(NameReport,"_prot",protQuant,'_expProtCut',cutoffprot,"_")
  NameReport<-paste0(NameReport,'_',proteinTag,'_',mitochTag,'_',corrMethod,'_',tagScaling)
  if(log2scaling) NameReport<-paste0(NameReport,"_pseud",pseudocount)

  if(remove0){
    NameReport<-paste0(NameReport,'noNullval')
  }else{
    NameReport<-paste0(NameReport,'withNullval')
  }
  NameReport<-gsub(' ','_',NameReport)

  dynTitle<-paste("Analysis for",studyTag)

  #more definition:
  loadRnaDF<-function(labelDFdata,quantM=quantMethod,normM=normMethod){
    suppressPackageStartupMessages(library(barzinePhdData))
    if(exists(paste0(labelDFdata,'.',quantM,'.',normM,'.pooled'),where='package:barzinePhdData'))
      DFdata<-rlang::eval_tidy(rlang::parse_expr(paste0("barzinePhdData::",labelDFdata,'.',quantM,'.',normM,'.pooled')))
    return(DFdata)
  }


  loadprotDF<-function(labelDFdata,quantM=protQuant,cDF=compDF){
    suppressPackageStartupMessages(library(barzinePhdData))
    if(exists(paste0(labelDFdata,'.',quantM),where='package:barzinePhdData'))
      DFdata<-rlang::eval_tidy(rlang::parse_expr(paste0("barzinePhdData::",labelDFdata,'.',quantM)))
    return(DFdata)
  }

  protein.coding<-rlang::eval_tidy(rlang::parse_expr("barzinePhdData::g.Pcoding"))
  if(analysis_level=="gene level"){
    mitoch<-rlang::eval_tidy(rlang::parse_expr("barzinePhdData::g.mitoch[,1]"))
  }else{
    mitoch<-rlang::eval_tidy(rlang::parse_expr("barzinePhdData::g.mitoch[,2]"))
  }

  if(setWorkingDir) setwd(workingDir)

  print(paste("... generating",NameReport))
    if(newDir){
      rmarkdown::render('templateChap6-tissues.Rmd',output_file=paste0(NameReport,'_tissue'),output_dir=paste("tissues/",NameReport))

    }else{
      rmarkdown::render('templateChap6-tissues.Rmd',output_file=paste0(NameReport,'_tissue'))

    }
}




report4chap6_genes<-function(set="1",
                              prefix='protRNAcomp_genes_',
                              quantMethod='htseq',
                              normMethod='fpkm',
                              protQuant='ppkm',
                              corrMethod='pearson',
                              log2scaling=TRUE,
                              pseudocount=1,
                              TissueCol=TissueCol,
                              datasetCol=datasetCol,
                              newDir=FALSE,
                              repTimes=10000,
                              nJOB=400,
                              setWorkingDir=TRUE,
                              workingDir='~/Documents/workspace/phd-analyses/chapter6'
){

  #Comparison between 3 or 2 studies: 3DF for pandey, uhlen and gtex,
  #                                   2DF for pandey and gtex,
  switch(as.character(set),
         "1"={compDF="3DF"
         labelDF1='pandey'
         labelDF2='uhlen'
         labelDF3='gtex'
         nameDFs=c('Pandey','Uhlén','GTEx')
         },
         "2"={compDF="2DF"
         labelDF1='pandey'
         labelDF2='uhlen'
         nameDFs=c('Pandey','Uhlén')
         },
         "3"={compDF="3DF"
         labelDF1='pandey'
         labelDF2='uhlen'
         labelDF3='gtex'
         nameDFs=c('Pandey','Uhlén','GTEx')
         })


  if(compDF=="3DF"){
    studyTag<-'three datasets, twelve tissues studies'
  }else{
    studyTag<-'two datasets, fifteen tissues studies'
  }

  #Quantification method
  if(missing("quantMethod")) quantMethod<-'htseq'
  #if(!quantMethod %in% c("htseq","cufflinks")) quantMethod<-'htseq'

  if(missing(protQuant)) protQuant='ppkm'
  protQuant<-tolower(protQuant)
  if(!protQuant%in%c('top3','ppkm')) protQuant<-'ppkm'

  #Normalisation method
  if(missing("normMethod")) normMethod<-'fpkm'


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

  preselected12T3DF<-c("Adrenal", "Colon", "Oesophagus", "Heart", "Kidney", "Liver",
                       "Lung", "Ovary", "Pancreas", "Prostate", "Testis", "Urinarybladder")
  preselected15T2DF<-c("Adrenal", "Colon", "Oesophagus", "Gallbladder", "Heart", "Kidney",
                       "Liver", "Lung", "Ovary", "Pancreas", "Placenta", "Prostate",
                       "Rectum", "Testis", "Urinarybladder")


  NameReport<-gsub(' ','_',prefix)
  NameReport<-paste0(NameReport,compDF,'_ens',ENSEMBL_VER,'_',analysis_level,'_rna_',quantMethod,'_',normMethod)
  NameReport<-paste0(NameReport,"_prot",protQuant)
  NameReport<-paste0(NameReport,'_',corrMethod,'_',tagScaling)
  if(log2scaling) NameReport<-paste0(NameReport,"_pseud",pseudocount)

  NameReport<-gsub(' ','_',NameReport)

  dynTitle<-paste("Analysis for",studyTag)

  #more definition:
  loadRnaDF<-function(labelDFdata,quantM=quantMethod,normM=normMethod){
    suppressPackageStartupMessages(library(barzinePhdData))
    if(exists(paste0(labelDFdata,'.',quantM,'.',normM,'.pooled'),where='package:barzinePhdData'))
      DFdata<-rlang::eval_tidy(rlang::parse_expr(paste0("barzinePhdData::",labelDFdata,'.',quantM,'.',normM,'.pooled')))
    return(DFdata)
  }

  loadprotDF<-function(labelDFdata,quantM=protQuant,cDF=compDF){
    suppressPackageStartupMessages(library(barzinePhdData))
    if(exists(paste0(labelDFdata,'.',quantM),where='package:barzinePhdData'))
      DFdata<-rlang::eval_tidy(rlang::parse_expr(paste0("barzinePhdData::",labelDFdata,'.',quantM)))
    return(DFdata)
  }


  if(setWorkingDir) setwd(workingDir)

  #print(paste("... generating",NameReport))
  if(newDir){
    rmarkdown::render('chap6Gene.Rmd',output_file=paste0(NameReport,'_gene'),output_dir=NameReport)

  }else{
    rmarkdown::render('chap6Gene.Rmd',output_file=paste0(NameReport,'_gene'))

  }
}

