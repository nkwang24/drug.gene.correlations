library(reshape2)
library(parallel)

gene.lm<-function(drug,drugname,type,exp_matrix,transpose)
{
  #output=as.data.frame(matrix(ncol=7, nrow=GeneNumber))
  #colnames(output)=c("Gene","Drug","Intercept","Intercept Pval", "Predictor1 Pval", 
  #"Predictor1 Coefficient", "Predictor1 tval")
  #output[,1]=GeneName
  #exp_matrix=RNAi.drug[1,]
  
  skip=F
  if (transpose==F)
  {
    ##add gene name to the end of the row in gene matrix
    genename<-exp_matrix[1]
    exp_matrix=as.numeric(exp_matrix[2:length(exp_matrix)])
  }
  else
  {
    ##gene name is the top of the list 
    genename<-as.character(exp_matrix[1])
    gene<-as.numeric(exp_matrix[2:length(exp_matrix)])
    exp_matrix=melt(gene,id.vars = "gene")
  }
  
  merge=data.frame(drug,exp_matrix,type,stringsAsFactors = F)
  colnames(merge)=c("Drug", "Gene", "Type")
  
  tryCatch({summary(lm(Drug ~ Gene+Type, merge)) },
           
           error=function(e){
             
             print(e);
             
             skip<<-T;
             
           })
  
  if(skip==T)
  {
    output=data.frame(Gene=genename,
                      Drug=as.character(drugname), 
                      intercept=NA, 
                      InterceptPval=NA, 
                      Log10Pval=NA,
                      Coefficient=NA, 
                      tval=NA)
  }
  else{
    results= summary(lm(Drug ~ Gene+Type, merge))$coefficients
    output<-data.frame(Gene=genename,
                       Drug=as.character(drugname), 
                       intercept=results[1,1], 
                       InterceptPval=results[1,4], 
                       Log10Pval=-log10(results[2,4])*sign(results[2,3]),
                       Coefficient=results[2,1], 
                       tval=results[2,3])
    return(output)}
}

drug.lm<-function(drug,exp_matrix,type,transpose)
{
  drugvect=as.numeric(drug[2:length(drug)])
  drugvect=melt(drugvect,id.vars="drug")
  drugname=as.character(drug[1])
  
  if (transpose==T)
  {
    cl <- makeCluster(10)
    clusterExport(cl,varlist=c('gene.lm',"melt",ls()),envir=environment())
    output.vect=parApply(cl = cl, exp_matrix, MARGIN=1, 
                         FUN = function(x) gene.lm(drug=drugvect,
                                                   drugname=drugname,
                                                   type=type,
                                                   exp_matrix=x,
                                                   transpose=T))
    stopCluster(cl)
    out<-do.call(rbind,output.vect)
    write.table(out,paste0("/home/graeberlab/Desktop/RNAi_CTRP_10_28/rnai_",drugname,".txt"),
                sep="\t",row.names=F,col.names = T,quote=F)
  }
  else
  {
    cl <- makeCluster(10)
    clusterExport(cl,varlist=c('gene.lm',"melt",ls()),envir=environment())
    output.vect=parApply(cl = cl, exp_matrix, MARGIN=2, 
                         FUN = function(x) gene.lm(drug=drugvect,
                                                   drugname=drugname,
                                                   type=type,
                                                   exp_matrix=x,
                                                   transpose=F))
    stopCluster(cl)
    out<-do.call(rbind,output.vect)
    write.table(out,paste0("/home/graeberlab/Desktop/RNAi_CTRP_10_28/crispr_",drugname,".txt"),
                sep="\t",row.names=F,col.names = T,quote=F)
  }
  return(out)
}

run.2multipleregression<-function(drug_matrix, exp_matrix , type, transpose)
{
  if (transpose==T)
  {
    cl <- makeCluster(10)
    clusterExport(cl,varlist=c('drug.lm','gene.lm',"melt",
                               "makeCluster","clusterExport","parApply","stopCluster","do.call",
                               ls()),envir=environment())
    output.vect=parApply(cl = cl, drug_matrix, MARGIN=1, 
                         FUN = function(x) drug.lm(drug=x,
                                                   exp_matrix=exp_matrix,
                                                   type=type,
                                                   transpose=T))
    stopCluster(cl)
    out<-do.call(rbind,output.vect)
  }
  else
  {
    
    cl <- makeCluster(10)
    clusterExport(cl,varlist=c('drug.lm','gene.lm',"melt",
                               "makeCluster","clusterExport","parApply","stopCluster","do.call",ls()),envir=environment())
    output.vect=parApply(cl = cl, drug_matrix, MARGIN=1, 
                         FUN = function(x) drug.lm(drug=x,
                                                   exp_matrix=exp_matrix,
                                                   type=type,
                                                   transpose=F))
    stopCluster(cl)
    out<-do.call(rbind,output.vect)
  }
  return(out)
  #write.table(out,paste0("/home/graeberlab/Desktop/08-30-2019/",drugname,".txt"),
  #sep="\t",row.names=F,col.names = T,quote=F)
}

##Example
rnai.ctrp.drug<-run.2multipleregression(ctrp.rani, RNAi.ctrp, as.character(info.rnai.ctrp$lineage), transpose = T)
write.table(rnai.ctrp.drug,"/home/graeberlab/Desktop/rnai.ctrp.drug_ex.txt",sep="\t",quote=F, row.names=F, col.names = T)





