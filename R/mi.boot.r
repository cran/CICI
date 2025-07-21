mi.boot <- function(x,fun,cond=NULL,pooled=FALSE,...){
  # checks
  if(any(unlist(lapply(x,class))!="gformula")){stop("Please provide a list of objects of class 'gformula'")}
  get.type <- function(xx){xx$setup$i.type}
  if(x[[1]]$setup$measure=="custom"){stop("Don't use 'mi.boot' on results produced by 'custom.measure'.\n Use mi.boot's options 'fun' and 'cond' to produce custom measures.\n")}
  if(is.list(x)==FALSE | (is.list(x)==TRUE & length(x)==1)){stop("Please provide a list (of length>1)")}
  i.type <- unlist(lapply(x,get.type))
  if(missing(fun)){stop("You have to supply 'fun'. \n If you are interested in the expectation, simply use 'fun=mean'.\n (and possibly pass on the option 'na.rm=T.')")}
  if(length(unique(i.type))>1){stop("You can only combine results of the same type")}else{i.type <- i.type[1]}
  if(is.null(x[[1]]$simulated.data)){stop("Set 'ret=TRUE' in gformula() to work with custom measures.")}
  # table
  mi.table <- x[[1]]$results
  # produce custom measure and calculate bootstrap s.e.
  if(x[[1]]$setup$B>0){x <- lapply(x,custom.measure,fun=fun,cond=cond,with.se=TRUE,verbose=FALSE,...)}else{
                       x <- lapply(x,custom.measure,fun=fun,cond=cond,with.se=FALSE,verbose=FALSE,...)
  }
  # psi
  getpsi <- function(xx){xx$results$psi}
  mi.psi <- apply(do.call("rbind",lapply(x,getpsi)),2,mean)
  mi.table$psi <- mi.psi
  # covariates if abar="natural"
  if(i.type=="natural" & is.null(x[[1]]$setup$Lnodes)==FALSE){
  rel.L.cols <-  colnames(x[[1]]$results)[!grepl(":",colnames(x[[1]]$results)) &  grepl("L_",colnames(x[[1]]$results))]
  if(length(rel.L.cols)>0){ 
  get.Ls <- function(xx){c(xx$results[,rel.L.cols])}
  mi.L <- apply(do.call("rbind",lapply(x,get.Ls)),2,mean)
  mi.table[,rel.L.cols] <- mi.L
  }else{x[[1]]$setup$Lnodes<-NULL;cat("natural course scenario for Lnodes not possible because values have not been provided\n")} # fix in future: with unequally distributed L's over time gformula struggles, and thus omit natural course for L's for now
  }
  # confidence intervals
  get.B <- function(xx){xx$setup$B}
  all.B <- unlist(lapply(x,get.B))
  if(all(all.B>1)){
  if(pooled==FALSE){
  getse <- function(xx){xx$results$se}
  mi.results <- mi.inference(est=lapply(x,getpsi),std.err=lapply(x,getse))
  mi.table$l95 <- mi.results$lower; mi.table$u95 <- mi.results$upper
      if(i.type=="natural" & is.null(x[[1]]$setup$Lnodes)==FALSE){
      getse.L <- function(xx){L.index <- grep(":se",colnames(x[[1]]$results)); c(xx$results[,L.index])  }
      mi.results.L <- mi.inference(est=lapply(x,get.Ls),std.err=lapply(x,getse.L))
      res.index <- grep(":l95",colnames(x[[1]]$results)); res.index2 <-  grep(":u95",colnames(x[[1]]$results))
      mi.table[,res.index] <- mi.results.L$lower
      mi.table[,res.index2]<- mi.results.L$upper
      }
  }else{
  get.bs <- function(xx){xx$b.results}
    if(i.type!="natural"){
    mi.results <- apply(do.call("rbind",lapply(x,get.bs)),2,quantile,probs=c(0.025,0.975))
    mi.table$l95 <- mi.results[1,]; mi.table$u95 <- mi.results[2,]
    }else{
    mi.results <- apply(do.call("rbind",lapply(x,get.bs)),2,quantile,probs=c(0.025,0.975))
    mi.table$l95[mi.table$a1=="natural"] <- mi.results[1,colnames(mi.results)%in%x[[1]]$setup$Ynodes]
    mi.table$u95[mi.table$a1=="natural"] <- mi.results[2,colnames(mi.results)%in%x[[1]]$setup$Ynodes]
    get.bs.obs  <- function(xx){xx$b.results2[[1]]}
    get.bs.diff <- function(xx){xx$b.results2[[2]]}
    mi.results2 <- apply(do.call("rbind",lapply(x,get.bs.obs)),2,quantile,probs=c(0.025,0.975))
    mi.results3 <- apply(do.call("rbind",lapply(x,get.bs.diff)),2,quantile,probs=c(0.025,0.975))
    mi.table$l95[mi.table$a1=="observed"] <- mi.results2[1,colnames(mi.results)%in%x[[1]]$setup$Ynodes]
    mi.table$u95[mi.table$a1=="observed"] <- mi.results2[2,colnames(mi.results)%in%x[[1]]$setup$Ynodes]
    mi.table$l95[mi.table$a1=="difference"] <- mi.results3[1,colnames(mi.results)%in%x[[1]]$setup$Ynodes]
    mi.table$u95[mi.table$a1=="difference"] <- mi.results3[2,colnames(mi.results)%in%x[[1]]$setup$Ynodes]
        if(is.null(x[[1]]$setup$Lnodes)==FALSE){
        res.index <- grep(":l95",colnames(x[[1]]$results)); res.index2 <-  grep(":u95",colnames(x[[1]]$results))
        mi.table[,c(res.index,res.index2)]<-NA
        cat("Still to implement CI's for L's \n")
        }
  }
  }
  }
  # weights & diagnostics
  if(any(grepl("crude_",colnames(x[[1]]$results)))){
  getweights <- function(xx){unlist(subset(xx$results,select=c("crude_weights","cond_weights")))}
  mi.table[,c("crude_weights","cond_weights")] <- apply(do.call("rbind",lapply(x,getweights)),2,mean)
  }
  mi.diagnostics<-NULL
  if(is.null(x[[1]]$diagnostics)==FALSE){
  getdiag1 <- function(xx){xx$diagnostics$crude_support}
  getdiag2 <- function(xx){xx$diagnostics$conditional_support}
  mi.crude <- x[[1]]$diagnostics$crude_support; mi.crude[1:nrow(mi.crude),1:ncol(mi.crude)] <- apply(do.call("rbind",lapply(lapply(x,getdiag1),unlist)),2,mean)
  mi.conditional <- x[[1]]$diagnostics$conditional_support; mi.conditional[1:nrow(mi.conditional),1:ncol(mi.conditional)] <- apply(do.call("rbind",lapply(lapply(x,getdiag2),unlist)),2,mean)
  mi.diagnostics <- list(
    crude_support = mi.crude,
    conditional_support = mi.conditional
  )
  }

  #return appropriate object of class gformula
  results <- list(results=mi.table,diagnostics=mi.diagnostics,setup=x[[1]]$setup)
  class(results) <- "gformula"
  results
}