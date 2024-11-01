custom.measure <- function(X,fun=NULL,cond=NULL,verbose=TRUE,with.se=FALSE,...){
if(is.null(fun)){stop("'fun' need to be specified")}
if(is.function(fun)==FALSE){stop("'fun' needs to be a function")}
if(is.null(X$simulated.data)){stop("Set 'ret=TRUE' in gformula() to work with custom measures.")}
if(verbose==TRUE){if(X$setup$i.type=="custom"){cat("Note: you chose custom interventions. \n 'custom.measure' asssumes that the intervention values at the first time point (time=1) are not identical.\n\n")}}
if(X$setup$measure=="custom"){stop("You can't apply 'custom measure' on an X onto which you already applied 'custom measure'.")}
if(X$setup$B==0 & with.se==T){stop("Standard error can only be calculated if B>0")}
if(X$setup$B>0){if(X$setup$B!=(length(X$simulated.data)-1)){X$setup$B<- (length(X$simulated.data)-1) }}
sdata <- X$simulated.data
odata <- X$observed.data
results.table <- X$results
if(length(grep("mult",X$setup$fams))>0 & X$setup$i.type=="natural"){
if(is.null(sum(sapply(odata[[1]],is.factor))>0)==FALSE){if(sum(sapply(odata[[1]],is.factor))>0){
for(b in 1:(X$setup$B+1)){odata[[b]][,sapply(odata[[b]],is.factor)] <- apply(sapply(( subset(odata[[b]],select=sapply(odata[[b]],is.factor))),as.character),2,as.numeric)}
}}
}
if(X$setup$i.type=="individual"){
  custom <- any(unlist(lapply(lapply(apply(X$results[, grep("^a", names(X$results)), drop=F],1,unique),na.omit),length))>1 )
  if(custom==TRUE){X$setup$i.type<-"custom"}else{X$setup$i.type<-"standard"}
  if(custom==TRUE){X$setup$abar <- do.call(rbind,lapply(lapply(X$setup$abar,colnames),as.numeric))}else{
                   X$setup$abar <- do.call(rbind,lapply(lapply(X$setup$abar,colnames),as.numeric))[,1]
  }
}

if(X$setup$i.type=="standard"){
for(i in 1:length(sdata[[1]])){
if(is.null(cond)==FALSE){sel <- with(as.data.frame(sdata[[1]][[i]]), eval(parse(text=cond),envir=environment()))}else{sel<-rep(TRUE,nrow(sdata[[1]][[i]]))}
results.table[results.table$a1==X$setup$abar[i],"psi"] <- apply(as.matrix(sdata[[1]][[i]][sel,X$setup$Ynodes],ncol=length(X$setup$Ynodes)),2,fun,...)
if(X$setup$B>0){
sbdata <- lapply(sdata, '[[', i); bfs <- matrix(NA,ncol=length(X$setup$Ynodes),nrow=X$setup$B)
for(b in 2:(X$setup$B+1)){bfs[b-1,] <- apply(data.matrix(sbdata[[b]][sel,X$setup$Ynodes]),2,fun,...) } 
results.table[results.table$a1==X$setup$abar[i],c("l95","u95")] <- t(apply(bfs,2,quantile,probs=c(0.025,0.975)))
if(with.se==TRUE){
  if(i==1){results.table$se<-NA}
  results.table[results.table$a1==X$setup$abar[i],c("se")] <- apply(bfs,2,sd)} 
}
}}

if(X$setup$i.type=="custom"){
for(i in 1:length(sdata[[1]])){
if(is.null(cond)==FALSE){sel <- with(as.data.frame(sdata[[1]][[i]]), eval(parse(text=cond),envir=environment()))}else{sel<-rep(TRUE,nrow(sdata[[1]][[i]]))}
results.table[results.table$a1==X$setup$abar[i,1],"psi"] <- apply(as.matrix(sdata[[1]][[i]][sel,X$setup$Ynodes],ncol=length(X$setup$Ynodes)),2,fun,...)
if(X$setup$B>0){
sbdata <- lapply(sdata, '[[', i); bfs <- matrix(NA,ncol=length(X$setup$Ynodes),nrow=X$setup$B)
for(b in 2:(X$setup$B+1)){bfs[b-1,] <- apply(data.matrix(sbdata[[b]][sel,X$setup$Ynodes]),2,fun,...) } 
results.table[results.table$a1==X$setup$abar[i],c("l95","u95")] <- t(apply(bfs,2,quantile,probs=c(0.025,0.975)))
if(with.se==TRUE){
  if(i==1){results.table$se<-NA}
  results.table[results.table$a1==X$setup$abar[i],c("se")] <- apply(bfs,2,sd)}  
}
}}


if(X$setup$i.type=="natural"){
if(is.null(cond)==FALSE){sel_nat <- with(as.data.frame(sdata[[1]][[1]]), eval(parse(text=cond),envir=environment()))}else{sel_nat<-rep(TRUE,nrow(sdata[[1]][[1]]))}
if(is.null(cond)==FALSE){sel_obs <- with(as.data.frame(odata[[1]]), eval(parse(text=cond),envir=environment()))}else{sel_obs<-rep(TRUE,nrow(odata[[1]]))}
results.table[results.table$a1=="natural","psi"] <- apply(as.matrix(sdata[[1]][[1]][sel_nat,X$setup$Ynodes],ncol=length(X$setup$Ynodes)),2,fun,...)
results.table[results.table$a1=="observed","psi"] <- try(apply(data.frame(odata[[1]][sel_obs,X$setup$Ynodes]),2,fun,...))
  if(is.null(X$setup$Lnodes)==FALSE){
  if(verbose==TRUE){cat("Note: specified function and condition will also be applied to counterfactual Lnodes. \n\n")}
      lblocks <- rep(NA,length(X$setup$Ynodes))
    for(i in 1:(length(X$setup$Ynodes))){
    blocks <- make.interval(which(colnames(odata[[1]])%in%X$setup$Ynodes))
    lblocks[i] <- sum(which(colnames(odata[[1]])%in%X$setup$Lnodes)%in%blocks[[i]])
    }
   if(is.wholenumber(length(X$setup$Lnodes)/length(X$setup$Ynodes)) | X$setup$n.t==1){
   results.table[results.table$a1=="natural",paste0("L_",1:(length(X$setup$Lnodes)/X$setup$n.t))]  <- matrix(apply(matrix(sdata[[1]][[1]][sel_nat,X$setup$Lnodes],ncol=length(X$setup$Lnodes)),2,fun,...),ncol=(length(X$setup$Lnodes)/X$setup$n.t),byrow=T)
   results.table[results.table$a1=="observed",paste0("L_",1:(length(X$setup$Lnodes)/X$setup$n.t))] <- matrix(try(apply(data.frame(odata[[1]][sel_obs,X$setup$Lnodes]),2,fun,...)),ncol=(length(X$setup$Lnodes)/X$setup$n.t),byrow=T)
   results.table[results.table$a1=="difference",paste0("L_",1:(length(X$setup$Lnodes)/X$setup$n.t))] <- results.table[results.table$a1=="natural",paste0("L_",1:(length(X$setup$Lnodes)/X$setup$n.t))] - 
                                                                                                            results.table[results.table$a1=="observed",paste0("L_",1:(length(X$setup$Lnodes)/X$setup$n.t))]                                                                                                       
   }else{
    if(length(lblocks)>1 & length(unique(lblocks[-1]))==1){
        results.table[results.table$a1=="natural" & results.table$time!=1,paste0("L_",1:max(lblocks))] <- matrix(apply(matrix(sdata[[1]][[1]][sel_nat,X$setup$Lnodes],ncol=length(X$setup$Lnodes)),2,fun,...),ncol=max(lblocks),byrow=T)
        results.table[results.table$a1=="observed"& results.table$time!=1,paste0("L_",1:max(lblocks))] <- matrix(try(apply(data.frame(odata[[1]][sel_obs,X$setup$Lnodes]),2,fun,...)),ncol=max(lblocks),byrow=T)
        if(verbose==TRUE){cat("Note: it is not entirely clear to which time points your Lnodes belong to. I guessed what could make sense, but please check.\n\n")}                                                                                          
   }else{
      if(max(lblocks,na.rm=T)==1){
        updat.index <- unlist(lapply(apply(t(outer(which(colnames(odata[[1]])%in%X$setup$Lnodes), which(colnames(odata[[1]])%in%X$setup$Ynodes), "<")),1,which),max))
        updt.Lnodes <- apply(as.matrix(sdata[[1]][[1]][sel_nat,X$setup$Lnodes],ncol=length(X$setup$Lnodes)),2,fun,...)[updat.index]
        updt.Lnodes2 <-  try(apply(data.frame(odata[[1]][sel_obs,X$setup$Lnodes]),2,fun,...))[updat.index]
        results.table[results.table$a1=="natural",paste0("L_",1)] <- updt.Lnodes
        results.table[results.table$a1=="observed",paste0("L_",1)] <- updt.Lnodes2 
   }else{if(verbose==TRUE){cat("Note: your number of Lnodes are not a multiple of your number of Ynodes.\n This is fine, but you seem to have more than one Lnode per time point, and maybe unequally distributed. \n It is difficult for me to guess which Lnodes belong 'together', and at which time point. \n Thus, no natural course values for your Lnodes are calculated. \n Use 'ret=TRUE' and calculate manually.\n")}}
   }}
   }
results.table[results.table$a1=="difference","psi"] <- results.table[results.table$a1=="natural","psi"] - results.table[results.table$a1=="observed","psi"]
  
    
if(X$setup$B>0){
sbdata <- lapply(sdata, '[[', 1); bfs <- ofs <- matrix(NA,ncol=length(c(X$setup$Ynodes,X$setup$Lnodes)),nrow=X$setup$B,dimnames=list(NULL,c(X$setup$Ynodes,X$setup$Lnodes)))
for(b in 2:(X$setup$B+1)){
    if(is.null(cond)==FALSE){sel <- with(as.data.frame(sbdata[[b]]), eval(parse(text=cond)))}else{sel<-rep(TRUE,nrow(sbdata[[b]]))}
    bfs[b-1,] <- apply(data.matrix(sbdata[[b]][sel,c(X$setup$Ynodes,X$setup$Lnodes)]),2,fun,...)  }
for(b in 2:(X$setup$B+1)){
    if(is.null(cond)==FALSE){sel <- with(as.data.frame(odata[[b]]), eval(parse(text=cond)))}else{sel<-rep(TRUE,nrow(odata[[b]]))}
    ofs[b-1,] <- apply(data.matrix( odata[[b]][sel,c(X$setup$Ynodes,X$setup$Lnodes)]),2,fun,...) 
}
dfs <- bfs-ofs
results.table[results.table$a1=="natural",c("l95","u95")] <- t(apply(subset(bfs,select=X$setup$Ynodes),2,quantile,probs=c(0.025,0.975)))
results.table[results.table$a1=="observed",c("l95","u95")] <- t(apply(subset(ofs,select=X$setup$Ynodes),2,quantile,probs=c(0.025,0.975)))
results.table[results.table$a1=="difference",c("l95","u95")] <-    t(apply(subset(dfs,select=X$setup$Ynodes),2,quantile,probs=c(0.025,0.975)))
  if(with.se==TRUE){
  results.table$se<-NA
  results.table[results.table$a1=="natural",c("se")] <- apply(subset(bfs,select=X$setup$Ynodes),2,sd)
  results.table[results.table$a1=="observed",c("se")] <- apply(subset(ofs,select=X$setup$Ynodes),2,sd)
  results.table[results.table$a1=="difference",c("se")] <-    apply(subset(dfs,select=X$setup$Ynodes),2,sd)
  }
  if(is.null(X$setup$Lnodes)==FALSE){
    if(is.wholenumber(length(X$setup$Lnodes)/length(X$setup$Ynodes)) | X$setup$n.t==1){
    results.table[results.table$a1=="natural",grep(":l95",colnames(results.table))] <-  matrix(apply(bfs[,X$setup$Lnodes],2,quantile,probs=c(0.025)),ncol=length(X$setup$Lnodes)/X$setup$n.t,byrow=T)
    results.table[results.table$a1=="natural",grep(":u95",colnames(results.table))] <-  matrix(apply(bfs[,X$setup$Lnodes],2,quantile,probs=c(0.975)),ncol=length(X$setup$Lnodes)/X$setup$n.t,byrow=T)
    results.table[results.table$a1=="observed",grep(":l95",colnames(results.table))] <-  matrix(apply(ofs[,X$setup$Lnodes],2,quantile,probs=c(0.025)),ncol=length(X$setup$Lnodes)/X$setup$n.t,byrow=T)
    results.table[results.table$a1=="observed",grep(":u95",colnames(results.table))] <-  matrix(apply(ofs[,X$setup$Lnodes],2,quantile,probs=c(0.975)),ncol=length(X$setup$Lnodes)/X$setup$n.t,byrow=T)
    results.table[results.table$a1=="difference",grep(":l95",colnames(results.table))] <-  matrix(apply(dfs[,X$setup$Lnodes],2,quantile,probs=c(0.025)),ncol=length(X$setup$Lnodes)/X$setup$n.t,byrow=T)
    results.table[results.table$a1=="difference",grep(":u95",colnames(results.table))] <-  matrix(apply(dfs[,X$setup$Lnodes],2,quantile,probs=c(0.975)),ncol=length(X$setup$Lnodes)/X$setup$n.t,byrow=T)
        if(with.se==TRUE){
            results.table[results.table$a1=="natural",paste0("L_",1:(length(X$setup$Lnodes)/X$setup$n.t),":se")]    <-  matrix(apply(bfs[,X$setup$Lnodes],2,sd),ncol=(length(X$setup$Lnodes)/X$setup$n.t),byrow=T)
            results.table[results.table$a1=="observed",paste0("L_",1:(length(X$setup$Lnodes)/X$setup$n.t),":se")]   <-  matrix(apply(ofs[,X$setup$Lnodes],2,sd),ncol=(length(X$setup$Lnodes)/X$setup$n.t),byrow=T)
            results.table[results.table$a1=="difference",paste0("L_",1:(length(X$setup$Lnodes)/X$setup$n.t),":se")] <-  matrix(apply(dfs[,X$setup$Lnodes],2,sd),ncol=(length(X$setup$Lnodes)/X$setup$n.t),byrow=T)}
    }else{
      if(length(lblocks)>1 & length(unique(lblocks[-1]))==1){
        results.table[results.table$a1=="natural" & results.table$time!=1,grep(":l95",colnames(results.table))] <-  matrix(apply(bfs[,X$setup$Lnodes],2,quantile,probs=c(0.025)),ncol=max(lblocks),byrow=T)
        results.table[results.table$a1=="natural" & results.table$time!=1,grep(":u95",colnames(results.table))] <-  matrix(apply(bfs[,X$setup$Lnodes],2,quantile,probs=c(0.975)),ncol=max(lblocks),byrow=T)
        results.table[results.table$a1=="observed" & results.table$time!=1,grep(":l95",colnames(results.table))] <-  matrix(apply(ofs[,X$setup$Lnodes],2,quantile,probs=c(0.025)),ncol=max(lblocks),byrow=T)
        results.table[results.table$a1=="observed" & results.table$time!=1,grep(":u95",colnames(results.table))] <-  matrix(apply(ofs[,X$setup$Lnodes],2,quantile,probs=c(0.975)),ncol=max(lblocks),byrow=T)
        results.table[results.table$a1=="difference" & results.table$time!=1,grep(":l95",colnames(results.table))] <-  matrix(apply(dfs[,X$setup$Lnodes],2,quantile,probs=c(0.025)),ncol=max(lblocks),byrow=T)
        results.table[results.table$a1=="difference" & results.table$time!=1,grep(":u95",colnames(results.table))] <-  matrix(apply(dfs[,X$setup$Lnodes],2,quantile,probs=c(0.975)),ncol=max(lblocks),byrow=T)
          if(with.se==TRUE){
            results.table[results.table$a1=="natural" & results.table$time!=1,paste0("L_",1:max(lblocks),":se")]      <-  matrix(apply(bfs[,X$setup$Lnodes],2,sd),ncol=max(lblocks),byrow=T)
            results.table[results.table$a1=="observed" & results.table$time!=1,paste0("L_",1:max(lblocks),":se")]     <-  matrix(apply(ofs[,X$setup$Lnodes],2,sd),ncol=max(lblocks),byrow=T)
            results.table[results.table$a1=="difference" & results.table$time!=1,paste0("L_",1:max(lblocks),":se")]   <-  matrix(apply(dfs[,X$setup$Lnodes],2,sd),ncol=max(lblocks),byrow=T)}
      }else{
        if(max(lblocks,na.rm=T)==1){
          results.table[results.table$a1=="natural", paste(paste0("L_",1),"l95",sep=":")] <-
          matrix(apply(bfs[,X$setup$Lnodes],2,quantile,probs=c(0.025)),ncol=1,byrow=T)[updat.index]
          results.table[results.table$a1=="natural", paste(paste0("L_",1),"u95",sep=":")] <-
          matrix(apply(bfs[,X$setup$Lnodes],2,quantile,probs=c(0.975)),ncol=1,byrow=T)[updat.index]
          results.table[results.table$a1=="observed", c(paste(paste0("L_",1),"l95",sep=":"))] <-
          apply(ofs[,X$setup$Lnodes],2,quantile,probs=c(0.025))[updat.index]
          results.table[results.table$a1=="observed", c(paste(paste0("L_",1),"u95",sep=":"))] <-
          apply(ofs[,X$setup$Lnodes],2,quantile,probs=c(0.975))[updat.index]
          col.order <- c(1:grep("psi",colnames(results.table)),which(colnames(results.table)%in%c("l95","u95")))
          col.order <- c(col.order,sort(grep(paste0("L_",1),colnames(results.table))))
          if(with.se==TRUE){col.order<-c(col.order, which(colnames(results.table)%in%c("se")))}
          results.table <- results.table[,col.order]
          results.table[results.table$a1=="difference", paste(paste0("L_",1),"l95",sep=":")] <-
          apply(bfs[,X$setup$Lnodes]-ofs[,X$setup$Lnodes],2,quantile,probs=c(0.025))[updat.index]
          results.table[results.table$a1=="difference", paste(paste0("L_",1),"u95",sep=":")] <-
          apply(bfs[,X$setup$Lnodes]-ofs[,X$setup$Lnodes],2,quantile,probs=c(0.975))[updat.index]
          if(with.se==TRUE){
            results.table[results.table$a1=="natural",paste0("L_",1,":se")]      <-  matrix(apply(bfs[,X$setup$Lnodes],2,sd),ncol=1,byrow=T)[updat.index]
            results.table[results.table$a1=="observed",paste0("L_",1,":se")]     <-  matrix(apply(ofs[,X$setup$Lnodes],2,sd),ncol=1,byrow=T)[updat.index]
            results.table[results.table$a1=="difference",paste0("L_",1,":se")]   <-  matrix(apply(dfs[,X$setup$Lnodes],2,sd),ncol=1,byrow=T)[updat.index]}
         }else{
            if(verbose==TRUE){cat("Note: your number of Lnodes are not a multiple of your number of Ynodes.\n This is fine, but you seem to have more than one Lnode per time point, and maybe unequally distributed. \n It is difficult for me to guess which Lnodes belong 'together', and at which time point. \n Thus, no natural course values for your Lnodes are calculated. \n Use 'ret=TRUE' and calculate manually.\n")}
              }}
    }
    }
}                                                                                                            


}

custom.setup <- X$setup
custom.setup$measure<-"custom"
bsa <- ofa <- dfa <- NULL; if(with.se==TRUE){bsa<-bfs
                if(X$setup$i.type=="natural"){ofa<-ofs; dfa<-dfs}}

results <- list(results=results.table,diagnostics=X$diagnostics,setup=custom.setup,b.results = bsa, b.results2 = list(ofa,dfa))
class(results) <- "gformula"
results
}