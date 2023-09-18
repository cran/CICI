fit.updated.formulas <- function(formulas,X){
all.models <- rep(list(NULL),length(formulas))
for(i in 1:length(all.models)){
  if(is.null(formulas[[i]])==FALSE){
    all.models[[i]] <- rep(list(NULL),length(formulas[[i]]))
    families <- extract.families(formulas[[i]],fdata=X)
    for(j in 1:length(all.models[[i]])){
      if(substr(families[j],1,4)!="mult"){gf<-formula(formulas[[i]][j])}else{gf<-eval(parse(text=formulas[[i]][j]))}
      all.models[[i]][[j]] <- try(gam(gf,data=X,family=families[j], control=list(keepData=T)),silent=T) # note: keepData argument sometimes produced hidden(?) error message related to cache
      if(class(all.models[[i]][[j]])[1]=="try-error"){all.models[[i]][[j]] <- "updated model could not be fitted"}
}}}
all.summaries <- all.models
robust.summary <- function(go){return(tryCatch(mgcv::summary.gam(go), error=function(e) "model could not be fitted"))}
for(i in 1:length(all.summaries)){if(is.null(all.summaries[[i]])==FALSE){all.summaries[[i]]<-try(lapply(all.summaries[[i]],robust.summary),silent=TRUE)}}
names(all.models) <- names(all.summaries) <- names(formulas)
return(list(fitted.models=all.models,all.summaries=all.summaries))
}