make.model.formulas <- function(X, Ynodes=NULL, Lnodes=NULL, Cnodes=NULL, Anodes=NULL, survival=FALSE, evaluate=FALSE){

if(survival==T & is.null(Cnodes)){stop("If you specify survival=T, then Cnodes can not be NULL")}
m.model.families <- assign.family(X)
A.model.formulas.m <- C.model.formulas.m <- L.model.formulas.m <- Y.model.formulas.m <-NULL
fitted.models <- fitted.model.summary <- NULL
fitted.model.summary.L <- fitted.model.summary.Y <- fitted.model.summary.A <- fitted.model.summary.C <- NULL
m.model.A.names <- m.model.C.names <-  m.model.Y.names <- m.model.L.names <- "ph"; ph <- NULL


if(is.null(Lnodes)==FALSE){
m.model.L.names  <- paste("m","_",colnames(X)[colnames(X)%in%c(Lnodes)],sep="")
m.model.L.Ynodes <- colnames(X)[colnames(X)%in%c(Lnodes)]
m.model.L.families <- m.model.families[colnames(X)%in%c(Lnodes)]
n.L.models.m <-  sum(colnames(X)%in%c(Lnodes))
L.model.formulas.m <- c(rep(NA,n.L.models.m))
for(j in 1:n.L.models.m){
L.data <- subset(X,select=c(1:which(colnames(X)%in%m.model.L.Ynodes[j])))
if(survival==T |(survival==F & is.null(Cnodes)==F)){
  siC <- which(colnames(L.data)%in%Cnodes[1:j]); if(length(siC)==0){siC<-NULL}
  if(survival==T){siY <- which(colnames(L.data)%in%Ynodes[1:j]); if(length(siY)==0){siY<-NULL}}else{siY<-NULL}
  si <- c(siY,siC)
  if(is.null(si)==FALSE){L.data <- subset(L.data, select=-si)}
}
L.model.formulas.m[j] <-   Reduce(paste, deparse(make.formula(L.data,approach="GLM",fam=m.model.L.families[j])))
if(evaluate==TRUE){assign(m.model.L.names[j],mgcv::gam(make.formula(L.data,approach="GLM",fam=m.model.L.families[j]),data=L.data,control=list(keepData=T),family=m.model.L.families[j]))}
}
}

if(is.null(Ynodes)==FALSE){
m.model.Y.names  <- paste("m","_",colnames(X)[colnames(X)%in%c(Ynodes)],sep="")
m.model.Y.Ynodes <- colnames(X)[colnames(X)%in%c(Ynodes)]
m.model.Y.families <- m.model.families[colnames(X)%in%c(Ynodes)]
n.Y.models.m <-  sum(colnames(X)%in%c(Ynodes))
Y.model.formulas.m <-  c(rep(NA,n.Y.models.m))
for(j in 1:n.Y.models.m){
Y.data <- subset(X,select=c(1:which(colnames(X)%in%m.model.Y.Ynodes[j])))
if(survival==T |(survival==F & is.null(Cnodes)==F)){
  if(survival==T){if(j>1){siY <- which(colnames(Y.data)%in%Ynodes[1:(j-1)])}else{siY<-NULL}}else{siY<-NULL} 
  siC <- which(colnames(Y.data)%in%Cnodes[1:j]); if(length(siC)==0){siC<-NULL}
  si <- c(siY,siC)
  if(is.null(si)==FALSE){Y.data <- subset(Y.data, select=-si)}
  }
Y.model.formulas.m[j] <-   Reduce(paste, deparse(make.formula(Y.data,approach="GLM",fam=m.model.Y.families[j])))
if(evaluate==TRUE){assign(m.model.Y.names[j],mgcv::gam(make.formula(Y.data,approach="GLM",fam=m.model.Y.families[j]),data=Y.data, control=list(keepData=T),family=m.model.Y.families[j]))}
}
}


if(is.null(Anodes)==FALSE){
m.model.A.names  <- paste("m","_",colnames(X)[colnames(X)%in%c(Anodes)],sep="")
m.model.A.Ynodes <- colnames(X)[colnames(X)%in%c(Anodes)]
m.model.A.families <- m.model.families[colnames(X)%in%c(Anodes)]
n.A.models.m <-  sum(colnames(X)%in%c(Anodes))
A.model.formulas.m <- c(rep(NA,n.A.models.m))
for(j in 1:n.A.models.m){
A.data <- subset(X,select=c(1:which(colnames(X)%in%m.model.A.Ynodes[j])))
if(survival==T |(survival==F & is.null(Cnodes)==F)){
  siC <- which(colnames(A.data)%in%Cnodes[1:j]); if(length(siC)==0){siC<-NULL}
  if(survival==T){siY <- which(colnames(A.data)%in%Ynodes[1:j]); if(length(siY)==0){siY<-NULL}}else{siY<-NULL}
  si <- c(siY,siC)
  if(is.null(si)==FALSE){A.data <- subset(A.data, select=-si)}
}
A.model.formulas.m[j] <-   Reduce(paste, deparse(make.formula(A.data,approach="GLM",fam=m.model.A.families[j])))
if(evaluate==TRUE){assign(m.model.A.names[j],mgcv::gam(make.formula(A.data,approach="GLM",fam=m.model.A.families[j]),data=A.data,family=m.model.A.families[j]))}
}
}

if(is.null(Cnodes)==FALSE){
m.model.C.names  <- paste("m","_",colnames(X)[colnames(X)%in%c(Cnodes)],sep="")
m.model.C.Ynodes <- colnames(X)[colnames(X)%in%c(Cnodes)]
m.model.C.families <- m.model.families[colnames(X)%in%c(Cnodes)]
n.C.models.m <-  sum(colnames(X)%in%c(Cnodes))
C.model.formulas.m <- c(rep(NA,n.C.models.m))
for(j in 1:n.C.models.m){
C.data <- subset(X,select=c(1:which(colnames(X)%in%Cnodes[j])))
if(survival==T){
  if(j>1){siC <- which(colnames(C.data)%in%Cnodes[1:(j-1)])}else{siC<-NULL} 
  siY <- which(colnames(C.data)%in%Ynodes[1:j]); if(length(siY)==0){siY<-NULL}
  si <- c(siY,siC)
  if(is.null(si)==FALSE){C.data <- subset(C.data, select=-si)}
}else{
  if(j>1){si <- which(colnames(C.data)%in%Cnodes[1:(j-1)])}else{si<-NULL}  
  if(is.null(si)==FALSE){C.data <- subset(C.data, select=-si)}
}
C.model.formulas.m[j] <-   Reduce(paste, deparse(make.formula(C.data,approach="GLM",fam=m.model.C.families[j])))
if(evaluate==TRUE){assign(m.model.C.names[j],mgcv::gam(make.formula(C.data,approach="GLM",fam=m.model.C.families[j]),data=C.data,family=m.model.C.families[j]))}
}}



model.names <- list(Lnames=L.model.formulas.m,Ynames=Y.model.formulas.m,Anames=A.model.formulas.m,Cnames=C.model.formulas.m)
for(i in 1:4){if(is.null(model.names[[i]])==FALSE){for(j in 1:length(model.names[[i]])){model.names[[i]][j]<-gsub(" ", "",model.names[[i]][j])}}}
if(evaluate==TRUE){fitted.models <- list(fitted.L=mget(m.model.L.names),fitted.Y=mget(m.model.Y.names),
fitted.A=mget(m.model.A.names),fitted.C=mget(m.model.C.names))}
if(is.null(Lnodes)==FALSE){if(evaluate==TRUE){fitted.model.summary.L <- lapply(fitted.models$fitted.L,summary)}}
if(is.null(Ynodes)==FALSE){if(evaluate==TRUE){fitted.model.summary.Y <- lapply(fitted.models$fitted.Y,summary)}}
if(is.null(Anodes)==FALSE){if(evaluate==TRUE){fitted.model.summary.A <- lapply(fitted.models$fitted.A,summary)}}
if(is.null(Cnodes)==FALSE){if(evaluate==TRUE){fitted.model.summary.C <- lapply(fitted.models$fitted.C,summary)}}
fitted.model.summary <- list(L=fitted.model.summary.L,Y=fitted.model.summary.Y,A=fitted.model.summary.A,C=fitted.model.summary.C)

return(list(model.names=model.names,fitted.models=fitted.models,fitted.model.summary=fitted.model.summary))
}


