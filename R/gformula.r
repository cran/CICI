# X <- Odat[,-1]
# Lnodes  <- colnames(X)[sort(c(grep("L1_",colnames(X)),grep("L2_",colnames(X))))]
# Ynodes  <- colnames(X)[sort(c(grep("Y_",colnames(X))))]
# Anodes  <- colnames(X)[sort(c(grep("A_",colnames(X))))]
# Cnodes  <- colnames(X)[sort(c(grep("C_",colnames(X))))]
# survivalY <- FALSE
# a.min <- round(min(X[,grep("A_",colnames(X))],na.rm=T)-1,digits=0);a.max <- round(max(X[,grep("A_",colnames(X))],na.rm=T)+1,digits=0);n.int<-10
# abar <- seq(a.min,a.max,length.out=n.int) # could also be a matrix of size n.int times t
# cbar <- "uncensored"
# Yform <- Aform <- Cform <- Lform <-  "GLM"
# calc.support = FALSE
# B = 0
# verbose=TRUE
# ret=F
# ncores=1
# prog=NULL
# seed = NULL
###
# optional: plot Strategy names for custom interventions at bottom?
# 6) calc.support for custom interventions: does not make sense. Not allow?
# 7) natural bounds 
# 8) shift 
# 9) SACE ; and competing() and revisit definitions of cbar
# 10) visit variables 

gformula <- function(X, Anodes, Ynodes, Lnodes = NULL, Cnodes = NULL,
                        abar =  NULL, cbar = "uncensored", survivalY = FALSE,
                        Yform = "GLM", Lform = "GLM", Aform = "GLM", Cform="GLM",
                        calc.support = FALSE, B = 0, ret=FALSE,  
                        ncores=1, verbose=TRUE, seed=NULL, prog=NULL,...){

### checks and setup ###
if(is.null(prog)==FALSE){write(matrix("started with setup..."),file=paste(prog,"/progress.txt",sep=""))}
model.families <- assign.family(X)

#
if(is.null(seed)==FALSE){set.seed(seed)}
if(cbar[1]!="uncensored" & cbar[1]!="natural" & calc.support==TRUE){calc.support<-FALSE; if(verbose==TRUE){cat("Note: For custom cbar interventions no support is calculated so far \n")}}
if(is.data.frame(X)==FALSE){stop("'X' (i.e. your data) needs to be a data.frame")}
if(missing(Anodes)){stop("'Anodes' is missing. Please specify your intervention node(s).")}  
if(missing(Ynodes)){stop("'Ynodes' is missing. Please specify your outcome node(s).")}
if(min(which(colnames(X)%in%Ynodes))<min(which(colnames(X)%in%Anodes))){stop("Ynodes can not occur before Anodes.\n Likely you specified pre-intervention variables as Ynode?")}
if(is.character(Anodes)==FALSE & is.null(Anodes)==FALSE){stop("'Anodes' needs to be a character vector containing the respective variable names.")}
if(is.character(Cnodes)==FALSE & is.null(Cnodes)==FALSE){stop("'Cnodes' needs to be a character vector containing the respective variable names.")}
if(is.character(Lnodes)==FALSE & is.null(Lnodes)==FALSE){stop("'Lnodes' needs to be a character vector containing the respective variable names.")}
if(is.character(Ynodes)==FALSE & is.null(Ynodes)==FALSE){stop("'Ynodes' needs to be a character vector containing the respective variable names.")}
if(is.null(Lnodes)==FALSE){if(min(which(colnames(X)%in%Lnodes))<min(which(colnames(X)%in%Anodes))){
    base.L <- colnames(X)[which(colnames(X)%in%Lnodes)[which(colnames(X)%in%Lnodes)<min(which(colnames(X)%in%Anodes))]]
    Lnodes <- Lnodes[which(!Lnodes%in%base.L)] ; if(length(Lnodes)==0){Lnodes<-NULL}
if(verbose==TRUE){cat(paste("Note: Lnodes are not supposed to appear before Anodes.\n The following Lnodes are treated as pre-intervention variables:",paste(base.L,collapse=" "),"\n"))}
}}
if(any(colnames(X)[min(which(colnames(X)%in%Anodes)):ncol(X)]%in%c(Anodes,Lnodes,Cnodes, Ynodes)==FALSE)){
  stop(paste("The following variables are part of your post-intervention data, but not specififed as L-, A-, C- or Y-node:", 
        paste(colnames(X)[min(which(colnames(X)%in%Anodes)):ncol(X)][colnames(X)[min(which(colnames(X)%in%Anodes)):ncol(X)]%in%c(Anodes,Lnodes,Cnodes, Ynodes)==FALSE], collapse=" ")))
}
if(any(substr(colnames(X),1,4)=="list")){stop("Variable names that start with 'list' are not allowed. Please rename.")}
if(any(c(Ynodes,Lnodes,Anodes,Cnodes)%in%colnames(X)==FALSE)){stop(paste("You have specified the following variable name(s) which are not part of the data:",
             paste(c(Ynodes,Lnodes,Anodes,Cnodes)[c(Ynodes,Lnodes,Anodes,Cnodes)%in%colnames(X)==FALSE],collapse=" ") ,"\n"))}
if(any(sapply(X,is.ordered))){stop("Ordered variables are not allowed")}
if(any(abar=="natural") & calc.support==TRUE){calc.support<-FALSE; if(verbose==TRUE){cat("Note: no support can be calculated under natural interventions.\n")}}
if(is.matrix(abar)==FALSE & (is.vector(abar)&is.numeric(abar))==FALSE & any(abar=="natural")==FALSE){stop("'abar' not correctly specified. Type ?gformula for help.")}
if(is.logical(survivalY)==FALSE){stop("survivalY needs to be TRUE or FALSE")}
if(survivalY==TRUE){if(verbose==TRUE){if(is.null(Cnodes)){cat("Warning: you indicate that you have survival data (survivalY=T), but you have no Cnodes specified. \n")}}}
if(is.character(Yform)==FALSE){stop("Yform needs to be a character vector")}
if(is.character(Cform)==FALSE){stop("Cform needs to be a character vector")}
if(is.character(Lform)==FALSE){stop("Lform needs to be a character vector")}
if(is.character(Aform)==FALSE){stop("Aform needs to be a character vector")}
if(!Yform[1]%in%c("GLM","GAM") & length(Yform)!=length(Ynodes)){stop(paste("You have provided",length(Yform),"model formula(s); but there are",length(Ynodes),"Ynodes"))}
if(!Lform[1]%in%c("GLM","GAM") & length(Lform)!=length(Lnodes)){stop(paste("You have provided",length(Lform),"model formula(s); but there are",length(Lnodes),"Lnodes"))}
if(!Aform[1]%in%c("GLM","GAM") & length(Aform)!=length(Anodes)){stop(paste("You have provided",length(Aform),"model formula(s); but there are",length(Anodes),"Anodes"))}
if(!Cform[1]%in%c("GLM","GAM") & length(Cform)!=length(Cnodes)){stop(paste("You have provided",length(Cform),"model formula(s); but there are",length(Cnodes),"Cnodes"))}
if(any(cbar!="uncensored") & any(cbar!="natural") & survivalY==TRUE){stop("custom cbar interventions for survival data currently not supported")} #check again 
if(is.logical(ret)==FALSE){stop("'ret' needs to either TRUE or FALSE")}
if(is.logical(verbose)==FALSE){stop("'verbose' needs to either TRUE or FALSE")}
if(is.numeric(ncores)==FALSE){stop("'ncores' needs to be numeric")}
if(is.character(prog)==FALSE & is.null(prog)==FALSE){stop("'prog' needs to be a character vector")}

if(any(model.families=="binomial")){bin.problem <- !sapply(subset(X,select=(model.families=="binomial")), is.binary)
    if(any(bin.problem)){if(verbose==TRUE){cat(paste("Binary variables have been recoded:",paste(names(bin.problem)[bin.problem],collapse=","),"\n"))}
    X[,names(bin.problem)[bin.problem]] <- data.frame(sapply(subset(X,select=names(bin.problem)[bin.problem]),binary.to.zeroone,verb=verbose)) }
}
if(any(substr(model.families,1,4)=="mult")){fac.problem <- !sapply(subset(X,select=(substr(model.families,1,4)=="mult")), right.coding)
    if(any(fac.problem)){if(verbose==TRUE){cat(paste("Categorical variables have been recoded:",paste(names(fac.problem)[fac.problem],collapse=","),"\n"))}
    X[,names(fac.problem)[fac.problem]] <- data.frame(sapply(subset(X,select=names(fac.problem)[fac.problem]),factor.to.numeric,verb=verbose)) }
}

if(verbose==TRUE){
if(any(abar=="natural")==FALSE){rel.nodes<-c(Ynodes,Lnodes)}else{rel.nodes<-c(Ynodes,Lnodes,Cnodes,Anodes)}
 gvar <- colnames(X)[model.families=="gaussian"]; bvar <- colnames(X)[model.families=="binomial"]; pvar <- colnames(X)[model.families=="poisson"]; mvar <- colnames(X)[substr(model.families,1,4)=="mult"] 
 cat(paste( "Note:\n linear regression used for:", paste(gvar[gvar%in%rel.nodes],"",collapse=""),"\n", "logistic regression used for:", paste(bvar[bvar%in%rel.nodes],"",collapse=""),"\n", 
 "Poisson regression used for:", paste(pvar[pvar%in%rel.nodes],"",collapse=""),"\n",  "Multinomial regression used for:", paste(mvar[mvar%in%rel.nodes],"",collapse=""),"\n"))
 baseline.var <-  colnames(X)[1:(max(1,which(colnames(X)%in%c((colnames(X)[colnames(X)%in%c(Anodes)])[1]))-1))]
  cat(paste(" Pre-Intervention variables are:",paste(baseline.var,collapse=" "),"\n\n"))
}
#



### matrices to store results ###
n.a <- length(Anodes); n.t <- length(Ynodes); n.l <-length(Lnodes); time.points <- 1:n.t
if(is.matrix(abar)==TRUE){interventions <- do.call(rbind, replicate(n.t, abar, simplify=FALSE)); i.type <- "custom"}else{
                            if(is.vector(abar)){interventions <- do.call(rbind, replicate(n.t, matrix(rep(abar,n.a),ncol=n.a), simplify=FALSE)); i.type<-"standard"}
                            if(any(abar=="natural")){interventions <- matrix(c(rep("natural",n.a),rep("observed",n.a),rep("difference",n.a)),ncol=n.a,byrow=T); i.type="natural"}
                            }
if(i.type!="natural"){n.int <- dim(interventions)[1]/n.t}else{n.int <- 2}
if(i.type=="natural" & calc.support==TRUE){calc.support<-FALSE; if(verbose==TRUE){cat("Note: For natural interventions no support is calculated \n")}}
if(is.matrix(cbar)==FALSE){if(cbar=="uncensored"){cbar <- matrix(0,nrow=nrow(X),ncol=length(Cnodes))}}


store.results           <- as.data.frame(matrix(NA,nrow=(n.t)*n.int+as.numeric(i.type=="natural")*n.t,ncol=1+n.a+1))
colnames(store.results) <- c("time",paste("a",1:n.a ,sep=""),"psi")
store.results$time      <- rep(time.points ,each=n.int+as.numeric(i.type=="natural"))
if(length(2:(n.a+1))>1){store.results[,2:(n.a+1)] <-  interventions}else{store.results[,2:(n.a+1)] <-  c(interventions)}

if(length(Ynodes)==length(Anodes)){needed <- !(store.results$time<matrix(rep(time.points,nrow(store.results)),nrow=nrow(store.results),byrow=T))}else{
    if(i.type!="natural"){needed <-  matrix(which(colnames(X)%in%Anodes),nrow=n.int*length(time.points),ncol=length(which(colnames(X)%in%Anodes)),byrow=T) <
               matrix(rep(which(colnames(X)%in%Ynodes),each=n.int),nrow=n.int*length(time.points),ncol=n.a)  }else{
               needed <-  matrix(which(colnames(X)%in%Anodes),nrow=(n.int+1)*length(time.points),ncol=length(which(colnames(X)%in%Anodes)),byrow=T) <
                          matrix(rep(which(colnames(X)%in%Ynodes),each=(n.int+1)),nrow=(n.int+1)*length(time.points),ncol=n.a) 
               }
}
store.results[,2:(n.a+1)][needed==FALSE] <- NA


### Parallelization & Setup
if(ncores>1){
  if(ncores > parallel::detectCores()){
     ncores <- parallel::detectCores();if(verbose==TRUE){cat(paste("Note: You only have",ncores,"threads which can be utilized. \n"))}
     }
     if(verbose==TRUE){cat(paste("Note: You initialized parallel computation using",ncores,"threads...initializing cluster now... \n"))}
     cl <- parallel::makeCluster(ncores); doParallel::registerDoParallel(cl)
     exp.var <- c("make.formula","calculate.support", "extract.families", "projection_linear", 
     "screen.cramersv","screen.glmnet.cramer", "censor", "adjust.sim.surv", "multiResultClass","multi.help","rmulti","rmean")
}else{exp.var=NULL; foreach::registerDoSEQ()}

if(i.type!="natural"){int.index <- (as.numeric(rownames(store.results[store.results$time==n.t,][1,]))):nrow(store.results)}else{int.index<- 1}#nrow(store.results)-1}
sim.data <- rep(list(rep(list(NA),length(int.index))),B+1); obs.data<-rep(list(NA),B+1);obs.data[[1]]<-X
if(B>0){analysis.b<- rep(list(NA),B)}else{analysis.b<-NULL}
if(is.null(prog)==FALSE){write(matrix("started with g-formula calculations in original data...\n"),file=paste(prog,"/progress.txt",sep=""),append=TRUE)}

# Prepare support measures, if relevant
if(calc.support==TRUE){
  support <- calculate.support(dat=X,A=Anodes,intervention=interventions[store.results$time==1,],...)
  updat.index2 <- unlist(lapply(apply(t(outer(which(colnames(X)%in%Anodes), which(colnames(X)%in%Ynodes), "<")),1,which),max))   
}

### ANALYSIS ###

analysis <- foreach(i = int.index, .export=exp.var) %dorng% try({     

mydat <- X

model.L.names  <- paste("m",i,"_",colnames(mydat)[colnames(mydat)%in%c(Lnodes)],sep="")
model.Y.names  <- paste("m",i,"_",colnames(mydat)[colnames(mydat)%in%c(Ynodes)],sep="")
model.L.Ynodes <- colnames(mydat)[colnames(mydat)%in%c(Lnodes)]
model.Y.Ynodes <- colnames(mydat)[colnames(mydat)%in%c(Ynodes)]
model.L.families <- model.families[colnames(mydat)%in%c(Lnodes)]
model.Y.families <- model.families[colnames(mydat)%in%c(Ynodes)]
n.L.models <-  sum(colnames(mydat)%in%c(Lnodes))
n.Y.models <-  sum(colnames(mydat)%in%c(Ynodes))

# Step 1: fit models
for(j in 1:n.L.models)try({
L.data <- data.frame(subset(mydat,select=c(1:which(colnames(mydat)%in%model.L.Ynodes[j]))))
assign(model.L.names[j],mgcv::gam(make.formula(L.data,approach=Lform,index=j,fam=model.L.families[j]),data=L.data,family=model.L.families[j]))   
},silent=TRUE)   
for(j in 1:n.Y.models)try({
Y.data <- data.frame(subset(mydat,select=c(1:which(colnames(mydat)%in%model.Y.Ynodes[j]))))
assign(model.Y.names[j],mgcv::gam(make.formula(Y.data,approach=Yform,index=j,fam=model.Y.families[j]),data=Y.data,family=model.Y.families[j]))    
},silent=TRUE)

if(i.type=="natural"){
  model.A.names  <- paste("m",i,"_",colnames(mydat)[colnames(mydat)%in%c(Anodes)],sep="")
  model.A.Ynodes <- colnames(mydat)[colnames(mydat)%in%c(Anodes)]
  model.A.families <- model.families[colnames(mydat)%in%c(Anodes)]
  n.A.models <-  sum(colnames(mydat)%in%c(Anodes))
    for(j in 1:n.A.models)try({
    A.data <- data.frame(subset(mydat,select=c(1:which(colnames(mydat)%in%model.A.Ynodes[j]))))
    assign(model.A.names[j],mgcv::gam(make.formula(A.data,approach=Aform,index=j,fam=model.A.families[j]),data=A.data,family=model.A.families[j]))   
    },silent=TRUE)
 if(is.null(Cnodes)==FALSE){
  model.C.names  <- paste("m",i,"_",colnames(mydat)[colnames(mydat)%in%c(Cnodes)],sep="")
  model.C.Ynodes <- colnames(mydat)[colnames(mydat)%in%c(Cnodes)]
  model.C.families <- model.families[colnames(mydat)%in%c(Cnodes)]
  n.C.models <-  sum(colnames(mydat)%in%c(Cnodes))
    for(j in 1:n.C.models)try({
    C.data <- data.frame(subset(mydat,select=c(1:which(colnames(mydat)%in%model.C.Ynodes[j]))))
    assign(model.C.names[j],mgcv::gam(make.formula(C.data,approach=Cform,index=j,fam=model.C.families[j]),data=C.data,family=model.C.families[j]))   
    },silent=TRUE)
 }}

# Step 2: intervene
gdata <- mydat
first.treatment <-  (colnames(mydat)[colnames(mydat)%in%c(Anodes)])[1]
all.Anodes  <-   which(colnames(mydat)%in%c(Anodes))
all.Cnodes  <-   which(colnames(mydat)%in%c(Cnodes))
gdata[,which(colnames(mydat)%in%first.treatment):ncol(gdata)] <- NA
if(i.type!="natural"){
gdata[,all.Anodes] <-   (store.results[i,2:(n.a+1)])[is.na(store.results[i,2:(n.a+1)])==F]
gdata[,all.Cnodes] <-    cbar}

# Step 3: simulate
sim.nodes <- c(Lnodes,Ynodes); if(i.type=="natural"){sim.nodes <- c(sim.nodes,Anodes,Cnodes)}
model.order <- paste("m",i,"_",colnames(mydat)[colnames(mydat)%in%sim.nodes],sep="")
var.order   <- colnames(mydat)[colnames(mydat)%in%sim.nodes]
for(j in 1:length(var.order)){
sim.family <- model.families[which(colnames(mydat)%in%sim.nodes)][j]
 if(sim.family%in%c("gaussian","binomial")){pm   <- predict(get(model.order[j]),newdata=gdata,type="response")}else{
                                            newd <- gdata[,!apply((apply(gdata,2,is.na)),2,all)]
                                            pm   <- predict(get(model.order[j]),newdata=newd,type="response")}
if(sim.family=="gaussian"){gdata[,var.order[j]] <- rnorm(length(pm),mean=pm,sd=sigma(get(model.order[j]))^2)}else{
  if(sim.family=="binomial"){gdata[,var.order[j]] <- rbinom(n=length(pm),1,prob=pm)}else{
    if(sim.family=="poisson"){gdata[,var.order[j]]  <- rpois(n=length(pm),lambda=pm)}else{
      if(substr(sim.family,1,4)=="mult"){if(any(is.na(pm))){pm[is.na(pm)]<-1e08;if(verbose==TRUE){cat("\n Note: multinomial prediction led to missing values which were replaced by '0' \n")}}
                                         gdata[,var.order[j]] <- rmulti(pm)}
    }}}
}
if(is.null(Cnodes)==FALSE){gdata <- t(apply(gdata,1,censor,C.index=which(colnames(X)%in%Cnodes)))}
if(survivalY==TRUE){gdata <- t(apply(gdata,1,censor,C.index=which(colnames(X)%in%Ynodes)))
                   gdata<-adjust.sim.surv(gdata,Yn=Ynodes)}

# Step 4: estimate psi
if(i.type!="natural"){res.nodes<-c(Ynodes)}else{res.nodes<-c(Ynodes,Lnodes)}
if(verbose==TRUE & i==min(int.index)){
if(any(apply(t(matrix(model.families,nrow=1,dimnames=list(NULL,colnames(X)))[,res.nodes]),2,substr,start=1,stop=4 )=="mult")){
  cat("Note: for reporting the results, your categorical variables have been made numerical and the mean is reported.\n Use 'ret=T' and 'custom.measures()' for meaningful output.\n\n")}
}
results1  <- apply(subset(gdata,select=colnames(mydat)%in%res.nodes),2,mean,na.rm=TRUE)

# return results
results2 <- gdata
results <- multiResultClass(); results$result1 <- results1; results$result2 <- results2
return(results)
},silent=TRUE)

##########
# Step 5: Bootstrapping
if(B>0){
if(verbose==TRUE){cat("starting with bootstrapping \n")};if(is.null(prog)==FALSE){write(matrix("started with bootstrapping...\n"),file=paste(prog,"/progress.txt",sep=""),append=TRUE)}
rng <- rngtools::RNGseq(B*length(int.index), seed); r <- NULL
b.index <- apply(matrix(rep(1:nrow(X),B),ncol=B), 2, sample, replace=TRUE)
boots <- foreach(b = 1:B) %:%
            foreach(i = int.index, r=rng[(b-1)*length(int.index) + 1:length(int.index)],
                    .export=exp.var, .errorhandling="pass") %dopar% {
                if(is.null(seed)==FALSE){rngtools::setRNG(r)}
                mydat <- X[b.index[,b],]
                if(is.null(prog)==FALSE & i==min(int.index)){try(write(matrix(paste("performing calculations on bootstrap sample",b)),file=paste(prog,"/progress.txt",sep=""),append=T))}
                if(verbose==TRUE){if(i==min(int.index)){cat(paste0("...",b));if(b%%10==0){cat("\n")} }}
                ##############################################
                model.L.names  <- paste("m",i,"_",colnames(mydat)[colnames(mydat)%in%c(Lnodes)],sep="")
                model.Y.names  <- paste("m",i,"_",colnames(mydat)[colnames(mydat)%in%c(Ynodes)],sep="")
                model.L.Ynodes <- colnames(mydat)[colnames(mydat)%in%c(Lnodes)]
                model.Y.Ynodes <- colnames(mydat)[colnames(mydat)%in%c(Ynodes)]
                model.L.families <- model.families[colnames(mydat)%in%c(Lnodes)]
                model.Y.families <- model.families[colnames(mydat)%in%c(Ynodes)]
                n.L.models <-  sum(colnames(mydat)%in%c(Lnodes))
                n.Y.models <-  sum(colnames(mydat)%in%c(Ynodes))
                # 1
                for(j in 1:n.L.models)try({
                L.data <- data.frame(subset(mydat,select=c(1:which(colnames(mydat)%in%model.L.Ynodes[j]))))
                assign(model.L.names[j],mgcv::gam(make.formula(L.data,approach=Lform,index=j,fam=model.L.families[j]),data=L.data,family=model.L.families[j],drop.unused.levels=FALSE))   
                },silent=TRUE)   # what if Lnodes=NULL?
                for(j in 1:n.Y.models)try({
                Y.data <- data.frame(subset(mydat,select=c(1:which(colnames(mydat)%in%model.Y.Ynodes[j]))))
                assign(model.Y.names[j],mgcv::gam(make.formula(Y.data,approach=Yform,index=j,fam=model.Y.families[j]),data=Y.data,family=model.Y.families[j],drop.unused.levels=FALSE))    
                },silent=TRUE)
                if(i.type=="natural"){
                model.A.names  <- paste("m",i,"_",colnames(mydat)[colnames(mydat)%in%c(Anodes)],sep="")
                model.A.Ynodes <- colnames(mydat)[colnames(mydat)%in%c(Anodes)]
                model.A.families <- model.families[colnames(mydat)%in%c(Anodes)]
                n.A.models <-  sum(colnames(mydat)%in%c(Anodes))
                for(j in 1:n.A.models)try({
                A.data <- data.frame(subset(mydat,select=c(1:which(colnames(mydat)%in%model.A.Ynodes[j]))))
                assign(model.A.names[j],mgcv::gam(make.formula(A.data,approach=Aform,index=j,fam=model.A.families[j]),data=A.data,family=model.A.families[j],drop.unused.levels=FALSE))   
                },silent=TRUE)
                if(is.null(Cnodes)==FALSE){
                model.C.names  <- paste("m",i,"_",colnames(mydat)[colnames(mydat)%in%c(Cnodes)],sep="")
                model.C.Ynodes <- colnames(mydat)[colnames(mydat)%in%c(Cnodes)]
                model.C.families <- model.families[colnames(mydat)%in%c(Cnodes)]
                n.C.models <-  sum(colnames(mydat)%in%c(Cnodes))
                for(j in 1:n.C.models)try({
                C.data <- data.frame(subset(mydat,select=c(1:which(colnames(mydat)%in%model.C.Ynodes[j]))))
                assign(model.C.names[j],mgcv::gam(make.formula(C.data,approach=Cform,index=j,fam=model.C.families[j]),data=C.data,family=model.C.families[j],drop.unused.levels=FALSE))   
                },silent=TRUE)
                }}
                # 2
                gdata <- mydat
                first.treatment <-  (colnames(mydat)[colnames(mydat)%in%c(Anodes)])[1]
                all.Anodes  <-   which(colnames(mydat)%in%c(Anodes))
                all.Cnodes  <-   which(colnames(mydat)%in%c(Cnodes))
                gdata[,which(colnames(mydat)%in%first.treatment):ncol(gdata)] <- NA
                if(i.type!="natural"){
                gdata[,all.Anodes] <-   (store.results[i,2:(n.t+1)])[is.na(store.results[i,2:(n.t+1)])==F]
                gdata[,all.Cnodes] <-    cbar}
                # 3
                sim.nodes <- c(Lnodes,Ynodes); if(i.type=="natural"){sim.nodes <- c(sim.nodes,Anodes,Cnodes)}
                model.order <- paste("m",i,"_",colnames(mydat)[colnames(mydat)%in%sim.nodes],sep="")
                var.order   <- colnames(mydat)[colnames(mydat)%in%sim.nodes]
                for(j in 1:length(var.order)){
                sim.family <- model.families[which(colnames(mydat)%in%sim.nodes)][j]
                if(sim.family%in%c("gaussian","binomial")){pm   <- predict(get(model.order[j]),newdata=gdata,type="response")}else{
                                            newd <- gdata[,!apply((apply(gdata,2,is.na)),2,all)]
                                            pm   <- predict(get(model.order[j]),newdata=newd,type="response")}
                if(sim.family=="gaussian"){gdata[,var.order[j]] <- rnorm(length(pm),mean=pm,sd=sigma(get(model.order[j]))^2)}else{
                if(sim.family=="binomial"){gdata[,var.order[j]] <- rbinom(n=length(pm),1,prob=pm)}else{
                if(sim.family=="poisson"){gdata[,var.order[j]]  <- rpois(n=length(pm),lambda=pm)}else{
                if(substr(sim.family,1,4)=="mult"){if(any(is.na(pm))){pm[is.na(pm)]<-1e08;if(verbose==TRUE){cat("\n Note: multinomial prediction led to missing values which were replaced by '0'")}}
                                                   gdata[,var.order[j]] <- rmulti(pm)}
                }}}
                }
                if(is.null(Cnodes)==FALSE){gdata <- t(apply(gdata,1,censor,C.index=which(colnames(X)%in%Cnodes)))}
                if(survivalY==TRUE){gdata <- t(apply(gdata,1,censor,C.index=which(colnames(X)%in%Ynodes)))
                                    gdata<-adjust.sim.surv(gdata,Yn=Ynodes)}
                # 4
                if(i.type!="natural"){res.nodes<-c(Ynodes)}else{res.nodes<-c(Ynodes,Lnodes)}
                results1  <- apply(subset(gdata,select=colnames(mydat)%in%res.nodes),2,mean,na.rm=TRUE) 
                # 
                results2 <- gdata
                results3 <- mydat
                results <- multiResultClass(); results$result1 <- results1; results$result2 <- results2; results$result3 <- results3
                return(results)
                #############################################
          }
          
}

##########

if(ncores>1){parallel::stopCluster(cl)}

if(ret==TRUE){sim.data[[1]] <- lapply(analysis, '[[', 2)}
analysis   <- do.call("rbind",lapply(analysis, '[[', 1))

if(i.type!="natural"){store.results$psi <- c(analysis)}else{
  store.results$psi[store.results$a1=="natural"] <- analysis[,Ynodes]
  store.results$psi[store.results$a1=="observed"] <- sapply(subset(X,select=Ynodes),rmean)
  if(is.null(Lnodes)==FALSE){
  lblocks <- rep(NA,length(Ynodes))
  for(i in 1:(length(Ynodes))){
  blocks <- make.interval(which(colnames(X)%in%Ynodes))
  lblocks[i] <- sum(which(colnames(X)%in%Lnodes)%in%blocks[[i]])
  }
   if(is.wholenumber(length(Lnodes)/length(Ynodes)) | n.t==1){
   store.results[store.results$a1=="natural",paste0("L_",1:(length(Lnodes)/n.t))] <- matrix(analysis[,Lnodes],byrow=T,ncol=length(Lnodes)/n.t)
   store.results[store.results$a1=="observed",paste0("L_",1:(length(Lnodes)/n.t))] <- matrix(sapply(subset(X,select=Lnodes),rmean),byrow=T,ncol=length(Lnodes)/n.t)}else{
      if(length(lblocks)>1 & length(unique(lblocks[-1]))==1){
        store.results[store.results$a1=="natural" & store.results$time!=1,paste0("L_",1:max(lblocks))] <- matrix(analysis[,Lnodes],byrow=T,ncol=max(lblocks))
        store.results[store.results$a1=="observed"& store.results$time!=1,paste0("L_",1:max(lblocks))] <- matrix(sapply(subset(X,select=Lnodes),rmean),byrow=T,ncol=max(lblocks))
        if(verbose==TRUE){cat("Note: it is not entirely clear to which time points your Lnodes belong to. I guessed what could make sense, but please check.\n\n")}
        }else{
        if(max(lblocks,na.rm=T)==1){
          updat.index <- unlist(lapply(apply(t(outer(which(colnames(X)%in%Lnodes), which(colnames(X)%in%Ynodes), "<")),1,which),max))
          updt.Lnodes <- analysis[,Lnodes][updat.index]; updt.Lnodes2 <-  sapply(subset(X,select=Lnodes),rmean)[updat.index]
          store.results[store.results$a1=="natural",paste0("L_",1)] <- updt.Lnodes
          store.results[store.results$a1=="observed",paste0("L_",1)] <- updt.Lnodes2
          }else{
            if(verbose==TRUE){cat("Note: your number of Lnodes are not a multiple of your number of Ynodes.\n This is fine, but you seem to have more than one Lnode per time point, and maybe unequally distributed. \n It is difficult for me to guess which Lnodes belong 'together', and at which time point. \n Thus, no natural course values for your Lnodes are calculated. \n Use 'ret=TRUE' and calculate manually.\n")}
              }}
   }}
  store.results[store.results$a1=="difference",grep("psi",colnames(store.results)):ncol(store.results)] <-
  store.results[store.results$a1=="natural",grep("psi",colnames(store.results)):ncol(store.results)] -
  store.results[store.results$a1=="observed",grep("psi",colnames(store.results)):ncol(store.results)] 
  }
  
if(B>0){
 for(b in 1:B){
    analysis.b[[b]] <- do.call("rbind",lapply(boots[[b]], '[[', 1))
    obs.data[[b+1]] <- boots[[b]][[1]]$result3
    if(ret==TRUE){for(i in 1:length(int.index)){sim.data[[b+1]][[i]]<- boots[[b]][[i]]$result2}}
  }
  boot.failure <- lapply(analysis.b,is.character)
  if(sum(unlist(boot.failure))>0){boots <- boots[-c(1:B)[unlist(boot.failure)]];analysis.b <- analysis.b[-c(1:B)[unlist(boot.failure)]]
                                  obs.data <- obs.data[-c(2:(B+1))[unlist(boot.failure)]]; sim.data <- sim.data[-c(2:(B+1))[unlist(boot.failure)]]
      if(verbose==TRUE){cat(paste("Caution:",sum(unlist(boot.failure)),"bootstrap sample(s) were removed due to errors \n"))}   }
  newB <- B-sum(unlist(boot.failure))
  if(i.type!="natural"){store.results[,c("l95","u95")] <-  t(apply(matrix(unlist(analysis.b),ncol=newB),1,quantile,probs=c(0.025,0.975)))}else{
  sim.results <- matrix(unlist(analysis.b),ncol=newB,dimnames=list(colnames(analysis.b[[1]]),NULL))
  store.results[store.results$a1=="natural",c("l95","u95")]  <-  t(apply(subset(sim.results,subset=rownames(sim.results)%in%Ynodes),1,quantile,probs=c(0.025,0.975)))
  store.results[store.results$a1=="observed",c("l95","u95")] <-  t(apply(matrix(unlist(lapply(obs.data[-1],lrmean,ind=Ynodes)),ncol=newB),1,quantile,probs=c(0.025,0.975)))
  store.results[store.results$a1=="difference",c("l95","u95")] <-  t(apply(sim.results[rownames(sim.results)%in%Ynodes,]-matrix(unlist(lapply(obs.data[-1],lrmean,ind=Ynodes)),ncol=newB),1,quantile,probs=c(0.025,0.975)))
  if(is.null(Lnodes)==FALSE){
    if(is.wholenumber(length(Lnodes)/length(Ynodes)) | n.t==1){
    store.results[store.results$a1=="natural", paste(paste0("L_",1:(length(Lnodes)/n.t)),"l95",sep=":")] <- 
    matrix(apply(sim.results[rownames(sim.results)%in%Lnodes,],1,quantile,probs=c(0.025)),ncol=(length(Lnodes)/n.t),byrow=T)
    store.results[store.results$a1=="natural", paste(paste0("L_",1:(length(Lnodes)/n.t)),"u95",sep=":")] <- 
    matrix(apply(sim.results[rownames(sim.results)%in%Lnodes,],1,quantile,probs=c(0.975)),ncol=(length(Lnodes)/n.t),byrow=T)
    store.results[store.results$a1=="observed", c(paste(paste0("L_",1:(length(Lnodes)/n.t)),"l95",sep=":"))] <- 
    matrix(t(apply(matrix(unlist(lapply(obs.data[-1],lrmean,ind=Lnodes)),ncol=newB),1,quantile,probs=c(0.025))),ncol=(length(Lnodes)/n.t),byrow=TRUE)
    store.results[store.results$a1=="observed", c(paste(paste0("L_",1:(length(Lnodes)/n.t)),"u95",sep=":"))] <- 
    matrix(t(apply(matrix(unlist(lapply(obs.data[-1],lrmean,ind=Lnodes)),ncol=newB),1,quantile,probs=c(0.975))),ncol=(length(Lnodes)/n.t),byrow=TRUE)
    col.order <- c(1:grep("psi",colnames(store.results)),which(colnames(store.results)%in%c("l95","u95")))
    for(i in 1:(length(Lnodes)/n.t)){col.order <- c(col.order,sort(grep(paste0("L_",i),colnames(store.results)))) }
    store.results <- store.results[,col.order]
    store.results[store.results$a1=="difference", paste(paste0("L_",1:(length(Lnodes)/n.t)),"l95",sep=":")] <- 
    matrix(apply(sim.results[rownames(sim.results)%in%Lnodes,]-matrix(unlist(lapply(obs.data[-1],lrmean,ind=Lnodes)),ncol=newB),1,quantile,probs=c(0.025)),ncol=(length(Lnodes)/n.t),byrow=T)
    store.results[store.results$a1=="difference", paste(paste0("L_",1:(length(Lnodes)/n.t)),"u95",sep=":")] <- 
    matrix(apply(sim.results[rownames(sim.results)%in%Lnodes,]-matrix(unlist(lapply(obs.data[-1],lrmean,ind=Lnodes)),ncol=newB),1,quantile,probs=c(0.975)),ncol=(length(Lnodes)/n.t),byrow=T)
    }else{
      if(length(lblocks)>1 & length(unique(lblocks[-1]))==1){
      store.results[store.results$a1=="natural" & store.results$time!=1, paste(paste0("L_",1:max(lblocks)),"l95",sep=":")] <- 
      matrix(apply(sim.results[rownames(sim.results)%in%Lnodes,],1,quantile,probs=c(0.025)),ncol=max(lblocks),byrow=T)
      store.results[store.results$a1=="natural" & store.results$time!=1, paste(paste0("L_",1:max(lblocks)),"u95",sep=":")] <- 
      matrix(apply(sim.results[rownames(sim.results)%in%Lnodes,],1,quantile,probs=c(0.975)),ncol=max(lblocks),byrow=T)
      store.results[store.results$a1=="observed" & store.results$time!=1, c(paste(paste0("L_",1:max(lblocks)),"l95",sep=":"))] <- 
      matrix(t(apply(matrix(unlist(lapply(obs.data[-1],lrmean,ind=Lnodes)),ncol=newB),1,quantile,probs=c(0.025))),ncol=max(lblocks),byrow=TRUE)
      store.results[store.results$a1=="observed" & store.results$time!=1, c(paste(paste0("L_",1:max(lblocks)),"u95",sep=":"))] <- 
      matrix(t(apply(matrix(unlist(lapply(obs.data[-1],lrmean,ind=Lnodes)),ncol=newB),1,quantile,probs=c(0.975))),ncol=max(lblocks),byrow=TRUE)
      col.order <- c(1:grep("psi",colnames(store.results)),which(colnames(store.results)%in%c("l95","u95")))
      for(i in 1:max(lblocks)){col.order <- c(col.order,sort(grep(paste0("L_",i),colnames(store.results)))) }
      store.results <- store.results[,col.order]
      store.results[store.results$a1=="difference" & store.results$time!=1, paste(paste0("L_",1:max(lblocks)),"l95",sep=":")] <- 
      matrix(apply(sim.results[rownames(sim.results)%in%Lnodes,]-matrix(unlist(lapply(obs.data[-1],lrmean,ind=Lnodes)),ncol=newB),1,quantile,probs=c(0.025)),ncol=max(lblocks),byrow=T)
      store.results[store.results$a1=="difference" & store.results$time!=1, paste(paste0("L_",1:max(lblocks)),"u95",sep=":")] <- 
      matrix(apply(sim.results[rownames(sim.results)%in%Lnodes,]-matrix(unlist(lapply(obs.data[-1],lrmean,ind=Lnodes)),ncol=newB),1,quantile,probs=c(0.975)),ncol=max(lblocks),byrow=T)
        }else{
          if(max(lblocks,na.rm=T)==1){
          store.results[store.results$a1=="natural", paste(paste0("L_",1),"l95",sep=":")] <- 
          matrix(apply(sim.results[rownames(sim.results)%in%Lnodes,],1,quantile,probs=c(0.025))[updat.index],ncol=1,byrow=T)
          store.results[store.results$a1=="natural", paste(paste0("L_",1),"u95",sep=":")] <- 
          matrix(apply(sim.results[rownames(sim.results)%in%Lnodes,],1,quantile,probs=c(0.975))[updat.index],ncol=1,byrow=T)
          store.results[store.results$a1=="observed", c(paste(paste0("L_",1),"l95",sep=":"))] <- 
          matrix(t(apply(matrix(unlist(lapply(obs.data[-1],lrmean,ind=Lnodes)),ncol=newB),1,quantile,probs=c(0.025))),1,byrow=TRUE)[,updat.index]
          store.results[store.results$a1=="observed", c(paste(paste0("L_",1),"u95",sep=":"))] <- 
          matrix(t(apply(matrix(unlist(lapply(obs.data[-1],lrmean,ind=Lnodes)),ncol=newB),1,quantile,probs=c(0.975))),1,byrow=TRUE)[,updat.index]
          col.order <- c(1:grep("psi",colnames(store.results)),which(colnames(store.results)%in%c("l95","u95")))
          col.order <- c(col.order,sort(grep(paste0("L_",1),colnames(store.results))))
          store.results <- store.results[,col.order]
          store.results[store.results$a1=="difference", paste(paste0("L_",1),"l95",sep=":")] <- 
          matrix(apply(sim.results[rownames(sim.results)%in%Lnodes,]-matrix(unlist(lapply(obs.data[-1],lrmean,ind=Lnodes)),ncol=newB),1,quantile,probs=c(0.025)),ncol=1,byrow=T)[updat.index,]
          store.results[store.results$a1=="difference", paste(paste0("L_",1),"u95",sep=":")] <- 
          matrix(apply(sim.results[rownames(sim.results)%in%Lnodes,]-matrix(unlist(lapply(obs.data[-1],lrmean,ind=Lnodes)),ncol=newB),1,quantile,probs=c(0.975)),ncol=1,byrow=T)[updat.index,]
      }}
    }
  }}
}
if(ret==FALSE){sim.data <- obs.data <- NULL}
if(i.type!="natural"){obs.data<-NULL}

# Step 6: calculate support if desired
if(calc.support==TRUE){
if(n.int<6 & n.int>2 & verbose==TRUE){cat(paste("Note: you have specified only",n.int,"interventions. The reported support diagnostics may not be reliable here. \n"))}
# if(length(Ynodes)==length(Anodes)){
# store.results$crude_weights <- c(support$crude_weights)
# store.results$cond_weights <- c(support$cond_weights)}else{
# store.results$crude_weights <- c(support$crude_weights[,updat.index2])
# store.results$cond_weights <- c(support$cond_weights[,updat.index2])
# }
diagn <- list(crude_support=support$crude_support,conditional_support=support$cond_support)
cn <- paste("a",1:n.a ,sep=""); if(i.type=="standard"){rn <- as.character(abar)}else{rn <- paste("Strategy",1:nrow(abar))}
diagn <- lapply(diagn,as.data.frame)
rownames(diagn$crude_support) <- rn; rownames(diagn$conditional_support) <- rn; colnames(diagn$crude_support) <- cn; colnames(diagn$conditional_support) <- cn
}else{diagn<-NULL}

if(is.null(prog)==FALSE){write(matrix("finished calculations.\n"),file=paste(prog,"/progress.txt",sep=""),append=TRUE)}

# Step 7: return results
res= list(results=store.results, 
          diagnostics=diagn, 
          simulated.data=sim.data, observed.data=obs.data,
          setup=list(i.type = i.type, n.t=n.t, B=B, fams=model.families, measure="default",
          Ynodes = Ynodes, Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes, abar=abar, support=calc.support, survival=survivalY)
          )

class(res) <- "gformula"
res

#
}