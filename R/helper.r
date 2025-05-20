######################
# helper functions   #
######################
ncat              <- function(myvec){length(unique(myvec[is.na(myvec)==F]))}
is.wholenumber    <- function(myvec, tol = .Machine$double.eps^0.5){if(is.numeric(myvec)){abs(myvec - round(myvec)) < tol}else{rep(FALSE,length(myvec))}}
test.binary       <- function(myvec){if(ncat(myvec)==2){return(TRUE)}else{return(FALSE)}}
test.poisson      <- function(myvec){if(ncat(myvec)%in%seq(3,20) & all(is.wholenumber(myvec[is.na(myvec)==F])) & is.numeric(myvec[is.na(myvec)==F])){return(TRUE)}else{return(FALSE)}}
test.multinomial  <- function(myvec){if(ncat(myvec)%in%seq(3,9)  & is.factor(myvec[is.na(myvec)==F])){return(TRUE)}else{return(FALSE)}}
is.binary         <- function(myvec){all(myvec[is.na(myvec)==F]%in%c(0,1))}
   
assign.family     <- function(mymat,specify=NULL){
                                     if(is.null(specify)==TRUE){
                                     bin.ind  <- pois.ind <- mult.ind <- vector(mode = "logical", length = ncol(mymat))
                                     ncat.mat <- apply(mymat,2,ncat)
                                     for(i in 1:ncol(mymat)){bin.ind[i]  <- test.binary(mymat[,i])
                                                             pois.ind[i] <- test.poisson(mymat[,i])
                                                             mult.ind[i] <- test.multinomial(mymat[,i])}
                                     fam <- rep("gaussian",ncol(mymat))
                                     fam[bin.ind]  <-  "binomial"
                                     fam[pois.ind] <-  "poisson"
                                     fam[mult.ind] <-  paste0("multinom(K=",ncat.mat[mult.ind]-1,")")
                                     }else{fam<-specify}
                                     return(fam)
}

make.formula <- function(dat,approach,index,fam,surv=F,survnodes=NULL){
    if(approach[1]%in%c("GLM","GAM")){
    outc <- colnames(dat)[length(colnames(dat))]
    dat2 <- subset(dat,select=-length(colnames(dat)))
    if(surv==TRUE){
      dat2 <- dat2[,!(names(dat2) %in% survnodes),drop=F]
    }
    if(length(colnames(dat))==1){covar <- "1"}else{
      if(approach=="GLM"){covar <- paste(colnames(dat2),collapse="+")}
      if(approach=="GAM"){
        cts.x <- apply(dat2, 2, function(x) (length(unique(x)) > 10))
        if(sum(!cts.x) > 0 & sum(cts.x) > 0){
        covar <-paste(paste(paste("s(",
            colnames(dat2[, cts.x, drop = FALSE]),
            ")", sep = ""), collapse = "+"),
            "+", paste(colnames(dat2[, !cts.x, drop = FALSE]),
                collapse = "+"))
        }
        if(sum(!cts.x) > 0 & sum(cts.x) == 0){
        covar <- paste(colnames(dat2[, !cts.x, drop = FALSE]),
                collapse = "+")
        }
        if(sum(!cts.x) == 0 & sum(cts.x) > 0){
        covar <- paste(paste(paste("s(",
            colnames(dat2[, cts.x, drop = FALSE]),
            ")", sep = ""), collapse = "+"))
         }
                        }
      }
      gam.formula <- paste(outc,"~",covar); gam.formula <- as.formula(gam.formula)
      if(substr(fam,1,4)=="mult"){
      gam.formula.2 <- paste("~",as.character(gam.formula)[3])
      gam.formula.3 <- paste0("as.numeric(as.character(",as.character(gam.formula)[2],"))","~",as.character(gam.formula)[3])
      gam.formula <- lapply(c(gam.formula.3,rep(list(gam.formula.2),as.numeric(substr(fam,12,12))-1)),as.formula)  }
      }else{gam.formula<-eval(parse(text=approach[index]))}
      return(gam.formula)
}

calculate.support <- function(dat,A,intervention,projection=projection_linear,...){
 if(length(A)==1){size<-1; dat$placeholder <- dat[,A]; A <- c(A,"placeholder"); intervention <- cbind(intervention,intervention)}else{size<-2} # this fixes the problem when we only have one A column
 dat2 <- dat3 <-dat
 if(dim(dat[,A])[2]!=ncol(intervention)){stop("Number of columns of A and intervention don't match")}
 my.cuts <- matrix(NA,ncol=ncol(intervention),nrow=nrow(intervention)-1)
 for(i in 1:dim(my.cuts)[2]){my.cuts[,i] <- (intervention[,i][-length(intervention[,i])] + intervention[,i][-1])/2}  # definition of intervals (and epsilon indirectly)
 for(i in 1:length(A)){
 dat[,A][,i] <- cut(subset(dat,select=A)[,i],breaks=c(-Inf,my.cuts[,i],Inf))  # needs adaption if epsilon would be controlled directly
 }
 cut.intervals <- rep(list(NULL),length(A))
 support.model.names  <-   paste("sm_",as.vector(outer(A, c(seq(1:length(intervention[,1]))), paste, sep="_")),sep="")
 followed_A <-  rep(list(rep(list(NULL),length(A))),nrow(intervention))
 followed_A_mat <- rep(list(NA),nrow(intervention))
 followed_A_consec <- rep(list(NA),length(A))
 for(i in 1:length(cut.intervals)){cut.intervals[[i]] <- cbind(c(-Inf,my.cuts[,i],Inf)[-length(c(-Inf,my.cuts[,i],Inf))] , c(-Inf,my.cuts[,i],Inf)[-1] )}
 for(i in 1:length(A)){
  for(j in 1:nrow(intervention)){
  followed_A[[j]][[i]] <- as.numeric((dat2[,A][,i] > cut.intervals[[i]][j,1]) & (dat2[,A][,i] <  cut.intervals[[i]][j,2]))  # crude
 }}
 for(i in 1:length(followed_A)){followed_A_mat[[i]]<-matrix(unlist(followed_A[[i]]),ncol=length(followed_A[[i]]))} # crude 
 followed_A_consec <- followed_A_mat
 itt <- function(myvec){myvec[is.na(myvec)]<-0;for(i in 2:length(myvec)){if(myvec[i-1]==0){myvec[i]<-0}}
                        return(myvec)}
 itt2 <- function(mymatrix){t(apply(mymatrix,1,itt))}
 followed_A_consec <- lapply(followed_A_consec,itt2) # -> if you don't follow rule at time t, you can't follow in future
 include_time0 <- function(mymat){return(cbind(1,mymat))}
 include<- lapply(followed_A_consec,include_time0)   # at time 0, before intervention, everyone follows the rule 
                                                     # CRUCIAL: who to include in support model. Currently, only those who follow "rule"
 for(i in 1:length(A)){
  for(j in 1:nrow(intervention)){
  dat3[,A][,i]    <- as.numeric((dat2[,A][,i] > cut.intervals[[i]][j,1]) & (dat2[,A][,i] <  cut.intervals[[i]][j,2]))
  support.model   <- suppressWarnings(try(glm(as.formula(paste(A[i],"~.")),data=dat3[as.logical(include[[j]][,i]),1:(which(colnames(dat3)%in%A[i]))],family=binomial),silent=TRUE))  # support given FULL past, adapt as needed if required
  assign(paste("sm_",A[i],"_",c(seq(1:length(intervention[,1])))[j],sep=""),support.model)
  dat3[,A][,i]  <- dat2[,A][,i]
 }}
 include_list <- lapply(lapply(include,as.data.frame),as.list)  # list, to define subset which follows the "rule" of interest
 all.predictions  <- all.predictions.2 <- rep(list(rep(list(rep(0,nrow(dat))),length(A))),nrow(intervention))
  for(i in 1:length(A)){
  for(j in 1:nrow(intervention)){
  pred.model   <- try(get(paste("sm_",A[i],"_",c(seq(1:length(intervention[,1])))[j],sep="")),silent=TRUE) # TO DO: this model, or variable screening?
  if(!inherits(pred.model, "try-error")){all.predictions[[j]][[i]][as.logical(include_list[[j]][[i]])]<-suppressWarnings(try(predict(pred.model,type="response",newdata=pred.model$data),silent=TRUE))}
  #all.predictions[[j]][[i]][followed_A_consec[[j]][,i]==0]<- 0 #check? Should predictions for those that don't follow rule be zero? -> 13.2.21: ?
  }}
  mymean <- function(mymat){apply(mymat,2,mean)}
  crude_support <- matrix(unlist(lapply(followed_A_consec,mymean)),nrow=nrow(intervention),ncol=ncol(intervention),byrow=T) # CRUDE SUPPORT, i.e. ga
  exp.support <- function(vec){mean(vec)}
  mnz.list      <- function(lis){lapply(lis,exp.support)}
  cond_support  <- suppressWarnings(matrix(unlist(lapply(all.predictions,mnz.list)),nrow=nrow(intervention),ncol=ncol(intervention),byrow=T)) # suppress Warnings, as min(empty)=INf -> warning
  cond_support[is.na(cond_support) | is.infinite(cond_support)] <-0
  cond_support <- round(t(apply(cond_support,1,cumprod)),digits=6)
  crude_weights <- apply(crude_support,2,projection,...)
  cond_weights <- apply(cond_support,2,projection,...)
  if(size==1){crude_support <- crude_support[,1,drop=F];cond_support <- cond_support[,1,drop=F]}
  return(list(crude_weights=crude_weights,cond_weights=cond_weights,crude_support=crude_support,cond_support=cond_support, gal = all.predictions))
}

projection_linear <- function(x,c1=0.1,c2=0.1){
 if(c1<0 | c1>1){stop("c1 needs to be in [0, 1]")}
 if(c2<0 | c2>1){stop("c2 needs to be in [0, 1]")}
 w <- x; w[x>=c2] <- 1; w[x<c2] <- c1 + ((1-c1)/(c2))*x[x<c2]; w[x==0] <- c1
 return(w)
 }


require.package <- function(package, message = paste("loading required package (", 
    package, ") failed; please install", sep = "")){
    if (!requireNamespace(package, quietly = FALSE)) {
        stop(message, call. = FALSE)
    }
    invisible(TRUE)
}

screen.cramersv <- function(dat, form, nscreen=4, cts.num=10, ...){
    if(length(all.vars(formula(form))[-1])>1){
    dat <- na.omit(dat[,all.vars(formula(form))] )
    var_cont <- apply(dat, 2, function(x) (length(unique(x)) > cts.num))
    cutf <- function(x){cut(x, unique(quantile(x, prob = c(0, 0.2, 0.4, 0.6, 0.8, 1))),include.lowest=T)}
    if(any(var_cont)){dat[, var_cont] <- apply(dat[, var_cont, drop = FALSE], 2, cutf)}
    Y <- dat[,all.vars(formula(form))[1]]; X <- dat[,all.vars(formula(form))[-1]]
    calc_cram_v <- function(x_var, y_var) cramer(table(y_var, x_var))
    cramers_v <- apply(X, 2, calc_cram_v, y_var = Y)
    whichVariable <- colnames(X)[unname(rank(-cramers_v, ties.method = "random") <= nscreen)]}else{whichVariable <- NULL}
    return(whichVariable)
}

screen.glmnet.cramer <- function(dat, form, alpha = 1, pw=T, nfolds = 10, nlambda = 150, ...){
    require.package("glmnet")
    if(substr(form,1,4)=="list"){form<-paste(strsplit(strsplit(form,"))")[[1]][1],"\\(")[[1]][4],"~",strsplit(strsplit(form,",")[[1]][1],"~")[[1]][2])} 
    if(length(all.vars(formula(form)))>2){
    dat <- na.omit(dat[,all.vars(formula(form))] )
    Y <- as.vector(as.matrix(dat[,all.vars(formula(form))[1]])); X <- as.data.frame(dat[,all.vars(formula(form))[-1]])
    savedat<-dat; saveX<-X
    myfamily <- assign.family(data.frame(dat[,all.vars(formula(form))[1]])) 
    if(substr(myfamily,1,4)=="mult"){myfamily<-"multinomial"}
    # needed for factor variables later on
    if (ncol(X) > 26 * 27) stop("Too many variables for this screening algorithm.\n Contact Michael Schomaker for solution.")
    let <- c(letters, sort(do.call("paste0", expand.grid(letters, letters[1:26]))))
    names(X) <- let[1:ncol(X)]
    # factors are coded as dummies which are standardized in cv.glmnet()
    # intercept is not in model.matrix() because its already in cv.glmnet()
    is_fact_var <- sapply(X, is.factor)
    X <- try(model.matrix(~ -1 + ., data = X), silent = FALSE)
    successfulfit <- FALSE
    cvIndex <- rep(1:nfolds,trunc(nrow(X)/nfolds)+1)[1:nrow(X)]
    fitCV <- try(glmnet::cv.glmnet(
      x = X, y = Y, lambda = NULL, type.measure = "deviance",
      nfolds = nfolds, family = myfamily, alpha = alpha,
      nlambda = nlambda, keep = T, foldid=cvIndex
    ), silent = TRUE)
    # if no variable was selected, penalization might have been too strong, try log(lambda)
     if(!inherits(fitCV,"try-error")){if (all(fitCV$nzero == 0) | all(is.na(fitCV$nzero))) {
      fitCV <- try(glmnet::cv.glmnet(
        x = X, y = Y, lambda = log(fitCV$glmnet.fit$lambda + 1), type.measure = "deviance",
        nfolds = nfolds, family = myfamily, alpha = alpha, keep = T, foldid=cvIndex
      ), silent = TRUE)
    }}
    if(inherits(fitCV,"try-error")){successfulfit <- FALSE}else{successfulfit <- TRUE}
    whichVariable <- NULL
    if(successfulfit==TRUE){
    coefs <- coef(fitCV$glmnet.fit, s = fitCV$lambda.min)
    if(myfamily!="multinomial"){if(all(coefs[-1]==0)){whichVariable<-NULL}}else{
       min1<-function(vec){vec[-1]}
       if(all(unlist(lapply(coefs,min1))==0)){whichVariable<-NULL}}
    #}else{  
    if(myfamily!="multinomial"){var_nms <- coefs@Dimnames[[1]]}else{
                                var_nms <- coefs[[1]]@Dimnames[[1]]
                                allzero<-function(vec){all(vec==0)}  
                                sel<-apply(do.call("cbind",coefs)[-1,],1,allzero)
                                coefs <- coefs[[1]][-1];coefs[sel]<-0;coefs[!sel]<-0.1
                                coefs <- c(0,coefs)
                                }
    # Instead of Group Lasso:
    # If any level of a dummy coded factor is selected, the whole factor is selected
    if (any(is_fact_var)) {
      nms_fac <- names(which(is_fact_var))
      is_selected <- coefs[-1] != 0 # drop intercept
      # model.matrix adds numbers to dummy coded factors which we need to get rid of
      var_nms_sel <- gsub("[^::a-z::]", "", var_nms[-1][is_selected])
      sel_fac <- nms_fac[nms_fac %in% var_nms_sel]
      sel_numer <- var_nms_sel[!var_nms_sel %in% sel_fac]
      all_sel_vars <- c(sel_fac, sel_numer)
      whichVariable <- names(is_fact_var) %in% all_sel_vars
      } else {
      # metric variables only
        whichVariable <- coefs[-1] != 0
      }
      whichVariable <- colnames(saveX)[whichVariable]
      }
    #}  
    if(is.null(whichVariable)){whichVariable<-screen.cramersv(dat=savedat,form=form)
    if(pw==T){cat("Lasso failed and screening was based on Cramer's V (for ",form,")\n")}}
    }else{whichVariable<-NULL}
    return(whichVariable)
} 



censor <- function(vec,C.index){
   start.cens <- (which(vec==1)[which(vec==1)%in%C.index])[1]
   if(is.na(start.cens)==FALSE){if(start.cens!=length(vec)){vec[(start.cens+1):length(vec)] <- NA}}
   return(vec)
 }
 
adjust.sim.surv <- function(mat,Yn){  #improve for-loop
  mymin <- function(vec){if(length(vec)>0){return(min(vec))}else{return(NA)}}
  find.first <- function(vec){mymin(which(vec==1))}
  first.event <- apply(mat[,Yn,drop=F],1,find.first)
  censor.at <- Yn[first.event]
  position <- rep(NA,nrow(mat))
  for(i in 1:nrow(mat)){if(length(which(censor.at[i]==colnames(mat)))>0){position[i]<-which(censor.at[i]==colnames(mat))}}
  for(i in 1:nrow(mat)){
    if(is.na(position[i])==FALSE){
      if(position[i]+1<=ncol(mat)){
        mat[i,(position[i]+1):ncol(mat)]<-NA
        mat[i,Yn[(min(first.event[i]+1,length(Yn))):length(Yn)]]<-1
        } 
    }}
return(mat)
}


multiResultClass <- function(result1=NULL,result2=NULL)
{
  me <- list(
    result1 = result1,
    result2 = result2
  )
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}

 
multi.help <- function(vec){which(t(rmultinom(1,1,vec))==1)-1}
rmulti <- function(probmat){apply(probmat,1,multi.help)}
prop <- function(vec,categ=0){vec <- na.omit(vec);sum(vec==categ)/length(vec)}
rmean <- function(vec){
vec <- na.omit(vec)
if(!is.factor(vec)){mean(vec)}else{mean(as.numeric(as.character(vec)))}}
lrmean <- function(mmat,ind){sapply(subset(mmat,select=ind),rmean)}

factor.to.numeric <- function(vec,verb){nf <- levels(na.omit(vec))
nums <- 0:(length(nf)-1)
code <- data.frame(nf,nums)
recoding <- paste(apply(code,1,paste,collapse=" replaced by "),collapse=" ; ")
vec2 <- rep(NA,length(vec))
for(i in 1:length(nums)){vec2[vec==code$nf[i]]<-code$nums[i]}
if(verb==TRUE){cat(paste(recoding,"\n"))}
return(as.numeric(vec2))
}

recode_to_factor <- function(vec, verb = TRUE) {
  vec <- as.factor(vec)
  nf <- levels(na.omit(vec))           
  nums <- 0:(length(nf) - 1)            
  code <- data.frame(nf = nf, nums = nums)  
  recoding <- paste(apply(code, 1, function(row) paste(row, collapse = " replaced by ")),
                    collapse = " ; ")
  vec2 <- as.numeric(vec) - 1         
  if(verb){cat(paste(recoding, "\n")) }
  return(factor(vec2, levels = nums))   
}


binary.to.zeroone<-function(vec,verb){nf <- unique(na.omit(vec))
nums <- c(0,1)
code <- data.frame(nf,nums)
recoding <- paste(apply(code,1,paste,collapse=" replaced by "),collapse=" ; ")
vec2 <- rep(NA,length(vec))
for(i in 1:length(nums)){vec2[vec==code$nf[i]]<-code$nums[i]}
if(verb==TRUE){cat(paste(recoding,"\n"))}
return(vec2)
}

right.coding<-function(vec){
uv <- unique(na.omit(vec))
pc <- factor(0:(length(uv)-1))
if(identical(levels(uv),levels(pc))){return(TRUE)}else{return(FALSE)}
}

extract.families <- function(forms=NULL,fdata){
 if(is.null(forms)==FALSE){
 fams <- matrix(assign.family(fdata),nrow=1,dimnames=list(NULL,colnames(fdata)))
 outcomes <- model.fams <-  rep(NA,length(forms))
 for(i in 1:length(forms)){outcomes[i]<- strsplit(forms[i],"~")[[1]][1]}
 for(i in 1:length(outcomes)){if(nchar(outcomes[i])>4){if(substr(outcomes[i],1,4)=="list"){outcomes[i]<-strsplit(strsplit(outcomes[i],"\\(")[[1]][4],"\\)")[[1]][1]}}}
 #outcomes <- gsub('.{1}$', '', outcomes)
 for(i in 1:length(forms)){model.fams[i] <- fams[which(colnames(fams)%in%outcomes[i])]}
 return(model.fams)
 }else{return(NULL)}
}

missing.data <- function(ml){any(lapply(ml,is.matrix)==FALSE)}

cramer<-function(mt){
    allterms <- matrix(NA,nrow=nrow(mt),ncol=ncol(mt))
    n <- sum(mt)
     for(i in 1:nrow(mt)){
      for(j in 1:ncol(mt)){
        allterms[i,j]<- ((mt[i,j]- ((sum(mt[i,])*sum(mt[,j]))/(n)) )^2)/(((sum(mt[i,])*sum(mt[,j]))/(n)))
     }}
     cramer <- sqrt(sum(allterms)/(n*(min(nrow(mt),ncol(mt))-1)))
     cramer
}

make.interval <- function(vec){
ni <- length(vec)
intvs <- rep(list(NA),ni)
nv <- c(0,vec)
for(i in 1:ni){intvs[[i]] <- seq(nv[i],nv[i+1])}
intvs
}

.onAttach <- function(libname = find.package("CICI"), pkgname = "CICI") {
packageStartupMessage("The manual for this package will be available soon. \n Type ?CICI for a first overview.")
}

mi.inference <- function (est, std.err, confidence = 0.95){
    qstar <- est[[1]]
    for (i in 2:length(est)) {
        qstar <- cbind(qstar, est[[i]])
    }
    qbar <- apply(qstar, 1, mean)
    u <- std.err[[1]]
    for (i in 2:length(std.err)) {
        u <- cbind(u, std.err[[i]])
    }
    u <- u^2
    ubar <- apply(u, 1, mean)
    bm <- apply(qstar, 1, var)
    m <- dim(qstar)[2]
    tm <- ubar + ((1 + (1/m)) * bm)
    rem <- (1 + (1/m)) * bm/ubar
    nu <- (m - 1) * (1 + (1/rem))^2
    alpha <- 1 - (1 - confidence)/2
    low <- qbar - qt(alpha, nu) * sqrt(tm)
    up <- qbar + qt(alpha, nu) * sqrt(tm)
    result <- list(est = qbar, std.err = sqrt(tm), df = nu, 
        lower = low, upper = up, r = rem)
    result
}

correct.models <- list(
  "L"=c("adherence.1 ~ comorbidity.0 + efv.0",
        "weight.1 ~ sex + log_age",
        "comorbidity.1 ~ log_age + weight.0 + comorbidity.0",
        "list(as.numeric(as.character(dose.1)) ~  I(sqrt(weight.1)) + dose.0,  ~   I(sqrt(weight.1)) + dose.0, ~ I(sqrt(weight.1)) + dose.0)",
        "adherence.2 ~ comorbidity.1 + adherence.1 + efv.1",
        "weight.2 ~ weight.1 + comorbidity.1",
        "comorbidity.2 ~ log_age + weight.1 + comorbidity.1",
        "list(as.numeric(as.character(dose.2)) ~ I(sqrt(weight.2)) + dose.1, ~ I(sqrt(weight.2)) + dose.1, ~  I(sqrt(weight.2)) + dose.1)",
        "adherence.3 ~ comorbidity.2 + adherence.2 + efv.2",
        "weight.3 ~ weight.2 + comorbidity.2",
        "comorbidity.3 ~ log_age + weight.2 + comorbidity.2",
        "list(as.numeric(as.character(dose.3)) ~ I(sqrt(weight.3)) + dose.2, ~ I(sqrt(weight.3)) + dose.2, ~ I(sqrt(weight.3)) + dose.2)",
        "adherence.4 ~ comorbidity.3 + adherence.3 + efv.3",
        "weight.4 ~ weight.3 + comorbidity.3",
        "comorbidity.4 ~ log_age + weight.2 + comorbidity.2",
        "list(as.numeric(as.character(dose.4)) ~ I(sqrt(weight.4)) + dose.3, ~ I(sqrt(weight.4)) + dose.3, ~ I(sqrt(weight.4)) + dose.3)"
  ),
  "A"=c("efv.0 ~ log_age + metabolic*dose.0",
        "efv.1 ~ log_age + dose.0 + metabolic*adherence.1",
        "efv.2 ~ log_age + dose.0 + metabolic*adherence.2*dose.1",
        "efv.3 ~ log_age + dose.0 + metabolic*adherence.3*dose.2",
        "efv.4 ~ log_age + dose.0 + metabolic*adherence.4*dose.3"
  ),
  "Y"=c("VL.0 ~ I(sqrt(efv.0))",
        "VL.1 ~ I(sqrt(efv.1)) + comorbidity.0",
        "VL.2 ~ I(sqrt(efv.2)) + comorbidity.1",
        "VL.3 ~ I(sqrt(efv.3)) + comorbidity.2",
        "VL.4 ~ I(sqrt(efv.4)) + comorbidity.3"
  )
) 





    

