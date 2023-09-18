model.formulas.update <- function(formulas,X,screening=screen.glmnet.cramer,with.s=FALSE,by=NA,...){
if(is.list(formulas)==FALSE){stop("Please provide an appropriate list (see help file)\n")}
selected.variables <- rep(list(NULL),length(formulas)); sf <- screening
for(i in 1:length(formulas)){
  if(is.null(formulas[[i]])==FALSE){
      selected.variables[[i]] <- rep(list(NA),length(formulas[[i]]))
         for(j in 1:length(formulas[[i]])){ 
         sv <- try(sf(X,formulas[[i]][j],...),silent=T)
         if(is.null(sv)==FALSE & !inherits(sv, "try-error")){selected.variables[[i]][[j]] <- sv}
         if(inherits(sv, "try-error")){selected.variables[[i]][[j]] <- NA}
  }}}
new.formulas <- formulas
for(i in 1:length(new.formulas)){
   if(is.null(new.formulas[[i]])==FALSE){
   for(j in 1:length(new.formulas[[i]])){
   covar <- "1"
   if(is.na(selected.variables[[i]][[j]][1])==FALSE){
   if(with.s==FALSE){covar <- paste(selected.variables[[i]][[j]],collapse=" + ")}else{
        X2 <- subset(X,select=selected.variables[[i]][[j]])
        cts.x <- apply(X2, 2, function(x) (length(unique(x)) > 10))
        if(sum(!cts.x) > 0 & sum(cts.x) > 0){
        covar <-paste(paste(paste("s(",
            colnames(X2[, cts.x, drop = FALSE]), ", by=", by,
            ")", sep = ""), collapse = " + "),
            "+", paste(colnames(X2[, !cts.x, drop = FALSE]),
                collapse = " + "))
        }
        if(sum(!cts.x) > 0 & sum(cts.x) == 0){
        covar <- paste(colnames(X2[, !cts.x, drop = FALSE]),
                collapse = " + ")
        }
        if(sum(!cts.x) == 0 & sum(cts.x) > 0){
        covar <- paste(paste(paste("s(",
            colnames(X2[, cts.x, drop = FALSE]), ", by=", by,
            ")", sep = ""), collapse = " + "))
         }}
   }
   if(substr(new.formulas[[i]][j],1,4)=="list"){new.formulas[[i]][j]<-paste(strsplit(strsplit(new.formulas[[i]][j],"))")[[1]][1],"\\(")[[1]][4],"~",strsplit(strsplit(new.formulas[[i]][j],",")[[1]][1],"~")[[1]][2])}
   new.formulas[[i]][j]  <- paste(all.vars(formula(new.formulas[[i]][j]))[1],"~",covar)
   if(substr(formulas[[i]][j],1,4)=="list"){
   gam.formula.2 <- paste("~",strsplit(new.formulas[[i]][j],"~")[[1]][2])
   gam.formula.3 <- paste0("as.numeric(as.character(",strsplit(new.formulas[[i]][j],"~")[[1]][1],"))",gam.formula.2)
   new.formulas[[i]][j] <- paste0("list(",paste(unlist(c(gam.formula.3,rep(list(gam.formula.2),length(strsplit(formulas[[i]][j],",")[[1]])-1))),collapse=","),")")
   }
   }
   }}
for(i in 1:length(new.formulas)){if(is.null(new.formulas[[i]])==FALSE){for(j in 1:length(new.formulas[[i]])){new.formulas[[i]][j]<-gsub(" ", "",new.formulas[[i]][j])}}}
return(new.formulas)
}