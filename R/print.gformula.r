print.gformula<-function(x, ...){
selcol <- c(which(colnames(x$results)=="psi"):dim(x$results)[2])
x$results[,selcol] <- apply(subset(x$results,select=selcol),2,round, digits=5)
print(x$results)
}
