plot.gformula <- function(x, msm.method=c("line","loess","gam","none"),
                             CI=FALSE, time.points=NULL,
                             cols=NULL, weight=NULL, survival=FALSE,
                             variable="psi", difference=FALSE, ...){

a1<-psi<-l95<-u95<-relvar<-Strategy<-sel.l95<-sel.u95<-NULL
gobject <- x
msm.method  <- match.arg(msm.method)
if(msm.method=="none" & is.null(weight)==FALSE){stop("For weighted MSM estimation, `msm.method' cannot be `none'")}
if(gobject$setup$i.type=="natural"){if(length(unique(gobject$results$time))==1){stop("Natural course scenario for 1 time point can not be plotted. \n  Simply look at the results table.")}}
if(CI==TRUE & gobject$setup$B==0){CI<-FALSE;cat("No confidence intervals printed because B=0")}
if(survival==T & gobject$setup$survival==F){survival<-F; cat("Note: survival curves can only displayed for survival setups. Thus survival is set back to FALSE." )}
if(survival==T & gobject$setup$n.t==1){survival<-F; cat("Note: survival curves can only displayed for >1 time points." )}
if(survival==T & length(time.points)==1){stop("Survival curves can only displayed for >1 time points." )}
if(survival==T & msm.method!="none"){cat("Note: no step functions provided. \n  Use msm.method='none' for plot without connecting lines. \n")}

if(all(gobject$setup$abar%in%c(0,1))){bin.int=TRUE}else{bin.int <- FALSE}
n.t <- gobject$setup$n.t
if(bin.int==TRUE & n.t==1){gobject$setup$i.type<-"standard"; msm.method<-"none"}
if(bin.int==TRUE & n.t>1){gobject$setup$i.type<-"custom"}
if(survival==T){gobject$setup$i.type<-"custom"}
  
results <- gobject$results
if(is.null(time.points)==FALSE){results<-results[results$time%in%time.points,]}
mycolors <- c("black", "orangered3","dodgerblue4", "springgreen3","gold","greenyellow","purple",sample(rainbow(25)))
#if(is.null(weight)==FALSE){if(weight=="crude"){weight<-results$crude_weights}else{if(weight=="cond"){weight<-results$cond_weights}}}

#####
if(gobject$setup$i.type=="standard"){
if(is.null(cols)==FALSE){if(length(cols)!=length(unique(results$time))){stop(paste("Provide",length(unique(results$time)),"colours under `cols' (not more, not less - or use default colours)"))}}
if(is.null(cols)){cols<-mycolors[1:length(unique(results$time))]}

gg1 <- ggplot(results, aes(x=a1,y=psi,col=as.factor(time)))  + geom_point(size=1.75) + theme_bw()
if(msm.method=="none"){gg2<-gg1}
if(msm.method=="line" & bin.int==FALSE){gg2<-gg1+ geom_line(linetype=3, linewidth=1.05)}
if(msm.method=="loess" | msm.method=="gam"){mymethod=paste("(estimated with ",msm.method,")",sep="")
                        gg2<- gg1+ geom_smooth(method = msm.method, se = FALSE, aes(weight=weight))}else{mymethod=NULL}
if(CI==TRUE & bin.int==FALSE){gg2 <- gg2 + geom_ribbon(aes(ymin = l95, ymax = u95, fill=as.factor(time)), alpha=0.1,show.legend = FALSE)}
if(CI==TRUE & bin.int==TRUE & n.t==1){gg2 <- gg2 + geom_pointrange(aes(ymin = l95, ymax = u95))}
gg3 <- gg2 + scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
       scale_x_continuous("Intervention", breaks=round(c(seq(min(results$a1),max(results$a1),length.out=length(results$a1[results$time==min(results$time)]))),digits=1) ) +
       scale_y_continuous(expression(psi)) +
       guides(fill = guide_legend(keywidth = 2, keyheight = 2, title="Time")) +
       ggtitle(paste("Dose-response curves",mymethod))  +
       theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=13),axis.title.y = element_text(size=13, angle = 90),
            axis.text.y = element_text(size=13), legend.text =  element_text(size=13), legend.title = element_text(size=13, face = "bold", hjust = 0),legend.position =   "bottom") +
       guides(col=guide_legend(title="Time"))
}
#####

#####
if(gobject$setup$i.type=="natural"){   
if(is.null(weight)==FALSE){warning("It is not meaningful to specify weights for the natural course scenario. \n Weights are being ignored.")}
if(variable=="psi"){label<-expression(psi)}else{label<-variable}
#
  if(difference==FALSE){
results <- results[results$a1!="difference",]
if(is.null(cols)==FALSE){if(length(cols)!=length(results$time[results$time==1])){stop(paste("Provide",length(results$time[results$time==1]),"colours under `cols' (not more, not less - or use default colours)"))}}
if(is.null(cols)){cols<-mycolors[1:length(results$time[results$time==1])]}
results$Strategy <- rep(c("natural \n intervention", "observed \n data"),length(unique(results$time)))
colnames(results)[colnames(results)==variable] <- "relvar"
if(variable=="psi"){colnames(results)[colnames(results)%in%c("l95","u95")]<-c("sel.l95","sel.u95")}else{
                    colnames(results)[grep(paste0(variable,":"),colnames(results))]<-c("sel.l95","sel.u95")}

gg1 <- ggplot(results, aes(x=time,y=relvar,col=as.factor(Strategy)))  + geom_point(size=1.75, position = position_dodge(width = 0+0.25*(CI==T & msm.method=="none"))) + theme_bw()
if(msm.method=="none"){gg2<-gg1}
if(msm.method!="none"){gg2<-gg1+ geom_line(linetype=3, linewidth=1.05)}
if(msm.method=="loess" | msm.method=="gam"){cat("Note: MSM methods 'loess' and 'gam' not supported for natural interventions. \n Too few observations available for modeling (typically). \n msm.method='line' is used. \n")}
if(CI==TRUE & msm.method!="none"){gg2 <- gg2 + geom_ribbon(aes(ymin = sel.l95, ymax = sel.u95, fill=as.factor(Strategy)), alpha=0.1,show.legend = FALSE)}
if(CI==TRUE & msm.method=="none"){gg2 <- gg2 + geom_pointrange(aes(ymin = sel.l95, ymax = sel.u95), linewidth=1.05, show.legend=FALSE, position = position_dodge(width = 0.25))}
gg3 <- gg2 + scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
       scale_x_continuous("Time",breaks=unique(results$time)) +
       scale_y_continuous(label) +
       ggtitle("Natural Course Scenario")  +
       theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=13),axis.title.y = element_text(size=13, angle = 90),
            axis.text.y = element_text(size=13), legend.text =  element_text(size=13), legend.title = element_text(size=13, face = "bold", hjust = 0),legend.position =   "right") +
       guides(col=guide_legend(title="Strategies"))
}else{
#
if(variable=="psi"){label<-expression(psi)}else{label<-expression(variable)}
colnames(results)[colnames(results)==variable] <- "relvar"
diff.results <- results[results$a1=="difference",]
if(variable=="psi"){colnames(diff.results)[colnames(diff.results)%in%c("l95","u95")]<-c("sel.l95","sel.u95")}else{
                    colnames(diff.results)[grep(paste0(variable,":"),colnames(diff.results))]<-c("sel.l95","sel.u95")}
if(is.null(cols)==FALSE){cols <- "black"}
gg1 <- ggplot(diff.results, aes(x=time,y=relvar))  + geom_point(size=1.75, position = position_dodge(width = 0+0.3*(CI==T & msm.method=="none"))) + theme_bw()
if(msm.method=="none"){gg2<-gg1}
if(msm.method!="none"){gg2<-gg1+ geom_line(linetype=3, linewidth=1.05)}
if(msm.method=="loess" | msm.method=="gam"){cat("Note: MSM methods 'loess' and 'gam' not supported for natural interventions. \n Too few observations available for modeling (typically). \n MSM method='line' is used. \n")}
if(CI==TRUE & msm.method!="none"){gg2 <- gg2 + geom_ribbon(aes(ymin = sel.l95, ymax = sel.u95), alpha=0.1,show.legend = FALSE)}
if(CI==TRUE & msm.method=="none"){gg2 <- gg2 + geom_pointrange(aes(ymin = sel.l95, ymax = sel.u95), linewidth=1.05, show.legend=FALSE, position = position_dodge(width = 0.25))}
gg3 <- gg2 + scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
       scale_x_continuous("Time",breaks=unique(results$time)) +
       scale_y_continuous("Difference") + geom_hline(yintercept = 0) +
       ggtitle("Difference between observed data \n and estimates under natural intervention")  +
       theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=13),axis.title.y = element_text(size=13, angle = 90),
            axis.text.y = element_text(size=13), legend.text =  element_text(size=13), legend.title = element_text(size=13, face = "bold", hjust = 0),legend.position =   "right") +
       guides(col=guide_legend(title="Strategies"))
}
}
#####

#####
if(gobject$setup$i.type=="custom"){
if(is.null(weight)==FALSE){warning("Weights for custom estimands are being ignored. \n As for custom estimands `time' is on the x-axis *weighted* smoothing is not ideal if t is small. \n Plot your own weighted graphs as required. \n")}
if(is.null(cols)==FALSE){if(length(cols)!=length(results$time[results$time==1])){stop(paste("Provide",length(results$time[results$time==1]),"colours under `cols' (not more, not less - or use default colours)"))}}
if(is.null(cols)){cols<-mycolors[1:length(results$time[results$time==1])]}
results$Strategy <- rep(paste("Strategy:",apply(round(results[results$time==max(results$time),grep("a",colnames(results))],digits=2),1,paste,collapse="; ")),length(unique(results$time)))

gg1 <- ggplot(results, aes(x=time,y=psi,col=as.factor(Strategy), fill=as.factor(Strategy)))  + geom_point(size=1.75,show.legend = T, position = position_dodge(width = 0+0.3*(CI==T & msm.method=="none"))) + theme_bw()
if(msm.method=="none"){gg2<-gg1}
if(msm.method!="none"){gg2<-gg1+ geom_line(linetype=3, linewidth=1.05, show.legend = FALSE)}
if(msm.method=="loess" | msm.method=="gam"){msm.method<-"line"; cat("Note: MSM methods 'loess' and 'gam' not supported for custom interventions (and survival setting). \n Too few observations available for modeling (typically). \n msm.method='line' is used. \n")}
if(CI==TRUE & bin.int==FALSE & msm.method=="line"){gg2 <- gg2 + geom_ribbon(aes(ymin = l95, ymax = u95), alpha=0.1,show.legend = FALSE)}
if(CI==TRUE & (bin.int==TRUE | msm.method=="none")){gg2 <- gg2 + geom_pointrange(aes(ymin = l95, ymax = u95), show.legend=FALSE, position = position_dodge(width = 0.3))}
gg3 <- gg2 + scale_color_manual(values = cols) +  scale_fill_manual(values = cols) +
       scale_x_continuous("Time",breaks=unique(results$time)) +
       scale_y_continuous(expression(psi)) +
       ggtitle("Effect of Strategies over Time")  +
       theme(axis.title.x = element_text(size=13), axis.text.x = element_text(size=13),axis.title.y = element_text(size=13, angle = 90),
            axis.text.y = element_text(size=13), legend.text =  element_text(size=13), legend.title = element_text(size=13, face = "bold", hjust = 0),legend.position =   "right") +
       guides(col=guide_legend(title="Strategies"), fill="none")
}

suppressMessages(plot(gg3))
}