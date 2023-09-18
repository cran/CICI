model.update <- function(gam.object,form){
if(gam.object$family$family=="multinom"){stop("Please update multinomial models manually.\n")}else{
return(update(gam.object,form,data=gam.object$data,family=gam.object$family$family))}
}

                    