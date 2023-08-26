### NOTE ###
#create new S3 class (building block of everything):
#S3 class to capture the structural break model, will have the following
#properties:
#   numbr: number of breaks (not maximum, but a specify number of breaks); 
#other S3 class can inherit a different meaning of numbr later 
#which is different from the maximum number of break you want
# numbr is numeric
#   mform: special type of class formula, which can understand the formula syntax
# dedicated to mbreaks package. (More detailed on mformula.R)
###


#The model should be:
#Input is a formula or a data


#initialization of class sbm (only care about estimation of structural break model given
# number of breaks, not concerning with break selection. Break selection model should be
# subclass inherits many of the methods below from this block model/class)
# this class will not use trimming (no relationship to testing?) -> to simplify, this one
# might not have testing functions available for number of break dates. to test for 
# number of breaks, we will have second class that is capable of this
# under control, we have: 
# robust; prewhit; hetdat; hetvar; hetomega; hetq?
#' @export
sbm <- function(mform, data, m=2, h=10, control=NULL,  ...){
  
  
  #check validity of formula in structural break model
  sbm.fml = mfml.check(mform)
  #extract necessary information from processed formula
  zvars = sbm.fml$zvars
  xvars = sbm.fml$xvars
  yvar = sbm.fml$yvar
  intercept = sbm.fml$intercept
  tform = sbm.fml$tform
  
  #default list of controls/options
  con = list('robust' = TRUE,'hetomega'=TRUE, 'hetdat' = TRUE, 'hetq'=TRUE,'hetvar'=TRUE,'prewhit'=FALSE)

  #check for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  
  #build data as matrix
  sbm.frame = sbm.build(sbm.fml,data) 
 
  y <- sbm.frame$y
  x <- sbm.frame$x
  z <- sbm.frame$z
  
  #compute other parameters to fit in arguments of aux functions
  p = dim(x)[2]
  q = dim(z)[2]
  
  #check validity of minimum segment length options &
  #check validity of number of breaks
  #check_segOpt(h,m,bigt) 
  
  #estimate the break date with given number of break
  if (p==0){
    datevec <- dating_purescSSR(y, z, m, h)$datevec 
  }else{
    datevec <- dating_partscSSR(y, z, x, m, h)$datevec 
  }
  est_m = m
  brdate <- as.matrix(datevec[,est_m])
  #estimate the sbm model parameters with estimated break date (for inference)
  sbm.est <- estim(y,z,x,est_m,brdate,con$robust,con$prewhit,con$hetomega,con$hetq,con$hetdat,con$hetvar)
  out <- c()
  
  
  #group data into one separate field model.sbm
  model.sbm <- c()
  model.sbm$y <- y
  model.sbm$x <- x
  model.sbm$z <- z
  attr(model.sbm, 'dep.var') <- colnames(y)
  attr(model.sbm, 'xvars') <- colnames(x)
  attr(model.sbm, 'zvars') <- colnames(z)
  #group description of the model into field description
  description <- c()
 
  tform.char <- as.character(tform)
  dform <- paste(tform.char[2],tform.char[1],tform.char[3])
  description$formula <- tform
  attr(description, 'formula') <- dform
  attr(description, 'MAXbreaks') <- m
  attr(description, 'min.segLength') <- h
  if (p==0){attr(description, 'type') <- 'Pure'}
  else{attr(description, 'type') <- 'Partial'}
 
 
  # break date estimates and CI stored as bound (for inference)
  breakdate <- c()
  breakdate$date <- brdate
  breakdate$bound <- sbm.est$bound
  attr(breakdate, 'numbreaks') <- est_m #this denotes differently
  #to reflect it is estimated by IC/test (for later more complex class)
  #could be more way of estimating confidence interval for break date later (for example bootstrap)
  attr(breakdate, 'procedure') <- 'Asymptotic'
  # regressors estimates and CI stored as stdev (for inference)
  coefficients <- c()
  sbm.est$coef <- as.vector(sbm.est$coef)
  coefficients$coef.z <- sbm.est$coef[1:((est_m+1)*q)]
  coefficients$SE.cor.z <- sbm.est$SE_correct[1:((est_m+1)*q)]
  if (p == 0){
    coefficients$coef.x <- NULL
    coefficients$SE.cor.x <- NULL}
  else{
    coefficients$coef.x <- sbm.est$coef[((est_m+1)*q+1):((est_m+1)*q+p)]
    coefficients$SE.cor.x <- sbm.est$SE_correct[((est_m+1)*q+1):((est_m+1)*q+p)]}
  
  attr(coefficients, 'xregs') <- colnames(x)
  attr(coefficients, 'zregs') <- colnames(z)
  attr(coefficients, 'numx') <- p
  attr(coefficients, 'numz') <- q
  #add backs all above fields
  out$description <- description
  out$model <- model.sbm
  out$breakdate <- breakdate
  out$coefficients <- coefficients
  
  #post-regression 
  out$deviance = sbm.est$SSR
  out$residuals = sbm.est$resid
  out$fitted.values = sbm.est$fitted
  # keep tracks of options used
  out$controls <- con
  
  class(out) <- 'sbm'
  #format the results
  out$table <- sbm.compile(out)
  out
}


#Lists of methods to do
#1) summary
#2) print (format as suggested)
#3) plot*** (important)
#4) coef (must be able to split out estimated coefs on z and x (optional))
#5) brdate (get the estimated break date out, in vector and other form for display)
#6) model.frame (need to rewrite this), which extracting model frame from a formula/fit
# (Unnecessary, can reuse base R)
#7) resid (return already computed residuals of the model)


#check object class
is.sbm <- function(x) inherits(x, 'sbm')

##generics declaration
summary <- function(x, ...){
  UseMethod('summary')
}
print <- function(x, ...){
  UseMethod('print')
}
plot <- function(x, ...){
  UseMethod('plot')
}
coef <- function(x, ...){
  UseMethod('coef')
}
residuals <- function(x,...){
  UseMethod('residuals')
}
fitted <- function(x,...){
  UseMethod('fitted')
}
brdate <- function(x,...){
  UseMethod('brdate',x)
}
getOption <- function(x,...){
  UseMethod('getOption',x)
}

#S3 methods:
#' @rawNamespace S3method(print,sbm)
#' export(print.sbm)
print.sbm <- function(x, ...){
  cat(paste(gettextf('\n%s structural break model with %d breaks',
              attr(x$description,'type'),attr(x$breakdate,'numbreaks')),'\n'))
  cat(paste('Model specification: ',attr(x$description,'formula'),"\n"))
  if (attr(x$breakdate,'numbreaks') < 1) {
   cat(paste('\nNo break estimated'))
  }else{
   print(x$table$date_tab,quote=FALSE)}
  cat("\n")
}

#' @export
summary.sbm <- function(x, ... ){
  cat(paste(gettextf('%s structural break model with %d estimated breaks',
               attr(x$description,'type'),attr(x$breakdate,'numbreaks')),'\n'))
  
  cat(paste('Model specification: ',attr(x$description,'formula'),"\n"))
  
  cat('\nEstimated date:\n')
  print(x$table$date_tab,quote=FALSE)
  
  cat('\nEstimated regime-specific coefficients:\n')
  print(x$table$zregs_tab,quote=FALSE)
  
  if(attr(x$coefficients, 'numx') == 0) {cat('\nNo full sample regressors\n')}
  else{
    cat('\nEstimated full-sample coefficients:\n\n')
    print(x$table$xregs_tab,quote='FALSE')}
  cat('\nMinimum SSR =',
    format(round(x$deviance,3),nsmall=3),'\n')
  invisible(x)
}

#' @export
brdate.sbm <- function(x,...) x$breakdate$date

#' @export
coef.sbm <- function(x,...){ 
  out <- c()
  numbreaks = attr(x$breakdate, 'numbreaks')
  p = attr(x$coefficients, 'numx')
  q = attr(x$coefficients, 'numz')
  cname.z = c()
  coef.z = matrix(0L,q,numbreaks+1)
  rname.z= attr(x$coefficients, 'zregs')
  for (i in 1:(numbreaks+1)){
    cname.z = cbind(cname.z,paste('Regime',i))
  }
  print(x$coefficients$coef.z)
  for (j in 1:q){
      coef.z[j,] = x$coefficients$coef.z[((j-1)*(numbreaks+1)+1):(j*(numbreaks+1))]}
  colnames(coef.z) <- cname.z
  rownames(coef.z) <- rname.z
  coef.x = x$coefficients$coef.x
  names(coef.x) <- attr(x$coefficients, 'xregs')
  out$coefficients.z <- coef.z
  out$coefficients.x <- coef.x
  out
}

#' @importFrom graphics lines plot legend segments abline
#' @export
plot.sbm <- function(x,caption = NULL,xlab = NULL,ylab = NULL,bound=0.95,null.model = TRUE,start = NULL,...){
  est_m <- attr(x$breakdate, 'numbreaks')
  #plot x-axis based on provided start date (need to implement later).
  if(!is.null(start)) {}
  #comparison between structural break vs no break
  y <- unclass(x$model$y)
  bigt <- dim(y)[1]
  fity <- x$fitted.values
  x_t <- seq(1,dim(y)[1],1)
  range_y <- max(y)-min(y)
  if(is.null(ylab)){ylab <- colnames(y)}
  if(is.null(xlab)){xlab <- 'Time'}
  graphics::plot(x_t,y,type='l',col="black", xlab=xlab,ylab=ylab, 
       ylim=c(min(y)-range_y/10,max(y)+range_y/10),lty=1)
  #plot fitted values series for break model
  graphics::lines(x_t, fity,type='l', col="blue",lty=2)
  
  if(null.model){ #turn on plot for no break model
    zreg = x$model$z
    if(attr(x$coefficients,'numx') == 0){xreg=matrix(0,bigt,0)}else{
    xreg = x$model$x}
    fixreg = cbind(xreg,zreg)
    fixbeta = olsqr(y,fixreg)
    fity_fix = fixreg%*%fixbeta
    #plot fitted values series for null model
    graphics::lines(x_t, fity_fix,type='l', col="dark red",lty=2)
  }
  
  if(null.model){
    graphics::legend('topleft',legend=c("true",
                                        paste(est_m,'break(s)'),"0 break"),
           lty=c(1,2,2), col=c("black","blue","red"), ncol=1,bty='n')
  }else{
    #0,max(y)+range_y/10
    graphics::legend('topleft',legend=c("true",paste(est_m,'break(s)')),
           lty=c(1,2), col=c("black","blue"), ncol=1,bty='n')
  }
  #plot estimated dates + CIs based on bound option
  
  bigT <-length(y)
  for (i in 1:est_m){
    graphics::abline(v=x$breakdate$date[i,1],lty=2)}
  if (is.null(bound)){}
  else if (bound == 0.95){
    for (i in 1:est_m){
      lb_i = x$breakdate$bound[i,1]
      ub_i = x$breakdate$bound[i,2]
      if (lb_i < 0){lb_i = 0}
      if(ub_i>bigT){ ub_i=bigT}
      graphics::segments(lb_i,min(y)*(12+i/est_m)/10,ub_i,min(y)*(12+i/est_m)/10,lty=1,col='red')
    }
  }
  else if (bound == 0.9){
  for (i in 1:est_m){
    lb_i = x$breakdate$bound[i,3]
    ub_i = x$breakdate$bound[i,4]
    if (lb_i < 0){lb_i = 0}
    if(ub_i>bigT){ ub_i=bigT}
    graphics::segments(lb_i,min(y)*(12+i/est_m)/10,ub_i,min(y)*(12+i/est_m)/10,lty=1,col='red')
  }
  }
  invisible()
}


#' @export 
getOption.sbm <- function(x,...){
  #default list of controls/options
  con = x$controls
  cat(gettextf('robust = %s; hetdat = %s; hetvar = %s; hetomega = %s; hetq = %s',
           con$robust,con$hetdat,con$hetvar,con$hetomega,con$hetq))
}

#' @export
residuals.sbm <- function(x,...) x$residuals

#' @export
fitted.sbm <- function(x,...) x$fitted.values


# helper functions for class sbm
#' @importFrom stats model.frame
#' @export
sbm.build <- function(fml,data){
  #extract necessary information from processed formula
  zvars = fml$zvars
  xvars = fml$xvars
  yvar =  fml$yvar
  intercept = fml$intercept
  tform = fml$tform
  ##### process and transform data to appropriate y,z,x matrices ######
  #fix option = na.fail for handling missing values for now
  #since this is time series, missing observations cannot be tolerant
  if (missing(data)){
    sbm.frame = model.frame(tform, na.action = na.fail)}
  else{
    sbm.frame = data
  }
  
  #reconstruction of x,y,z in matrix form
  y_ind = match(yvar,colnames(sbm.frame))
  y = as.matrix(sbm.frame[,y_ind,drop=FALSE])
  bigt = dim(y)[1]
  #check if data has any observations
  if (bigt == 0) {stop('No observations in the data')}
  #check if any z regressors
  if (length(zvars)==0){ stop('No z regressors. Use `lm` instead')}
  
  #construct constants (as specified in formula)
  if (attr(intercept,'z.const') == 1){
    z.c = matrix(1,bigt,1)
    colnames(z.c)[1] <- '(Intercept)'
  }
  if (attr(intercept, 'x.const') == 1){
    x.c = matrix(1,bigt,1)
    colnames(x.c)[1] <- '(Intercept)'
  }
  
  #add additional regressors from frame
  z_ind = match(zvars,colnames(sbm.frame[,-1]))
  x_ind = match(xvars,colnames(sbm.frame[,-1]))
  
  #form regressor matrices for non-const vars (x is optional)
  z.v=as.matrix(sbm.frame[,z_ind[which(!is.na(z_ind))],drop=FALSE])
  x.v=as.matrix(sbm.frame[,x_ind[which(!is.na(x_ind))],drop=FALSE])
  
  if (attr(intercept,'z.const') == 1){
    if(identical(z_ind,integer(0))){z = z.c}
    else{z = cbind(z.c,z.v)}
  }else{z = z.v}
  if (attr(intercept,'x.const') == 1){
    if(identical(x_ind,integer(0))){x = x.c}
    else{x = cbind(x.c,x.v)}
  }else{x = x.v}
  list(y=y,z=z,x=x)
}

#' @export
sbm.compile <- function(x, digits = 3,...){
  ## format date estimation
  bound95 = c()
  bound90 = c()
  coln = c()
  numbreaks = attr(x$breakdate, 'numbreaks')
  for (i in 1:numbreaks){
    coln = c(coln, paste('Break',i,sep=''))
    bound95 = c(bound95,paste('(',x$breakdate$bound[i,1],',',x$breakdate$bound[i,2],')',sep=''))
    bound90 = c(bound90,paste('(',x$breakdate$bound[i,3],',',x$breakdate$bound[i,4],')',sep=''))
  }
  date_tab = data.frame(date = t(x$breakdate$date),stringsAsFactors = FALSE)
  date_tab = rbind(date_tab,bound95)
  date_tab = rbind(date_tab,bound90)
  colnames(date_tab) = coln
  rownames(date_tab) = c('Date','95% CI','90% CI')
  
  ## format z regressors coefficients (varying estimates over regimes)
  q = attr(x$coefficients, 'numz')
  rnameRS = attr(x$coefficients, 'zregs')
  cnameRS = c()
  coefRS = c()
  for (i in 1:(numbreaks+1)){
    cnameRS = cbind(cnameRS,paste('Regime',i))
  }
  for (j in 1:q){
    coefRSj = c()
    for (i in 1:(numbreaks+1)){
      coef_val = format(round(x$coefficients$coef.z[(j-1)*(numbreaks+1)+i],digits),nsmall=digits)
      coef_std = paste('(',format(round(x$coefficients$SE.cor.z[(j-1)*(numbreaks+1)+i],digits),nsmall=digits),')',sep='')
      coefRSj = cbind(coefRSj,paste(coef_val,coef_std,sep = ' '))
    }
    coefRS = rbind(coefRS,coefRSj)
  }
  
  rnameRSf = paste(rnameRS,'(SE)',sep=' ')
  coef_tabRS=data.frame(coef = coefRS)
  rownames(coef_tabRS) = rnameRSf
  colnames(coef_tabRS) = cnameRS
  
  ## format x regressors coefficients (constant estimates over regimes)
  p = attr(x$coefficients, 'numx')
  if(p == 0){coef_tabRW = NULL}
  else{
    rnameRW = attr(x$coefficients, 'xregs')
    cnameRW = 'All regimes'
    coefRW = c()
    for (j in 1:p){
      coef_val = format(round(x$coefficients$coef.x[j],digits),nsmall=digits)
      coef_std = paste('(',format(round(x$coefficients$SE.cor.x[j],digits),nsmall=digits),')',sep='')
      coefRW = cbind(coefRW,paste(coef_val,coef_std,sep=' '))}
    
    rnameRWf = paste(rnameRW,'(SE)',sep=' ')
    coef_tabRW=data.frame(coef = coefRW)
    colnames(coef_tabRW) = rnameRWf
    rownames(coef_tabRW) = cnameRW}
  
  table = list('date_tab' = date_tab,'zregs_tab' = coef_tabRS, 'xregs_tab' = coef_tabRW)
  return(table)
}


