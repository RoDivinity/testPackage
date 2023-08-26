#####. NOTE  ######
# This is rewritten R functions to decompose formula class
# dedicated only for mbreaks package (the syntax is required as below)
# The formula must have the following form:
# y ~ (?1)+z1+z2 | x1+x2
# The formula must have scalar dependent variable RHS, invoked by name y
# The LHS is more complicated:
# It could have either 1 or 2 parts, separated by operator |
# LHS before | is variable with coefficients changed across regimes, or z vars
# Implicitly, constant is included. If users want to impose unchanged mean 
# across regimes, must specify -1 as in usual formula
# LHS after | is optional. They specify regressors with unchanged coefficients
##########





###### checking validity of formula
# process formula when specifying mbreaks model
# Input: a R formula specific to mbreak models as explained in header NOTE
# Require lhs, and rhs
# RHS includes: z regressors (coefficients allowed to change across regimes, IMPLICITLY has constant) 
# and x regressors separated by | (coefficients are not changed across regimes, excluding constant)
# Output: a list containing: y name, z name and x name (if any)
###

#several checks are implemented to ensure conformity:
# check for duplicate regressors (avoid perfect collinearity)
# check for constant in z, x or none
# check for dependent variable


mfml.check <- function(mform) {
  envForm <- attr(mform,'.Environment')
  #first check: if expression for evaluation is a formula
  if(!inherits(mform,"formula")){stop('The expression for model is not a formula.')}
  
  #second check: at least y ~ z formula
  if (length(mform) == 2){ stop('Missing dependent variable y')}
  
  #split formula & store rhs/lhs
  rhs <- if(length(mform) > 2) mform[[3]] #store rhs
  lhs <- if(length(mform) > 2) mform[[2]] #store lhs
  
  if(is.null(rhs)) {stop('There is no z and x regressors')} #optional check 
  yvar <- lhs
  
  #third check: correct syntax in RHS (only 0 or 1 separator | is allowed)
  p.rhs <- mfml.check_operator(rhs)
  zvars.pr <- p.rhs$zvars
  xvars.pr <- p.rhs$xvars
  zvars.call <- p.rhs$zvars.call
  xvars.call <- p.rhs$xvars.call
  
  #fourth check: constant specification
  c.rhs <- mfml.check_const(zvars.pr,xvars.pr)
  zvars.final <- c.rhs$zvars
  xvars.final <- c.rhs$xvars
  intercept <- c.rhs$intercept
  
  #fifth check: duplicate regressors
  mfml.check_overlap(zvars.final,xvars.final)
  
  #sixth check: no intercepts in both x and z and no regressors specify (-1|-1)
  mfml.check_empty(zvars.final,xvars.final,intercept)
  
  if(is.null(xvars.call)){flagNoX <- TRUE}else{flagNoX <- FALSE}
  #reconstruct the formula via calls so we can invoke data frame later
  if (flagNoX){
    mform.final <- . ~ .
    mform.final[[2]] <- yvar
    mform.final[[3]] <- zvars.call
  }
  else{
    mform.final <- . ~ . + .
    mform.final[[3]][[2]] <- zvars.call
    mform.final[[3]][[3]] <- xvars.call
    mform.final[[2]] <- yvar
  }
  attr(mform.final,'.Environment') <- envForm
  list('yvar'=as.character(yvar),'zvars' = zvars.final,'xvars' = xvars.final,intercept=intercept,'tform' = mform.final)
}


mfml.check_operator <- function (rhs, op = '|'){
  #case 1: no x regressors => no | operator
  if(length(rhs) < 3 || rhs[[1]] != op){
    zvars.call <- rhs
    xvars.call <- NULL
    zvars <- mfml.get_vars(zvars.call)
    xvars <- list()
    message('No x regressors in formula')} #the formula contains only z regrssors
  #case 2: at least 1 | operator
  else {
    #split potential z and x regressors
    zvars.call = rhs[[2]] #store z regressors
    xvars.call = rhs[[3]] #store x regressors
    if(length(zvars.call)>1 && zvars.call[[1]] == op)
      {stop('Too many separator `|` in RHS regressors. You can only have x and z regressors')}
    if(length(xvars.call)>1 && xvars.call[[1]] == op)
      {stop('Too many separator `|` in RHS regressors. You can only have x and z regressors')}
    zvars <-mfml.get_vars(zvars.call)
    xvars <-mfml.get_vars(xvars.call)
  }
  list(zvars = zvars,xvars = xvars, zvars.call = zvars.call, xvars.call = xvars.call)
}



###
# Need to check overlapping x and z regressors (separate help function for easier
# maintainance)
mfml.check_overlap <-function(xregs,zregs){
  #check if any x regressors
  if (is.null(xregs)){
  }
  else{
    #get the type regressors with more variables
    xlen = length(xregs)
    zlen = length(zregs)
    if (xlen > zlen){
      lvars = xregs
      svars = zregs}
    else{
      lvars = zregs
      svars = xregs}
    if (any(as.vector(lvars)%in%as.vector(svars))){
      ov.flag = as.vector(lvars)%in%as.vector(svars)
      ov.vars = lvars[which(ov.flag==TRUE)]
      stop(gettextf('Variables %s are both z and x type regressors',
                    toString(ov.vars)))}
  }
}

#check if regression model has constant terms in z or x or none 
# To specify constant mean across regimes, 
# must explicitly specify 1 in x regressors
mfml.check_const <- function(zvars,xvars){
  zvars <- as.vector(zvars)
  xvars <- as.vector(xvars)
  intercept <- list()
  #case 1: specify -1 in zterms => no constant in z regressors
  if ('-1' %in% zvars){
    #remove '-1' in z terms & set property of z regressors
    zvars <- zvars[-which(zvars=='-1')]
    attr(intercept,'z.const') <- 0
    #case 1a: specify -1 in xterms => no constant in both z and x regressors
    if('-1' %in% xvars){
      #remove '-1' in x terms & set property of x regressors
      xvars <- xvars[-which(xvars=='-1')]
      attr(intercept,'x.const') <- 0}
    else{attr(intercept,'x.const') <- 1}}
  #case 2: implicitly include intercept in z regressors
  else {
    #constant is prioritize in z regressors
    attr(intercept,'z.const') <- 1
    attr(intercept,'x.const') <- 0 
  }
  list(zvars = zvars,xvars = xvars,intercept = intercept)
}



#check if formula is empty (only when no z and x regressors, 
# and users specify -1 in both x and z)
mfml.check_empty <- function(zvars,xvars,intercept){
  if(identical(zvars,character(0)) && identical(xvars,character(0)) ){
    if (attr(intercept,'z.const')==0 && attr(intercept,'x.const')==0){
      stop('Model is empty with no regressors')
    }
  }
}

#collect all variables in part of formula
mfml.get_vars <- function(vars.fml,sep = '+'){
  vars <- c()
  if (length(vars.fml) == 1) {
    vars = vars.fml}
  else{
    #extract each part in call from formula sequentially
    while(length(vars.fml)>1 && vars.fml[[1]] == sep){
      vars <- c(vars.fml[[3]],vars)
      vars.fml <- vars.fml[[2]]
    }
    #add last term on the right
    vars <- c(vars.fml, vars)
  }
  vars <- as.character(vars) #convert list of calls to characters
  vars <- unique(vars) #remove duplicate calls
}


# 
# #test case
# testf = y~x1+x2+x3+x4
# testf_1sep = y~x1+x2|x3+x4
# testf_2sep = y~x1|x2|x3+x4
# testf_3sep = y~x1|x2|x3|x4
# testf_dsep = y~x1&x2|x3+x4
# testf_dsep2 = y~x1|x2+x3&x4
# testf_dsep3 = y~x1|x2|x3&x4
# 
# # #debug code cont
# x1 = runif(10)
# x2 = runif(10)
# x3 = runif(10)
# x4 = runif(10)
# y = runif(10)
# testnew = sbm.check_formula(y~x1+x2+x3+x4|x2+x1+x3)
# testnew2 = sbm.check_formula(y~x1+x2|x3)
# zv = testnew$zvars
# xv = testnew$xvars
