#####################################################
#####################################################
# Collection of functions to support Gibbs sampling
# originally from the spatial econometrics library
# adapted by MP Roeling from Matlab code into R
# december 2017, statements in functions are original
#####################################################

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
erfinv <- function(x) qnorm((x + 1)/2)/sqrt(2)

stdn_cdf = function(x){
  # cdf = .5*(ones(size(x))+erf(x/sqrt(2)))
  cdf = 0.5 * (rep(1, length(x)) + erf(x/sqrt(2)))
  return(cdf)
}

stdn_inv = function(x){  
  if(nargs != 1){
    stop('Wrong # of arguments to stdn_inv')
  }
  ninv = sqrt(2) * erfinv(2 * x - 1) 
  return(ninv)
}

norm_cdf = function(x, m, v){
# PURPOSE: computes the cumulative normal distribution 
#          for each component of x with mean m, variance v
  # Written by TT (Teresa.Twaroch@ci.tuwien.ac.at) on Jun 3, 1993
# Updated by KH (Kurt.Hornik@ci.tuwien.ac.at) on Oct 26, 1994
# Copyright Dept of Probability Theory and Statistics TU Wien
# Updated by James P. Lesage, jpl@jpl.econ.utoledo.edu 1/7/97

  if(is.data.frame(x)){
    r = nrow(as.data.frame(x))
    c = ncol(as.data.frame(x))
  } else {
    r = length(x)  
    c = 1
  }

  if (r*c == 0){
    stop('norm_cdf: x must not be empty')
  }
  
  cdf = rep(0, r)
  if (nargs() == 1){
    m = rep(0, r)
    v = rep(1, r)
    cdf[1:r] = stdn_cdf((x[1:r] - m[1:r]) * sqrt (v[1:r]))
  } else {
    cdf[1:r] = stdn_cdf((x[1:r] - m[1:r]) * sqrt (v[1:r]))
  }
  return(cdf)
}


normt_rnd = function(mu,sigma2,left,right){
  # PURPOSE: random draws from a normal truncated to (left,right) interval
  # ------------------------------------------------------
    # USAGE: y = normt_rnd(mu,sigma2,left,right)
  # where:   mu = mean (nobs x 1)
  #      sigma2 = variance (nobs x 1)
  #        left = left truncation points (nobs x 1),   999 for +Infinity
  #       right = right truncation points (nobs x 1), -999 for -Infinity
  # ------------------------------------------------------
    # RETURNS: y = (nobs x 1) vector
  # ------------------------------------------------------
    # NOTES: use y = normt_rnd(mu,sigma2,left,999)
  #        to produce a left-truncated draw
  #        use y = normt_rnd(mu_vector,sigma2_vector,-999*ones(nobs,1),right_vector)
  #        to produce a vector of right-truncated draws
  # ------------------------------------------------------
    # SEE ALSO: normlt_rnd (left truncated draws), normrt_rnd (right truncated)
  #
  
  # adopted from truncnorm3.m file by
  # Justin Tobias
  # James P. LeSage, Dept of Finance & Economics
  # Texas State Univeristy-San Marcos
  # 601 University Drive
  # San Marcos, TX 78666
  # jlesage@spatial-econometrics.com
  # last updated 10/2007
  #
  # For information on Bayesian Econometric Methods: 
    # Volume 7 of Econometric Exercises Series
  # Cambridge University Press by Gary Koop, Dale Poirer and Justin Tobias
  # see:
    # www.econ.iastate.edu/faculty/tobias/Bayesian_exercises.html
  
  
  if(nargs() != 4){
    stop('normt_rnd: wrong # of input arguments');
  }
  
  stderrs = sqrt(sigma2)
  
  points_left = stderrs == -999
  points_right = stderrs == 999
  
  a_term = norm_cdf((left-mu)/stderrs)
  a_term[points_left] = 0
  
  b_term = norm_cdf( (right-mu)/stderrs);
  b_term[points_right] = 1;
  
  uniforms = runif(length(mu))
  
  p = a_term + uniforms*(b_term - a_term)
  
  result = mu + stderrs*norm_inv(p)
  return(result)
}

norm_inv = function(x, m, v){
  # PURPOSE: computes the quantile (inverse of the CDF) 
  #          for each component of x with mean m, variance v
  #---------------------------------------------------
    # USAGE: invp = norm_inv(x,m,v)
  # where: x = variable vector (nx1)
  #        m = mean vector (default=0)
  #        v = variance vector (default=1)
  #---------------------------------------------------
    # RETURNS: invp (nx1) vector
  #---------------------------------------------------
    # SEE ALSO: norm_d, norm_rnd, norm_inv, norm_cdf
  #---------------------------------------------------
    
    # Written by KH (Kurt.Hornik@ci.tuwien.ac.at) on Oct 26, 1994
  # Copyright Dept of Probability Theory and Statistics TU Wien
  # Converted to MATLAB by JP LeSage, jpl@jpl.econ.utoledo.edu
  
  if(nargs > 3){
    stop('Wrong # of arguments to norm_inv');
  }
  
  r = dim(x1)[1]
  c = dim(x1)[2]
  s = r * c
  
  if (nargs() == 1){
    m = rep(0, 1)
    v = rep(1, 1)
  }
  
  x = matrix(x1, 1, s, byrow = TRUE)
  m = matrix(m, 1, s, byrow = TRUE)
  v = matrix(v, 1, s, byrow = TRUE)
  
  invp = rep(0, 1)
  
  invp = m + sqrt(v) * stdn_inv(x) # erfinv gaat fout

  invp = matrix(invp, r, c, byrow = TRUE)
  return(invp)
}



pr_like = function(b,y,x){
  # PURPOSE: evaluate probit log-likelihood
  #-----------------------------------------------------
  # USAGE:    like = pr_like(b,y,x,flag) 
  # where:     b = parameter vector (k x 1)
  #            y = dependent variable vector (n x 1)
  #            x = explanatory variables matrix (n x m)
  #-----------------------------------------------------
  # NOTE: this function returns a scalar
  #-----------------------------------------------------
  # SEE also: hessian, gradnt, gradt
  #-----------------------------------------------------
  # REFERENCES: Green, 1997 page 883
  #-----------------------------------------------------
  
  # written by:
  # James P. LeSage, Dept of Economics
  # University of Toledo
  # 2801 W. Bancroft St,
  # Toledo, OH 43606
  # jpl@jpl.econ.utoledo.edu
  
  # error check
  if (nargs() != 3){
    stop('wrong # of arguments to pr_like')
  }
  m = length(b)
  junk = 1

  i = rep(1, length(y))
  
  cdf = norm_cdf(x%*%b)
  
  tmp = cdf[cdf <= 0]
  n1 = length(tmp)
  if (n1 != 0){
    cdf[cdf <= 0] = 0.00001*rep(1, length(tmp))
  }
  
  tmp = cdf[cdf >= 1]
  n1 = length(tmp)
  if (n1 != 0){
    cdf[cdf <= 0] = 0.99999**rep(1, length(tmp))
  }
  
  out = c(y)*log(cdf)+(i-c(y))*log(i-cdf)
  like = sum(out)
  return(like)
}


gamm_rnd = function(n,a){
  # PURPOSE: a vector of random draws from the gamma distribution
  #---------------------------------------------------
    # USAGE: r = gamm_rnd(n,A)
  # where: n = the row size of the n x 1 vector drawn 
  #        a = a parameter such that the mean of the gamma = a
  #            and the variance of the gamma = a
  #        notes: x = gamm_rnd(n,a*0.5)*2,equals chisq a random deviate
  #        For different parameters, A,B use:
    #	B*gamm_rnd(n,A) to produce a vector of random deviates from the gamma
  #	distribution with shape parameter A and scale parameter B.  
  #   The distribution then has mean A*B and variance A*B^2.
  #---------------------------------------------------
    # RETURNS:
    #        r = an n x 1 vector of random numbers from
  #        the gamma(A) distribution      
  # --------------------------------------------------
    # SEE ALSO: gamm_inv, gamm_pdf, gamm_cdf
  #---------------------------------------------------
    
    # modified slightly by J. LeSage
  # to avoid an error in Matlab version 7.01
  
  #RGAMMA   Random numbers from the gamma distribution
  #
  #         x = rgamma(n,a)
  
  # GNU Public Licence Copyright (c) Anders Holtsberg 10-05-2000.
  
  # This consumes about a third of the execution time compared to 
  # the Mathworks function GAMRND in a third the number of
  # codelines. Yihaaa! (But it does not work with different parameters)
  #
  # The algorithm is a rejection method. The logarithm of the gamma 
  # variable is simulated by dominating it with a double exponential.
  # The proof is easy since the log density is convex!
    # 
  # Reference: There is no reference! Send me an email if you can't
  # figure it out.
  
  
  if(any(any(a<=0))){
    stop('Parameter a is wrong')
  }
  
  y0 = log(a)-1/sqrt(a)
  c = a - exp(y0)
  m = ceiling(n*(1.7 + 0.6*(min(min(a))<2)))
  
  y = log(runif(m))*sign(runif(m)-0.5)/c + log(a)
  f = a*y-exp(y) - (a*y0 - exp(y0))
  g = c*(abs((y0-log(a))) - abs(y-log(a)));
  reject = (log(runif(m)) + g) > f
  y = y[!reject] 
  x = rep(0, n)
  
  if(length(y) >= n){
     x = exp(y[1:n])
  }else{
     tmp = rgamma(n - length(y), a)
     x = c(exp(y), tmp)
  }
  return(x)
}

probit_g = function(y,x,ndraw,nomit,prior=NA,seed=2017){
  # PURPOSE: MCMC sampler for the Bayesian heteroscedastic Probit model  
  #          y = X B + E, E = N(0,V), 
  #          V = diag(v1,v2,...vn), r/vi = ID chi(r)/r, r = Gamma(m,k)
  #          B = N(c,T)
  # --------------------------------------------------------------
  # USAGE: results =  probit_g(y,x,ndraw,nomit,prior,start)
  # where: y = nobs x 1 independent variable vector
  #        x = nobs x nvar explanatory variables matrix
  #       ndraw = # of draws
  #       nomit = # of initial draws omitted for burn-in
  #       prior = a structure for prior information input
  #               prior.beta, prior means for beta,  c above (default=0)
  #               priov.bcov, prior beta covariance, T above (default=1e+12)
  #               prior.rval, r prior hyperparameter, default=4
  #               prior.m,    informative Gamma(m,k) prior on r
  #               prior.k,    informative Gamma(m,k) prior on r
  #       seed = (optional) string for the random number generator seed
  #              e.g., seed = num2str(1234);
  #---------------------------------------------------------------
  #----------------------------------------------------------------
  # NOTE: use either improper prior.rval 
  #       or informative Gamma prior.m, prior.k, not both of them
  #---------------------------------------------------------------
  # References: James H. Albert and Siddhartha Chib
  #             Bayesian Analysis of Binary and Polychotomous
  #             Response Data JASA, June 1993, pp. 669            
  #----------------------------------------------------------------
  
  # y = yc
  n = dim(x)[1] 
  k = dim(x)[2]
  
  # error checking on input
  # check for all 1's or all 0's
  if(length(unique(y)) != 2){
    stop("Not sure about your input dependent, 
          do they have only 2 categories without missing values?")
  }
  
  if(any(is.na(prior))){stop("no prior provided")}
  
  b0 = rep(1,k)
  sflag = 0
  
  if(nargs() == 6){ # user-supplied a seed 
    set.seed(seed)  # set seed 
    
    sflag = 1
    IN = rep(1, n)
    
    fields = prior
    nf = length(fields)
    mm = 0
    rval = 4        # rval = 4 is default
    nu = 0          # default diffuse prior for sige
    d0 = 0
    c = rep(0, k)
    T = diag(k)*1e+12
    
    for(i in 1:nf){
      if(names(fields[i]) == "rval"){
        rval = as.numeric(fields[i])
      }else if(names(fields[i]) == "m"){
        # if you provide a value of m you should also provide a value of k
        mm = as.numeric(fields[i])
        kk = as.numeric(fields[i+1])
        rval = gamm_rnd(1,1,mm,kk)    # initial value for rval
      }else if(names(fields[i]) == "beta"){
        c = as.numeric(fields[i])
      }else if(names(fields[i]) == "bcov"){
        T = as.numeric(fields[i])
      }else if(names(fields[i]) == "nu"){
        nu = as.numeric(fields[i])
      }else if(names(fields[i]) == "d0"){
        d0 = as.numeric(fields[i])
      }
    }
    
  } else if(nargs() == 5){ # probit maximum likelihood starting values
    b0 = rep(1,k) 
    V = rep(1,n)
    IN = rep(1, n)  # initial value for V 
    
    fields = prior
    nf = length(fields)
    mm = 0
    rval = 4        # rval = 4 is default
    nu = 0          # default diffuse prior for sige
    d0 = 0
    c = rep(0, k)
    T = diag(k)*1e+12
    
    for(i in 1:nf){
      if(names(fields[i]) == "rval"){
        rval = as.numeric(fields[i])
      }else if(names(fields[i]) == "m"){
        # if you provide a value of m you should also provide a value of k
        mm = as.numeric(fields[i])
        kk = as.numeric(fields[i+1])
        rval = gamm_rnd(1,1,mm,kk)    # initial value for rval
      }else if(names(fields[i]) == "beta"){
        c = as.numeric(fields[i])
      }else if(names(fields[i]) == "bcov"){
        T = as.numeric(fields[i])
      }else if(names(fields[i]) == "nu"){
        nu = as.numeric(fields[i])
      }else if(names(fields[i]) == "d0"){
        d0 = as.numeric(fields[i])
      }
    }
    
  } else if(nargs() == 4){ # use default prior
    
    mm = 0 
    rval = 4        # rval = 4 is default
    nu = 0          # default diffuse prior for sige
    d0 = 0
    c = rep(0, k)
    T = diag(k)*1e+12  
    V = rep(1, n)   # initial value for V 
    IN = rep(1, n)  
    
  } else {
    stop("Wrong # of arguments to probit_g")
  }
  
  checkk = length(c)
  if(checkk != k){
    stop('probit_g: prior means are wrong')
  }
  
  checkk = dim(T)[1]
  if(checkk != k){
    stop('probit_g: prior bcov are wrong')
  }
  
  Q = solve(T)
  Qpc = Q%*%c
  
  bsave = matrix(0,ndraw-nomit,k)    # allocate storage for results
  ymean = rep(0, n)
  rsave = rep(0, ndraw-nomit)
  vmean = rep(0, n)
  yhat = rep(0, n)
  
  sv = 1.0         
  yin = y             # save original y-values
  bhat = b0           # starting value for beta
  
  pb <- txtProgressBar(min = 0, max = ndraw, style = 3)
  
  for(i in 1:ndraw){
    xstar = x*c(sqrt(V))
    ystar = y*c(sqrt(V))
    
    lp=xstar%*%bhat 
    
    # sample Z from normal truncated right at 0 if yin = 0  
    # sample Z from normal truncated left at 0 if yin = 1
    # mean of truncated is the predicted value of yi (XiTB) stored in lp
    for(j in 1:length(yin)){
      if(yin[j] == 0){
        y[j] = rtruncnorm(1, a=-Inf, b = 0, mean = lp[j], sd = 1)
      } else if(yin[j] == 1){
        y[j] = rtruncnorm(1, a=0, b = Inf, mean = lp[j], sd = 1)
      }
    }
    
    # update beta
    xpxi = solve(t(xstar)%*%xstar + Q) 
    xpy = t(xstar)%*%ystar + Qpc
    bhat = xpxi%*%xpy
    for(j in 1:k){
      bhat[j] = bhat[j] +  rnorm(1, xpxi[j,j], 1)
    }
    
    # update V
    e = y - x%*%bhat  
    
    # i am assuming we need 100 values with df rval+1
    chiv = rchisq(n, rval+1)
    vi = ((e*e) + IN*rval) / chiv
    V = IN/vi
    if(mm != 0){
      rval = gamm_rnd(1, mm, kk)  # update rval
    }             
    
    if(i > nomit){ #if we are past burn-in, save the draws
      bsave[i-nomit,] = t(bhat)
      ymean = ymean + lp
      vmean = vmean + vi
      yhat = yhat + stdn_cdf(c(y))
      
      if(mm != 0){
        rsave[i-nomit] = rval
      }
    } # end of if i > nomit
    setTxtProgressBar(pb, i)
  } # End the sampling
  close(pb)
  
  vmean = vmean/(ndraw-nomit)
  ymean = ymean/(ndraw-nomit)
  yhat = yhat/(ndraw-nomit)
  
  bmean = colMeans(bsave)
  
  # compute McFadden R-squared
  tmp = yin[yin == 1] # find ones
  P = length(tmp)
  cnt0 = n-P
  cnt1 = P
  P = P/n             # proportion of 1's
  like0 = n*(P*log(P) + (1-P)*log(1-P))     # restricted likelihood
  like1 = pr_like(bmean,yin,x);            # unrestricted Likelihood
  r2mf = 1-(abs(like1)/abs(like0));         # McFadden pseudo-R2 
  # compute Estrella R-squared
  term0 = (2/n)*like0
  term1 = 1/(abs(like1)/abs(like0))^term0
  rsqr = 1-term1                            # Estrella R2
  
  # return results
  results <- NULL
  results$meth  = "probit_g";
  results$r2mf = r2mf;
  results$rsqr = rsqr;
  results$bdraw = bsave;
  results$pmean = c;
  results$pstd  = sqrt(diag(T));
  results$vmean = vmean;
  results$ymean = ymean;
  results$yhat = yhat;
  if(mm != 0){
    results$rdraw = rsave;
    results$m     = mm;
    results$k     = kk;
  }else{
    results$r     = rval;
    results$rdraw = rsave;
  }
  results$nobs  = n;
  results$nvar  = k;
  results$y     = yin;
  results$x     = x;
  results$ndraw = ndraw;
  results$nomit = nomit;
  results$pflag = "plevel";
  return(results)
}

probit_g_mp = function(y,x,ndraw,nomit,prior=NA,seed=2017){
  # PURPOSE: MCMC sampler for the Bayesian heteroscedastic Probit model  with imputation
  #          y = X B + E, E = N(0,V), 
  #          V = diag(v1,v2,...vn), r/vi = ID chi(r)/r, r = Gamma(m,k)
  #          B = N(c,T)
  # --------------------------------------------------------------
  # USAGE: results =  probit_g(y,x,ndraw,nomit,prior,start)
  # where: y = nobs x 1 independent variable vector
  #        x = nobs x nvar explanatory variables matrix
  #       ndraw = # of draws
  #       nomit = # of initial draws omitted for burn-in
  #       prior = a structure for prior information input
  #               prior.beta, prior means for beta,  c above (default=0)
  #               priov.bcov, prior beta covariance, T above (default=1e+12)
  #               prior.rval, r prior hyperparameter, default=4
  #               prior.m,    informative Gamma(m,k) prior on r
  #               prior.k,    informative Gamma(m,k) prior on r
  #       seed = (optional) string for the random number generator seed
  #              e.g., seed = num2str(1234);
  #---------------------------------------------------------------
  #----------------------------------------------------------------
  # NOTE: use either improper prior.rval 
  #       or informative Gamma prior.m, prior.k, not both of them
  #---------------------------------------------------------------
  # References: James H. Albert and Siddhartha Chib
  #             Bayesian Analysis of Binary and Polychotomous
  #             Response Data JASA, June 1993, pp. 669            
  #----------------------------------------------------------------
  
  # y = yc
  n = dim(x)[1] 
  k = dim(x)[2]
  
  # error checking on input
  # check for all 1's or all 0's
  if(length(unique(y)) != 2){
    stop("Not sure about your input dependent, 
         do they have only 2 categories without missing values?")
  }
  
  if(any(is.na(prior))){stop("no prior provided")}
  
  b0 = rep(1,k)
  sflag = 0
  
  if(nargs() == 6){ # user-supplied a seed 
    set.seed(seed)  # set seed 
    
    sflag = 1
    IN = rep(1, n)
    
    fields = prior
    nf = length(fields)
    mm = 0
    rval = 4        # rval = 4 is default
    nu = 0          # default diffuse prior for sige
    d0 = 0
    c = rep(0, k)
    T = diag(k)*1e+12
    
    for(i in 1:nf){
      if(names(fields[i]) == "rval"){
        rval = as.numeric(fields[i])
      }else if(names(fields[i]) == "m"){
        # if you provide a value of m you should also provide a value of k
        mm = as.numeric(fields[i])
        kk = as.numeric(fields[i+1])
        rval = gamm_rnd(1,1,mm,kk)    # initial value for rval
      }else if(names(fields[i]) == "beta"){
        c = as.numeric(fields[i])
      }else if(names(fields[i]) == "bcov"){
        T = as.numeric(fields[i])
      }else if(names(fields[i]) == "nu"){
        nu = as.numeric(fields[i])
      }else if(names(fields[i]) == "d0"){
        d0 = as.numeric(fields[i])
      }
    }
    
  } else if(nargs() == 5){ # probit maximum likelihood starting values
    b0 = rep(1,k) 
    V = rep(1,n)
    IN = rep(1, n)  # initial value for V 
    
    fields = prior
    nf = length(fields)
    mm = 0
    rval = 4        # rval = 4 is default
    nu = 0          # default diffuse prior for sige
    d0 = 0
    c = rep(0, k)
    T = diag(k)*1e+12
    
    for(i in 1:nf){
      if(names(fields[i]) == "rval"){
        rval = as.numeric(fields[i])
      }else if(names(fields[i]) == "m"){
        # if you provide a value of m you should also provide a value of k
        mm = as.numeric(fields[i])
        kk = as.numeric(fields[i+1])
        rval = gamm_rnd(1,1,mm,kk)    # initial value for rval
      }else if(names(fields[i]) == "beta"){
        c = as.numeric(fields[i])
      }else if(names(fields[i]) == "bcov"){
        T = as.numeric(fields[i])
      }else if(names(fields[i]) == "nu"){
        nu = as.numeric(fields[i])
      }else if(names(fields[i]) == "d0"){
        d0 = as.numeric(fields[i])
      }
    }
    
  } else if(nargs() == 4){ # use default prior
    
    mm = 0 
    rval = 4        # rval = 4 is default
    nu = 0          # default diffuse prior for sige
    d0 = 0
    c = rep(0, k)
    T = diag(k)*1e+12  
    V = rep(1, n)   # initial value for V 
    IN = rep(1, n)  
    
  } else {
    stop("Wrong # of arguments to probit_g")
  }
  
  checkk = length(c)
  if(checkk != k){
    stop('probit_g: prior means are wrong')
  }
  
  checkk = dim(T)[1]
  if(checkk != k){
    stop('probit_g: prior bcov are wrong')
  }
  
  Q = solve(T)
  Qpc = Q%*%c
  
  bsave = matrix(0,ndraw-nomit,k)    # allocate storage for results
  ymean = rep(0, n)
  rsave = rep(0, ndraw-nomit)
  vmean = rep(0, n)
  yhat = rep(0, n)
  
  sv = 1.0         
  yin = y             # save original y-values
  bhat = b0           # starting value for beta
  
  #(generated with n = 101 in probit_gd but not in the actual data used to run)
  impmatrix = matrix(0, ndraw-nomit, 6)
  #ximpute = c(1.0000, -1.2766, 0.0400) # x values of the guy that has to be imputed with y = 0
  #yimpute = 0;
  ximpute = c(1.0000, -0.6131, -1.0305) # x values of the guy that has to be imputed with y = 1 
  yimpute = 1;
  
  pb <- txtProgressBar(min = 0, max = ndraw, style = 3)
  
  for(i in 1:ndraw){
    xstar = x*c(sqrt(V))
    ystar = y*c(sqrt(V))
    
    lp=xstar%*%bhat 
    
    # sample Z from normal truncated right at 0 if yin = 0  
    # sample Z from normal truncated left at 0 if yin = 1
    # mean of truncated is the predicted value of yi (XiTB) stored in lp
    for(j in 1:length(yin)){
      if(yin[j] == 0){
        y[j] = rtruncnorm(1, a=-Inf, b = 0, mean = lp[j], sd = 1)
      } else if(yin[j] == 1){
        y[j] = rtruncnorm(1, a=0, b = Inf, mean = lp[j], sd = 1)
      }
    }
    
    # now, lets draw the truncated value for y_mis, which has true value 0
    # in the first occasion;
    if (yimpute == 0){
      xstarimp = ximpute*sqrt(1)
      lpimp = xstarimp%*%bhat; # calculate the mean for the truncated (XtB)
      yimptrue = rtruncnorm(1, a=-Inf, b = 0, mean = lpimp, sd = 1) # right truncated at 0 is proper 
      yimpfalse = rtruncnorm(1, a=0, b = Inf, mean = lpimp, sd = 1)
    }else{
      xstarimp = ximpute*sqrt(1)
      lpimp = xstarimp%*%bhat; # calculate the mean for the truncated (XtB)
      yimptrue = rtruncnorm(1, a=0, b = Inf, mean = lpimp, sd = 1) # left truncated at 0 is proper 
      yimpfalse = rtruncnorm(1, a=-Inf, b = 0, mean = lpimp, sd = 1)
    }
    
    # update beta
    xpxi = solve(t(xstar)%*%xstar + Q) 
    xpy = t(xstar)%*%ystar + Qpc
    bhat = xpxi%*%xpy
    for(j in 1:k){
      bhat[j] = bhat[j] +  rnorm(1, xpxi[j,j], 1)
    }
    
    # update V
    e = y - x%*%bhat  
    
    # calculate e for imputed individual;       
    eimptrue = yimptrue - xstarimp%*%bhat 
    eimpfalse = yimpfalse - xstarimp%*%bhat
    
    # i am assuming we need 100 values with df rval+1
    chiv = rchisq(n, rval+1)
    vi = ((e*e) + IN*rval) / chiv
    V = IN/vi
    if(mm != 0){
      rval = gamm_rnd(1, mm, kk)  # update rval
    }             
    
    if(i > nomit){ #if we are past burn-in, save the draws
      bsave[i-nomit,] = t(bhat)
      ymean = ymean + lp
      vmean = vmean + vi
      yhat = yhat + stdn_cdf(c(y))
      
      # save imputation values to matrix
      impmatrix[i-nomit, 1] = lpimp;
      impmatrix[i-nomit, 2] = yimptrue;
      impmatrix[i-nomit, 3] = yimpfalse;
      impmatrix[i-nomit, 4] = xstarimp%*%bhat;
      impmatrix[i-nomit, 5] = eimptrue;
      impmatrix[i-nomit, 6] = eimpfalse;
      
      if(mm != 0){
        rsave[i-nomit] = rval
      }
    } # end of if i > nomit
    setTxtProgressBar(pb, i)
  } # End the sampling
  close(pb)
  
  vmean = vmean/(ndraw-nomit)
  ymean = ymean/(ndraw-nomit)
  yhat = yhat/(ndraw-nomit)
  
  bmean = colMeans(bsave)
  
  # compute McFadden R-squared
  tmp = yin[yin == 1] # find ones
  P = length(tmp)
  cnt0 = n-P
  cnt1 = P
  P = P/n             # proportion of 1's
  like0 = n*(P*log(P) + (1-P)*log(1-P))     # restricted likelihood
  like1 = pr_like(bmean,yin,x);            # unrestricted Likelihood
  r2mf = 1-(abs(like1)/abs(like0));         # McFadden pseudo-R2 
  # compute Estrella R-squared
  term0 = (2/n)*like0
  term1 = 1/(abs(like1)/abs(like0))^term0
  rsqr = 1-term1                            # Estrella R2
  
  # return results
  results <- NULL
  results$meth  = "probit_g";
  results$r2mf = r2mf;
  results$rsqr = rsqr;
  results$bdraw = bsave;
  results$pmean = c;
  results$pstd  = sqrt(diag(T));
  results$vmean = vmean;
  results$ymean = ymean;
  results$yhat = yhat;
  results$imputation = impmatrix;
  if(mm != 0){
    results$rdraw = rsave;
    results$m     = mm;
    results$k     = kk;
  }else{
    results$r     = rval;
    results$rdraw = rsave;
  }
  results$nobs  = n;
  results$nvar  = k;
  results$y     = yin;
  results$x     = x;
  results$ndraw = ndraw;
  results$nomit = nomit;
  results$pflag = "plevel";
  return(results)
}
