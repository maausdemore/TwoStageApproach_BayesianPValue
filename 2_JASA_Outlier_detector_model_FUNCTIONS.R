

##### 2_JASA_Outlier_detector_model_FUNCTIONS.R for submission
##### The following functions are to be used with the file 2_JASA_Outlier_detector_model_FUNCTIONS.R.

###################################################################################################################
###                                               Begin Document                                                ###
###################################################################################################################

###################################################################################################################
###                                                LIBRARY CHECK                                                ###
###################################################################################################################
### Determine if appropriate packages are installed:
### If not, install them; If so, load them.
###################################################################################################################
##### Install packages that are not installed on the system...
## For creating pseudo-spectra
if(!"fda" %in% installed.packages()){install.packages("fda")}
if(!"Matrix" %in% installed.packages()){install.packages("Matrix")}
## For estimating model parameters
if(!"mvtnorm" %in% installed.packages()){install.packages("mvtnorm")}
if(!"MCMCpack" %in% installed.packages()){install.packages("MCMCpack")}
## For parallel computing
if(!"parallel" %in% installed.packages()){install.packages("parallel")}
if(!"proxy" %in% installed.packages()){install.packages("proxy")}
## For baseline adjustment of spectra
if(!"baseline" %in% installed.packages()){install.packages("baseline")}
## For plotting
if(!"ggplot2" %in% installed.packages()){install.packages("ggplot2")}
if(!"gridExtra" %in% installed.packages()){install.packages("gridExtra")}
## For smoothing XRF spectra
if(!"signal" %in% installed.packages()){install.packages("signal")}
###################################################################################################################
### Load packages
library(fda)
library(Matrix)
library(mvtnorm)
library(MCMCpack)
library(parallel)
library(proxy)
library(baseline)
library(ggplot2)
library(gridExtra)
###################################################################################################################

###################################################################################################################
##### FUNCTIONS: TABLE OF CONTENTS
###################################################################################################################
# 1) Functions to obtain H-Vals..............................................................................57-373
# 2) Functuons to create FTIR spectra.......................................................................377-464
# 3) Functions to calculate score between spectra...........................................................468-606
# 4) Functions to obtain c(alpha)...........................................................................610-645
# 5) Functions to obtain power of metric....................................................................649-749
# 6) Functions to obtain RMP of metric......................................................................753-812
###################################################################################################################

###################################################################################################################
### 1. Functions to obtain H-Vals
###################################################################################################################

###################################################################################################################
### FUNCTION: P.fun
### INPUT: n;
### OUTPUT: NxN design matrix, P, for outlier detector model
###################################################################################################################
P.fun <- function(n){
   n.exp <- t(combn(n,2))
   P <- apply(n.exp, 1, function(x,n){tmp <- rep(0,n); tmp[x] <- 1; return(tmp)}, n)
   return(t(P))
}
###################################################################################################################

###################################################################################################################
### FUNCTION: estimates.for.outlier.detection.fun
### INPUT: a vector of scores, S; a design matrix P for within-source only; the number of control objects, N;
### OUTPUT: estimates of parameters theta, sigma a, and sigma e, for the outlier detector model; sums of squares
###         estimates for SSa and SSe
###################################################################################################################
#### estimate parameters of within-source score model (and outlier detector)
estimates.for.outlier.detection.fun <- function(S,P,N){
   ## get number of comparisons
   n <- choose(N,2)

   ## get mean estimate
   theta.hat <- mean(S)

   ## get SSt
   #SSt <- as.vector(S.vect)%*%(diag(N)-rep(1,N)%*%t(rep(1,N))/N)%*%as.vector(S.vect)
   SSt <- sum((S - theta.hat)^2)

   ## get SSa
   tmp.PS <- colSums(P*S)
   SSa <- max((N-1)^2/(N-2)*sum((tmp.PS[tmp.PS!=0]/(N-1) - theta.hat)^2),0)
   MSa <- SSa/(N-1)

   ## get SSe
   #print(c(SSt,SSa))
   SSe <- max(SSt-SSa,0)
   MSe <- SSe/(n-N)

   ## get sig.e.hat
   sig.e.hat <- max(MSe,0)

   ## get sig.a.hat
   sig.a.hat <- max((MSa - MSe)/(N-2),0)

   return(list(theta.hat=theta.hat,sig.a.hat=sig.a.hat,sig.e.hat=sig.e.hat,SSa=SSa,SSe=SSe))
}
###################################################################################################################

###################################################################################################################
### FUNCTION: get.partial.sigma.inv.outlier.detection.fun
### INPUT: the number of control objects, N; a design matrix P=NULL;
### OUTPUT: inverse of the covariance matrix for a set of scores; returns partial matrices for section 3.2.3 in
###         Doug draft paper for VSP fibers (given new values for sigma a and sigma e, allows for calculating the
###         inverse of Sigma in an efficient manner)
###################################################################################################################
get.partial.sigma.inv.outlier.detection.fun <- function(N, P=NULL){
   ## get P matrix if necessary
   if (is.null(P)){
      P <- P.fun(N)
   }

   ## get number of comparisons considered
   n <- choose(N,2)

   ## part 1 for sigma
   part1 <- matrix(1,ncol=1,nrow=n)%*%matrix(1,ncol=n,nrow=1)/n

   ## part 2 for sigma
   part2 <- (N-1)^2/(N-2)*(1/(N-1)*P-1/n*matrix(1,nrow=n,ncol=1)%*%matrix(1,nrow=1,ncol=N))%*%(1/(N-1)*t(P)-1/n*matrix(1,nrow=N,ncol=1)%*%matrix(1,nrow=1,ncol=n))

   ## part 3 for sigma
   part3 <- diag(n) - part1 - part2

   return(list(part1=part1, part2=part2, part3=part3))
}
###################################################################################################################


###################################################################################################################
### FUNCTION: sample.var.outlier.detection.fun
### INPUT: the number of samples to be obtained, N.samples; the sums of squares for a, SSa; the sums of squares for
###        e, SSe; the degrees of freedom for SSa, df.a=n0-1; the degrees of freedom for SSe, N-n0; the value of
###        alpha for the prior distribution of sig.a, alpha.SSa=1; the value of alpha for the prior distribution of
###        sig.e, alpha.SSe=1; the value of beta for the prior distribution of sig.a, beta.SSa=1; the value of beta
###        for the prior distribution of sig.e, beta.SSe=1; a parameter to determine if samples from the posterior
###        distributions of sig.a and sig.e will be returned, return.sample=TRUE
### OUTPUT: a sample of sigma.a and sigma.e for outlier detection model; posterior means for eta.a and eta.e if
###         return.sample=FALSE
###################################################################################################################
sample.var.outlier.detection.fun <- function(N.samples, SSa, SSe, df.a = n0-1, df.e=choose(n0,2)-n0, alpha.SSa=1, alpha.SSe=1, beta.SSa=1, beta.SSe=1, return.sample=TRUE){
   ## get alpha value for posterior distribution
   post.alpha.SSa <- alpha.SSa + df.a/2
   post.alpha.SSe <- alpha.SSe + df.e/2

   # get beta value for posterior distribution
   post.beta.SSa <- beta.SSa+SSa/2
   post.beta.SSe <- beta.SSe+SSe/2

   # get mean of the posterior distributions
   post.mean.eta.e <- post.beta.SSe/(post.alpha.SSe-1)
   post.mean.eta.a <- post.beta.SSa/(post.alpha.SSa-1)

   if (return.sample){
      ## get eta samples from posterior distribution
      eta.a.samps <- rinvgamma(n=N.samples, shape=post.alpha.SSa, scale=post.beta.SSa)
      eta.e.samps <- rinvgamma(n=N.samples, shape=post.alpha.SSe, scale=post.beta.SSe)

      ## get the variance from eta
      design.matrix <- matrix(c(df.a-1,1,0,1), nrow=2, ncol=2, byrow=TRUE)
      sigma.samps <- solve(design.matrix)%*%rbind(eta.a.samps,eta.e.samps)

      ## kill samples that have one of the variances < 0
      tmp <- which(sigma.samps[1,] < 0)
      if (length(tmp)>0){
         sigma.samps[1,tmp] <- sigma.samps[2,tmp]/2
      }
      return(sigma.samps)
      } else {
      return(c(post.mean.eta.a,post.mean.eta.e))
   }
}
###################################################################################################################

###################################################################################################################
### FUNCTION: sample.mu.outlier.detection.fun
### INPUT: a matrix of samples of sigma a and sigma e, sigma; a vector of scores, S.vect; the number of control
###        objects, N; a partial matrix for easy computation of Sigma, sig.N.inv.parts; a prior for the variance
###        componenets, lambda=1; a prior for the mean, mu.0=0; a parameter to determing if the samples will be
###        returned in the output, return.sample=TRUE;
### OUTPUT: a sample from the distribution of mu; posterior means for mu if return.sample=FALSE
###################################################################################################################
sample.mu.outlier.detection.fun <- function(sigma, S.vect, N, sig.N.inv.parts, lambda=1, mu.0=0, return.sample=TRUE){
   ## initalise some stuff
   n <- choose(N,2)
   sigma.a <- sigma[1]
   sigma.e <- sigma[2]

   ## get the inverse sigma
   Sig.mat.inv <- 1/(2*(N-1)*sigma.a+sigma.e)*sig.N.inv.parts$part1 + 1/((N-2)*sigma.a+sigma.e)*sig.N.inv.parts$part2 + 1/sigma.e*sig.N.inv.parts$part3

   ## get the one vector
   one.vect <- rep(1, n)

   ## get variance of the posterior distribution for mu
   variance <- (one.vect%*%Sig.mat.inv%*%one.vect+(lambda*2*sigma.a+sigma.e)^-1)^-1

   ## get mean of the posterior distribution for mu
   mean.num <- S.vect%*%Sig.mat.inv%*%one.vect+(lambda*2*sigma.a+sigma.e)^-1*mu.0
   mean.post <- mean.num*variance

   ## get posterior sample for mu
   sampl <- rnorm(1, mean.post, sqrt(variance))

   ## if we want a sample
   if (return.sample){
      return(sampl)
   } else {
      return(mean.post)
   }
}
###################################################################################################################

###################################################################################################################
### FUNCTION: get.lik.from.score.outlier.detection.fun
### INPUT: a vector of scores, S.vect; a sample from the posterior distribution of the mean, mu; a sample from the
###        posterior distribution of sigma a, sig.e; a sample from the posterior distribution of sigma e, sig.e;
###        the number of control objects, N; a partial matrix for easy computation of Sigma, sig.inv.parts;
### OUTPUT: the log density for a vector of scores
###################################################################################################################
get.lik.from.score.outlier.detection.fun <- function(S.vect, mu, sig.a, sig.e, N, sig.inv.parts){
   ## initialise variables
   n <- choose(N,2)
   lambda1 <- 2*(N-1)*sig.a + sig.e
   lambda2 <- (N-2)*sig.a + sig.e
   lambda3 <- sig.e

   ## get the determinant and normalisation constant
   detr <- log(lambda1) + (N-1)*log(lambda2) + (n-N)*log(lambda3) + n*log(2*pi)

   ## get the SS
   SS <- (n*(mean(S.vect)-mu)^2 )/lambda1 + (t(S.vect)%*%sig.inv.parts$part2%*%S.vect)/lambda2 +  (t(S.vect)%*%sig.inv.parts$part3%*%S.vect)/lambda3

   ## get the log likelihood of the density
   ldens <- (detr + SS)*-0.5

   return(ldens)
}
###################################################################################################################

###################################################################################################################
### FUNCTION: sample.Sm.given.Sn.outlier.detection.fun
### INPUT: a vector of within-control-source scores, Sn; the number of samples to obtain, N.Sm.star; a sample from
###        the posterior distribution of the mean, mu; a sample from the posterior distribution of sigma a, sig.a;
###        a sample from the posterior distribution of sigma e, sig.e; the number of control objects, N; the number
###        of trace objects, M; the design matrix of all comparisons being made, PPt.M; a partial matrix for easy
###        computation of Sigma, sig.N.inv.parts;
### OUTPUT: a sample of Sm given Sn for determinging the H-Val
###################################################################################################################
sample.Sm.given.Sn.outlier.detection.fun <- function(Sn, N.Sm.star, mu, sig.a, sig.e, N, M, PPt.M, sig.N.inv.parts){
   ## define some variables
   n <- choose(N,2)
   m <- choose(N+M,2) - n

   ## get the full covariance matrix
   sig.full <- PPt.M*sig.a + diag(choose(N+M,2))*sig.e

   ## get the partial bits
   sig.m.m <- sig.full[1:m,1:m]
   sig.n.m <- sig.full[(m+1):(m+n),1:m]
   sig.m.n <- t(sig.n.m)

   ## get the inverse of sig.n.n
   sig.n.n.inv <- 1/(2*(N-1)*sig.a+sig.e)*sig.N.inv.parts$part1 + 1/((N-2)*sig.a+sig.e)*sig.N.inv.parts$part2 + 1/sig.e*sig.N.inv.parts$part3

   ## get the conditional mean
   if (length(mu)==1){
      mu.cond <- rep(mu,m) + sig.m.n%*%sig.n.n.inv%*%(Sn-rep(mu,n))
   } else {
      mu.cond <- mu[1:m] + sig.m.n%*%sig.n.n.inv%*%(Sn-mu[(m+1):length(mu)])
   }

   ## get the conditional cov
   cov.cond <- sig.m.m - sig.m.n%*%sig.n.n.inv%*%sig.n.m

   ## get the choleski decomposition
   chol.cov.cond <- Matrix::chol(cov.cond)

   ## sample a bunch of z~N(0,1)
   z <- rnorm(m*N.Sm.star,0,1)
   z <- matrix(z,m,N.Sm.star)

   ## get a sample from the Sm|Sn
   Sm.given.Sn <- as.vector(mu.cond)+t(chol.cov.cond)%*%z

   return(Sm.given.Sn)
}
###################################################################################################################

###################################################################################################################
### FUNCTION: indicator.outlier.detection.fun
### INPUT: the vector of trace-comparisons, Sm; the vector of within-control-source comparisons, Sn; a sample of
###        scores obtained by using the estimates of the parameters, Sm.star; a sample from the posterior
###        distribution of the mean, mu; a sample from the posterior distributions of sigma a and sigma e,
###        sig.samples; the number of control objects, N; the number of trace objects, M; the design matrix of
###        comparisons being made, PPt; a partial matrix for easy computation of Sigma, sig.NM.inv.parts; a partial
###        matrix for easy computation of Sigma, sig.N.inv.parts;
### OUTPUT: one or zero, depending on if the conditional log-likelihood of the observed vector of scores is larger
###        than the conditional log-likelihood of the sampled vector of scores
###################################################################################################################
indicator.outlier.detection.fun <- function(Sm, Sn, Sm.star, mu, sig.samples, N, M, PPt, sig.NM.inv.parts, sig.N.inv.parts){
   ## assign samples to sig.a and sig.e
   sig.a <- sig.samples[1]
   sig.e <- sig.samples[2]

   ## get the log likelihood of the conditional distribution Sm|Sn for the observed dtaa
   loglike.Sn <- get.lik.from.score.outlier.detection.fun(Sn, mu=mu, sig.a=sig.a, sig.e=sig.e, N=N, sig.inv.parts=sig.N.inv.parts)
   cond.loglike.real <- get.lik.from.score.outlier.detection.fun(c(Sm,Sn), mu=mu, sig.a=sig.a, sig.e=sig.e, N=N+M, sig.inv.parts=sig.NM.inv.parts) - loglike.Sn

   ## get the log likelihood of the conditinal distribution Sm|Sn for the resampled data
   cond.loglike.star <- get.lik.from.score.outlier.detection.fun(c(Sm.star, Sn), mu=mu, sig.a=sig.a, sig.e=sig.e, N=N+M, sig.inv.parts=sig.NM.inv.parts) - loglike.Sn

   ## get 1 if loglike(obs) >= loglike(star)
   return(ifelse(cond.loglike.real >= cond.loglike.star, 1, 0))
}
###################################################################################################################

###################################################################################################################
### FUNCTION: H.fun
### INPUT: the vector of trace-comparisons, Sm; the vector of within-control-source comparisons, Sn; the number of
###        trace objects, M;the number of control objects, N; the number of samples to observe from the posterior
###        distributions; P.M=NULL; P.N=NULL; sig.N.inv.parts=NULL; sig.NM.inv.parts=NULL;
### OUTPUT: Bayesian p-value: the probability of observing a sample more extreme than the one we've observed
###################################################################################################################
H.fun <- function(Sm, Sn, M, N, N.samples, P.M=NULL, P.N=NULL, sig.N.inv.parts=NULL, sig.NM.inv.parts=NULL){
   ## get number of control-source comparisons
   n <- choose(N,2)

   ## create the matrices if they dont exist
   if (is.null(P.M)){
      P.M <- P.fun(N+M)
   }
   PPt.M <- P.M%*%t(P.M)
   if (is.null(P.N)){
      P.N <- P.fun(N)
   }
   if (is.null(sig.N.inv.parts)){
      sig.N.inv.parts <- get.partial.sigma.inv.outlier.detection.fun(N,P.N)
   }
   if (is.null(sig.NM.inv.parts)){
      sig.NM.inv.parts <- get.partial.sigma.inv.outlier.detection.fun(N+M,P.M)
   }

   ## get the parameters from the control objects from within-scores ONLY
   param <- estimates.for.outlier.detection.fun(Sn,P.N,N)

   ## sample values for Psi = (theta, sig.a, sig.e)
   ## sample sig.a, sig.e and mu from the posterior of P(mu, sig.a, sig.e | SS)
   cte <- 0.001
   sig.samples <- sample.var.outlier.detection.fun(N.samples, SSa=param$SSa, SSe=param$SSe, df.a = N-1, df.e=n-N, alpha.SSa=2+cte, alpha.SSe=2+cte, beta.SSa=(1+cte)*param$SSa/(N-1), beta.SSe=(1+cte)*param$SSe/(n-N))
   ## sample mu
   mu.tilda <- apply(sig.samples,2,sample.mu.outlier.detection.fun,  S.vect=Sn, N=N, sig.N.inv.parts=sig.N.inv.parts, lambda=1, mu.0=param$theta.hat)

   ## sample Sm.star
   Sm.star <- sapply(1:N.samples, function(x,Sn,mu.tilda,sig.samples,N,M,PPt.M,sig.N.inv.parts){sample.Sm.given.Sn.outlier.detection.fun(Sn, N.Sm.star=1, mu=mu.tilda[x], sig.a=sig.samples[1,x], sig.e=sig.samples[2,x], N, M, PPt.M, sig.N.inv.parts)},Sn,mu.tilda,sig.samples,N,M,PPt.M,sig.N.inv.parts)

   ## calculate the h.val
   h.val <- mean(sapply(1:N.samples,function(x,Sm,Sn,Sm.star,mu.tilda,sig.samples,N,M,PPt.M,sig.NM.inv.parts, sig.N.inv.parts){indicator.outlier.detection.fun(Sm, Sn, Sm.star=Sm.star[,x], mu=mu.tilda[x], sig.samples=sig.samples[,x], N, M, PPt.M, sig.NM.inv.parts, sig.N.inv.parts)},Sm,Sn,Sm.star,mu.tilda,sig.samples,N,M,PPt.M,sig.NM.inv.parts, sig.N.inv.parts))

   return(h.val)
}
###################################################################################################################



###################################################################################################################
### 2. Functions to create FTIR spectra
###################################################################################################################

###################################################################################################################
### FUNCTION: create.FTIR.objects.from.real.spectra.outlier.detection.fun
### INPUT: the number of basis functions used to create the pseudo-spectra, p.basis; the number of points that the
###        spectra will be evaluated at, p.eval=500; the number of objects from the contorl source, N; the number
###        of objects from the trace source, M; the matrix of training spectra from a single source where the first
###        column contains values for the x-axis, spectra; a parameter indicating if we want the M trace samples to
###        be from the same source as the N control samples, same=TRUE; p.pca=7; pca.diff=0; the vector of points
###        at which the pseudo-spectra will be evaluated, xs.eval=NULL;
### OUTPUT: a matrix of pseudo-spectra obtained by learning the coefficients of basis functions from the real data
###################################################################################################################
create.FTIR.objects.from.real.spectra.outlier.detection.fun <- function(p.basis, p.eval=500, N, M, spectra, same=TRUE, p.pca=7, pca.diff=0, xs.eval=NULL){
   ## get the range of the xs axis of the spectra
   range.spectra <- range(spectra[,1])

   ## sort the specrtra matrix along the x axis by ascending order
   ix.spectra <- order(spectra[,1])
   spectra <- spectra[ix.spectra,]


   ## evaluate the functions at selected points
   if (!is.null(xs.eval)){
      ix.spectra <- ix.spectra[1:length(xs.eval)]-min(ix.spectra[1:length(xs.eval)])+1
      xs.eval <- xs.eval[ix.spectra]
   } else {
      ix.spectra <- ix.spectra[1:p.eval]-min(ix.spectra[1:p.eval])+1
      xs.eval <- seq(spectra[1,1],spectra[dim(spectra)[1],1],length=p.eval)
   }

   ##### model using fda
   ## create the bases
   bases <- create.bspline.basis(rangeval=range(spectra[,1]),nbasis=p.basis, norder=4)
   ## evaluate the bases functions at some points along the range of values of interest
   basis.eval.points <- seq(range(spectra[,1])[1], range(spectra[,1])[2],length=dim(spectra)[1])
   basismat <- eval.basis(basis.eval.points, bases)

   ## get the coefficients for the basis functions
   basis.coef <- solve(crossprod(basismat), crossprod(basismat, data.matrix(spectra[,2:dim(spectra)[2]])))
   spectra.fd <- fd(basis.coef,bases)

   ## use the coefficients to create a distribution and resample
   mean.coef <- rowMeans(basis.coef)
   cov.coef <- (cov(t(basis.coef)))
   if (sum(diag(cov.coef))<0.001){
      cov.coef <- cov.coef*0.001/sum(diag(cov.coef))
   }
   if (same | (pca.diff==0)){
      ## resample coef from that distribution
      basiscoef.sample <- t(rmvnorm(N+M,mean.coef,cov.coef))
      ## get the functional objects
      xfd = fd(basiscoef.sample, bases)
      tmp.spectra.mat <- cbind(xs.eval,as.matrix((eval.fd(xs.eval,xfd))))
   } else {
      ## resample N coef from the distribution for the control samples
      N.basiscoef.sample <- t(rmvnorm(N,mean.coef,cov.coef))

      ## perform PCA
      pca.spectra <- pca.fd(spectra.fd, p.pca)

      ## change the scores before the reconstruction
      weights <- rep(0,p.pca)
      weights[2] <- pca.diff
      reconstructed.coef <- rowMeans((pca.spectra$harmonics$coefs)%*%(t(pca.spectra$scores)+weights)) + pca.spectra$meanfd$coefs

      ## resample M coef from the distribution for the trace samples
      M.basiscoef.sample <- t(rmvnorm(M,reconstructed.coef,cov.coef))

      ## put all the coefficients together
      basiscoef.sample <- cbind(M.basiscoef.sample,N.basiscoef.sample)

      ## get the functional objects
      xfd = fd(basiscoef.sample, bases)
      tmp.spectra.mat <- cbind(xs.eval,as.matrix((eval.fd(xs.eval,xfd))))
   }

   tmp.spectra.mat <- tmp.spectra.mat[ix.spectra,]

   ## prepare a vector of names for outlier detector
   e.names <- rep(c('e.u','e.s'),times=c(M,N))
   e.names <- paste(e.names,c((1:M), (1:N)),sep=".")
   colnames(tmp.spectra.mat) <- c("xs.eval",e.names)

   return(tmp.spectra.mat)
}
###################################################################################################################



###################################################################################################################
### 3. Functions to calculate scores between objects
###################################################################################################################

###################################################################################################################
### FUNCTION: FTIR.dist.SCALAR.fun
### INPUT: a pair of spectral objects to be compared, spectra; the number of control objects, N; the number of
###        trace objects, M; the kernel function to be used to compare the pair of spectra, kern.fun; a list of
###        additional parameters to be used in kern.fun, dots; the functions to be sourced, my.functions;
### OUTPUT: a list of (1) within-control-source scores, Sn, and (2) trace-source scores, Sm
###################################################################################################################
FTIR.dist.SCALAR.fun <- function(spectra, N, M, kern.fun, dots, my.functions){
   ## get combinations
   my.grid <- t(combn(M+N,2))
   ## get indices for within-control and trace object scores
   ix.known.source.scores <- which(my.grid[,1]>M)
   between.source.scores <- which(my.grid[,1]<=M & my.grid[,2]>M)

   ## get gram matrix for all data
   K.mat <- gram.mat.fun(spectra, kern.fun, dots)
   S.vect <- K.mat[lower.tri(K.mat)]

   ## Assign within-control scores and trace object scores
   Sn <- S.vect[ix.known.source.scores]
   Sm <- S.vect[-ix.known.source.scores]

   return(list(Sn=Sn, Sm=Sm))
}
###################################################################################################################

###################################################################################################################
### FUNCTION: gram.mat.fun
### INPUT: a set of spectra, dat; the kernel function to be used to compare the pair of spectra, kern.fun; a list
###        of additional parameters to be used in kern.fun, dots=NULL;
### OUTPUT: the "Gram matrix" of all pairwise comparisons
###################################################################################################################
gram.mat.fun <- function(dat, kern.fun, dots=NULL){
   # create distance metric
   if (!pr_DB$entry_exists("kern.fun")){
      pr_DB$set_entry(FUN = kern.fun, names = c("kern.fun"))
   }else {
      pr_DB$modify_entry(FUN = kern.fun, names = c("kern.fun"))
   }

   # get data in matrix format
   dat <- as.matrix(dat)

   # create Gram matrix:
   # get off diagona values of Gram matrix
   K.mat <- as.matrix(proxy::dist(dat,method='kern.fun',dots=dots,diag=TRUE,by_rows=FALSE,upper = TRUE))
   # get diagonal values of Gram matrix
   #tmp.Kmat <- apply(dat,2,FUN=function(x,kern.fun,dots){kern.fun(x,x,dots=dots)},kern.fun=kern.fun,dots=dots)
   # create diagonal of Gram matrix
   #diag(K.mat) <- tmp.Kmat

   # return Gram matrix
   return((K.mat))
}
#################################################################################################################

###################################################################################################################
### FUNCTION: kern.CEDRIC.fun
### INPUT: the first spectra, spectra1; the second spectra to be compared to the first spectra, spectra2; the list
###        of additional parameters to be used to calculate the score between spectra1 and spectra2, dots;
### OUTPUT: the "Gram matrix" of all pairwise comparisons
###################################################################################################################
kern.CEDRIC.fun <- function(spectra1, spectra2, dots){
   ## get indices to consider
   ix.pair <- index.pair.fun(spectra1, spectra2, wav.num=dots$xs.eval)

   ## cross correlation kernel
   tmp.ccf <- ccf(spectra1[ix.pair],spectra2[ix.pair],lag.max=dots$max.lag,plot=FALSE)$acf
   tmp.ccf.4 <- (2*dots$max.lag+1 - sqrt(tmp.ccf%*%tmp.ccf))*10

   ## prepare the normalised euclidean kernel
   spectra1.norm <- (spectra1[ix.pair]-min(spectra1[ix.pair]))
   spectra1.norm <- spectra1.norm/(sum(spectra1.norm*mean(abs(diff(dots$wav.num)))))*length(spectra1.norm)
   spectra2.norm <- (spectra2[ix.pair]-min(spectra2[ix.pair]))
   spectra2.norm <- spectra2.norm/(sum(spectra2.norm*mean(abs(diff(dots$wav.num)))))*length(spectra2.norm)
   ## get euclidean norm
   tmp.norm.euc <- sqrt(t(spectra1.norm-spectra2.norm)%*%(spectra1.norm-spectra2.norm))/length(spectra1.norm)*10000

   return(log(tmp.ccf.4*tmp.norm.euc)*1000000)
}
###################################################################################################################

###################################################################################################################
### FUNCTION: index.indiv.fun
### INPUT: a baseline corrected spectra, corrected.ghoul; the x-axis of points at which corrected.ghoul is
###        to be evaluated, wav.num;
### OUTPUT: a set of indices of the points to be kept for the input baseline corrected spectra
###################################################################################################################
index.indiv.fun <- function(corrected.ghoul, wav.num){
   ## get first derivate
   first.deriv <- diff(corrected.ghoul)/mean(diff(wav.num))*100
   n <- round(length(corrected.ghoul)/50)
   quantile.param <- 0.6

   ## get the moving average for the intensity (we keep it large to 'protect' the peaks)
   ma.inty <- as.vector(filter(corrected.ghoul,rep(1/n,n)))
   ma.inty[is.na(ma.inty)] <- 0

   ## get the moving average for the derivative (we keep it large to 'protect' the peaks)
   ma.deriv <- as.vector(filter(c(0,first.deriv),rep(1/n,n)))
   ma.deriv[is.na(ma.deriv)] <- 0

   ## select the epsilon for the intensity threshold
   epsilon.inty <- quantile(ma.inty,quantile.param,na.rm=TRUE)
   ## select the epsilon for the derivative threshold
   epsilon.deriv <- quantile(abs(ma.deriv),quantile.param,na.rm=TRUE)
   ## get the index of the points we want to keep for deriv
   ix.deriv <- which(abs(ma.deriv)>epsilon.deriv)

   ## get the index of the points we want to keep for inty
   ix.inty <- which(ma.inty>epsilon.inty)
   ix.keep.ghoul <- sort(union(ix.deriv,ix.inty))

   return(ix.keep.ghoul)
}
###################################################################################################################

###################################################################################################################
### FUNCTION: index.pair.fun
### INPUT: a baseline corrected spectra, corrected.ghoul.1; a second baseline corrected spectra to be compared to
###        the corrected.ghoul.1, corrected.ghoul.2; the x-axis of points at which corrected.ghoul is to be
###        evaluated, wav.num;
### OUTPUT: a set of indices of the points to be kept for a pair of baseline corrected spectra
###################################################################################################################
index.pair.fun <- function(corrected.ghoul.1, corrected.ghoul.2, wav.num){
   ## get indices to keep for first spectra
   ix.keep.ghoul.1 <- index.indiv.fun(corrected.ghoul.1, wav.num)
   ## get indices to keep for second spectra
   ix.keep.ghoul.2 <- index.indiv.fun(corrected.ghoul.2, wav.num)
   ## get union of two sets of indices
   ix.keep.ghoul <- sort(union(ix.keep.ghoul.1,ix.keep.ghoul.2))

   return(ix.keep.ghoul)
}
###################################################################################################################



###################################################################################################################
### 4. Functions to obtain c(alpha)
###################################################################################################################

###################################################################################################################
### FUNCTION: HvalCluster.unconditional.c.alpha.FTIR.fun
### INPUT: a set of spectral objects, spectra; the number of trace objects, M; the number of control objects, N; the
###        number of samples to be obtained from the posterior distributions of the parameters mu, sigma a, and
###        sigma e, N.samples; the design matrix of all pairwise comparisons, P.M; the design matrix of all within-
###        control-source comparisons, P.N; a partial matrix for easy computation of Sigma, sig.N.inv.parts; a
###        partial matrix for easy computation of Sigma, sig.NM.inv.parts; the kernel function to be used to compare
###        pairs of spectra, kern.fun; the FTIR distance function to be used to compute the scores between a set of
###        spectra, FTIR.dist;=FTIR.dist.SCALAR.fun; the functions to be sourced, my.functions; ...
### OUTPUT: the H-Val associated with a set of N+M spectra
###################################################################################################################
HvalCluster.unconditional.c.alpha.FTIR.fun <- function(spectra, M, N, N.samples, P.M, P.N, sig.N.inv.parts, sig.NM.inv.parts, kern.fun, FTIR.dist = FTIR.dist.SCALAR.fun, my.functions, ...){
   ## get dots info
   dots <- list(...)
   if(length(dots)!=0){
      tmp <- as.list(substitute(list(...)))[-1L]
      for(i in 1:length(dots)){
         if(nchar(names(dots)[[i]])==0){
            names(dots)[[i]] <- as.character(tmp[[i]])
         }
      }
   }

   ## get scores
   my.scores <- FTIR.dist(spectra, N, M, kern.fun, dots, my.functions)

   ## get h.val
   h.val <- H.fun(my.scores$Sm, my.scores$Sn, M, N, N.samples, P.M=P.M, P.N=P.N, sig.N.inv.parts=sig.N.inv.parts, sig.NM.inv.parts=sig.NM.inv.parts)

   return(h.val)
}
###################################################################################################################



###################################################################################################################
### 5. Functions to obtain power of metric
###################################################################################################################

###################################################################################################################
### FUNCTION: POWER.fun
### INPUT: a pair of spectral sources to be differentiated, one.combn; the number of trace objects, M; the number
###        of control objects, N; the number of samples to be obtained from the posterior distributions of mu, sig
###        a and sig.e, N.samples; the number of simulations to be run, N.sims; the design matrix of all pairwise
###        comparisons, P.M; the design matrix of all within-control-source comparisons, P.N; a partial matrix for
###        easy computation of Sigma, sig.N.inv.parts; a partial matrix for easy computation of Sigma,
###        sig.NM.inv.parts; the number of basis funtions used to create a set of pseudo-spectra, p.basis; the
###        number of points at which the pseudo-spectra will be evaluated, p.eval; the list of true spectra,
###        spectra.abs.ls; the kernel function to be used to determine the value of the score between a pair of
###        spectra, kern.fun; the FTIR distance function to be used to compute the scores between a set of spectra,
###        FTIR.dist; the functions to be sourced, my functions; the list of additional parameters to be used in
###        kern.fun, dots;
### OUTPUT: the power associated with a set of samples from a pair of sources
###################################################################################################################
POWER.fun <- function(one.combn, M, N, N.samples, N.sims, P.M, P.N, sig.N.inv.parts, sig.NM.inv.parts, p.basis, p.eval, spectra.abs.ls, kern.fun, FTIR.dist, my.functions, dots){
   ## define which source is control and which is trace
   one.combn <- as.numeric(one.combn)
   trace.source <- one.combn[1]
   control.source <- one.combn[2]

   ## create a list of length N.sims of samples from the fixed trace source
   trace.spectra.ls <- lapply(rep(trace.source, N.sims), function(x, p.basis, p.eval, N, M, spectra.abs.ls, wav.num, xs.eval) {create.FTIR.objects.from.real.spectra.outlier.detection.fun(p.basis=p.basis, p.eval=p.eval, N=N, M=M, spectra=cbind(wav.num, spectra.abs.ls[[x]]), xs.eval=xs.eval)[,-1]}, p.basis=p.basis, p.eval=p.eval, N=1, M=M-1, spectra.abs.ls=spectra.abs.ls, wav.num=dots$wav.num, xs.eval=dots$xs.eval)

   ## create a list of length N.sims of samples from the current control source
   control.spectra.ls <- lapply(rep(control.source, N.sims), function(x, p.basis, p.eval, N, M, spectra.abs.ls, wav.num, xs.eval) {create.FTIR.objects.from.real.spectra.outlier.detection.fun(p.basis, p.eval, N, M, spectra=cbind(wav.num, spectra.abs.ls[[x]]), xs.eval=xs.eval)[,-1]}, p.basis=p.basis, p.eval=p.eval, N=N-1, M=1, spectra.abs.ls=spectra.abs.ls, wav.num=dots$wav.num, xs.eval=dots$xs.eval)

   ## combine the columns of the two lists
   combined.spectra.ls <- lapply(1:N.sims, function(x, trace.spectra.ls, control.spectra.ls){cbind(trace.spectra.ls[[x]], control.spectra.ls[[x]])}, trace.spectra.ls, control.spectra.ls)
   combined.spectra.ls <- lapply(combined.spectra.ls, function(x, M, N){colnames(x) <- c(paste("e.u.", 1:M, sep=""), paste("e.s.", 1:N, sep="")); return(x)}, M, N)

   ## get the scores
   my.scores <- lapply(combined.spectra.ls, function(x, N, M, kern.fun, dots, my.functions){FTIR.dist(x, N, M, kern.fun, dots, my.functions)}, N, M, kern.fun, dots, my.functions)

   ## separate out between scores (we want to return these later)
   control.mean.spectra <- rowMeans(spectra.abs.ls[[control.source]])[dots$ix]
   trace.mean.spectra <- rowMeans(spectra.abs.ls[[trace.source]])[dots$ix]
   scores.between <- kern.fun(control.mean.spectra,trace.mean.spectra,dots=dots)

   ## get H.vals
   h.val <- lapply(my.scores, function(x, M, N, N.samples, P.M=P.M, P.N=P.N, sig.N.inv.parts=sig.N.inv.parts, sig.NM.inv.parts=sig.NM.inv.parts){H.fun(x$Sm, x$Sn, M, N, N.samples, P.M, P.N, sig.N.inv.parts, sig.NM.inv.parts)}, M, N, N.samples, P.M=P.M, P.N=P.N, sig.N.inv.parts=sig.N.inv.parts, sig.NM.inv.parts=sig.NM.inv.parts)

   ## create a matrix of H.vals and scores.between to return
   results.mat <- list(comparison=one.combn, h.val=unlist(h.val), scores.between=scores.between)

   return(results.mat)
}

###################################################################################################################
### FUNCTION: POWER.overlaid.fun
### INPUT: a pair of spectral sources to be differentiated, one.combn; the number of trace objects, M; the number
###        of control objects, N; the number of samples to be obtained from the posterior distributions of mu, sig
###        a and sig.e, N.samples; the design matrix of all pairwise comparisons, P.M; the design matrix of all
###        within-control-source comparisons, P.N; a partial matrix for easy computation of Sigma, sig.N.inv.parts;
###        a partial matrix for easy computation of Sigma, sig.NM.inv.parts; the number of basis funtions used to
###        create a set of pseudo-spectra, p.basis; the number of points at which the pseudo-spectra will be evaluated,
###        p.eval; the list of true spectra, spectra.abs.ls; the kernel function to be used to determine the value of
###        the score between a pair of spectra, kern.fun; the FTIR distance function to be used to compute the scores
###        between a set of spectra, FTIR.dist; the functions to be sourced, my functions; the list of additional
###        parameters to be used in kern.fun, dots;
### OUTPUT: the power associated with a set of samples from a pair of sources
###################################################################################################################
POWER.overlaid.fun <- function(one.combn, M, N, N.samples, P.M, P.N, sig.N.inv.parts, sig.NM.inv.parts, p.basis, p.eval, spectra.abs.ls, kern.fun, FTIR.dist, my.functions, dots){
   ## define which source is control and which is trace
   one.combn <- as.numeric(one.combn)
   trace.source <- one.combn[1]
   control.source <- one.combn[2]

   ## create a list of length N.sims of samples from the fixed trace source
   trace.spectra.ls <- create.FTIR.objects.from.real.spectra.outlier.detection.fun(p.basis=p.basis, p.eval=p.eval, N=1, M=M-1, spectra=cbind(dots$wav.num, spectra.abs.ls[[trace.source]]), xs.eval=dots$xs.eval)[,-1]

   ## create a list of length N.sims of samples from the current control source
   control.spectra.ls <- create.FTIR.objects.from.real.spectra.outlier.detection.fun(p.basis, p.eval, N=N-1, M=1, spectra=cbind(dots$wav.num, spectra.abs.ls[[control.source]]), xs.eval=dots$xs.eval)[,-1]

   ## combine the columns of the two lists
   combined.spectra.ls <- cbind(trace.spectra.ls, control.spectra.ls)
   colnames(combined.spectra.ls) <- c(paste("e.u.", 1:M, sep=""), paste("e.s.", 1:N, sep=""))

   ## get the scores
   my.scores <- FTIR.dist(combined.spectra.ls, N, M, kern.fun, dots, my.functions)

   ## separate out between scores (we want to return these later)
   control.mean.spectra <- rowMeans(spectra.abs.ls[[control.source]])[dots$ix]
   trace.mean.spectra <- rowMeans(spectra.abs.ls[[trace.source]])[dots$ix]

   ## get average between score
   scores.between <- kern.fun(control.mean.spectra,trace.mean.spectra,dots=dots)

   ## get H.val (low H.val -> reject H0)
   h.val <- H.fun(my.scores$Sm, my.scores$Sn, M, N, N.samples, P.M, P.N, sig.N.inv.parts, sig.NM.inv.parts)

   ## create a matrix of H.vals and scores.between to return
   results.vect <- c(h.val=unlist(h.val), scores.between=scores.between)

   return(results.vect)
}
###################################################################################################################



###################################################################################################################
### 6. Functions to obtain RMP of metric
###################################################################################################################

###################################################################################################################
### FUNCTION: RMP.fun
### INPUT: the pair of spectral sources to be considered, spectra.ix; the list of true spectra, spectra.abs.ls; the
###        number of basis functions to be used to create the pseudo-spectra, p.basis; the number of points at which
###        the pseudo-spectra will be evaluated, p.eval; the number of control objects, N; the number of trace
###        objects, M; the design matrix of all within-control-source comparisnons, P.N; the design matrix of all
###        pairwise comparisons, P.M; a partial matrix for easy computation of Sigma, sig.N.inv.parts; a partial
###        matrix for easy computation of Sigma, sig.NM.inv.parts; the number of samples to be obtained from the
###        posterior distributions of mu, sigma a and sigma e, N.samples; the functions to be sourced, my.functions;
###        the kernel function to be used to determine the score between two spectra, kern.fun; the FTIR distance
###        function to be used to compute the scores between a set of spectra, FTIR.dist; the number of Simulations
###        to be run, n.sims=1; the list of additional parameters to be used in kern.fun, dots;
### OUTPUT: the RMP associated with a pair of spectral sources
###################################################################################################################
RMP.fun <- function(spectra.ix, spectra.abs.ls, trace.objects = NULL, p.basis, p.eval, N, M, P.N, P.M, sig.N.inv.parts, sig.NM.inv.parts, N.samples, my.functions,  kern.fun, FTIR.dist, n.sims=1, dots){
   ## which source are we considering?
   which.source <- spectra.ix
   comparison.sources <- (1:length(spectra.abs.ls))[-which.source]

   ## prepare the return
   h.val <- matrix(0,n.sims,length(spectra.abs.ls)-1)
   colnames(h.val) <- comparison.sources

   ## begin sims
   for (i in 1:n.sims){

      ## determine if trace objects are fixed are random and assign accordingly
      if (is.null(trace.objects)){
         trace.objects.tmp <- create.FTIR.objects.from.real.spectra.outlier.detection.fun(p.basis=p.basis, p.eval=p.eval, N=1, M=M-1, spectra=cbind(dots$wav.num, spectra.abs.ls[[which.source]]), xs.eval=dots$xs.eval)[,-1]
      } else {
         trace.objects.tmp <- trace.objects[[which.source]]
      }

      ## create a cluster to generate trace objects
      n.nodes <- min(detectCores()-1, 7)
      cl <- makeCluster(n.nodes)
      ## source functions in each worker
      clusterCall(cl, fun=function(my.functions){source(my.functions)},my.functions)
      ## get a set of comparison objects
      comparison.spectra.ls <- parLapplyLB(cl, comparison.sources, fun = function(x, p.basis, p.eval, N, M, spectra.abs.ls, wav.num, xs.eval){create.FTIR.objects.from.real.spectra.outlier.detection.fun(p.basis=p.basis, p.eval=p.eval, N=N, M=M, spectra=cbind(wav.num, spectra.abs.ls[[x]]),xs.eval, same=TRUE)[,-1]}, p.basis = p.basis, p.eval=p.eval, N=N-1, M=1, spectra.abs.ls=spectra.abs.ls, wav.num = dots$wav.num, xs.eval=dots$xs.eval)

      ## combine two sets of spectra
      RMP.spectra.ls <- lapply(1:length(comparison.sources), function(x, trace.objects.tmp, comparison.spectra.ls){cbind(trace.objects.tmp, comparison.spectra.ls[[x]])}, trace.objects.tmp, comparison.spectra.ls)

      ## get the scores
      my.scores <- lapply(RMP.spectra.ls, function(x, N, M, kern.fun, dots, my.functions){FTIR.dist(x, N=N, M=M, kern.fun=kern.fun, dots=dots, my.functions=my.functions)}, N, M, kern.fun, dots, my.functions)

      ## get H.val (low H.val -> reject H0)
      HVals <- parLapplyLB(cl, my.scores, function(x, M, N, N.samples, P.M, P.N, sig.N.inv.parts, sig.NM.inv.parts){H.fun(x$Sm, x$Sn, M, N, N.samples, P.M, P.N, sig.N.inv.parts, sig.NM.inv.parts)}, M, N, N.samples, P.M, P.N, sig.N.inv.parts, sig.NM.inv.parts)
      stopCluster(cl)
      h.val[i,] <- unlist(HVals)
   }

   return(h.val)
}
###################################################################################################################

###################################################################################################################
###                                                End Document                                                 ###
###################################################################################################################
