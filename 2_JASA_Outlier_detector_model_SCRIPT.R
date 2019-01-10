

##### 2_JASA_Outlier_detector_model_SCRIPT.R for submission
##### The following code reproduces the data portrayed in the Section 6 - WORKED EXAMPLE in the submitted paper,
##### "Two-Stage Approach for the Inference of the Source of High-Dimensional and Complex Chemical Data in Forensic
##### Science". This script is to be used with the file 2_JASA_Outlier_detector_model_FUNCTIONS.R.

###################################################################################################################
###                                               Begin Document                                                ###
###################################################################################################################

###################################################################################################################
##### Table of Contents
###################################################################################################################
# 1) Observe Distribution of Scores........................................................................36 - 146
# 2) Observe Distribution of Test Statistic...............................................................150 - 299
# 3) Power Analysis.......................................................................................303 - 494
# 4) RMP Analysis.........................................................................................498 - 716
# 5) Additional Plots in Paper............................................................................720 - 854
###################################################################################################################

##### NOTE: It is necessary that the user define three file paths throughout the document. This is the only aspect
##### of the document that requires any changes on the part of the user. There are three file paths per section
##### that require your personal file paths. They are:
#####
##### (1) my.path:      This is the file path associated with the location on your computer where your data is
#####                   saved;
##### (2) my.functions: This is the file path associated with the lcoation on your computer where your functions
#####                   are saved;
##### (3) my.results:   This is the file path associated with the location on your computer where you would like
#####                   to save your results.
#####
##### These paths should be defined on lines 46 - 48 in section 1, lines 160 - 162 in section 2, lines 309 - 311 in
##### section 3, lines 504 - 506 in section 4, and lines 724 - 726 in section 5

###################################################################################################################
### 1. Observe Distribution of Scores
###################################################################################################################

############################## 1. Set up workspace

##### Start from fresh workspace
rm(list=ls())

##### Define paths
my.path <- "~/Dropbox/2_PhD/1_Project/PhD_project/2_data/0_Cyril_Sources_168_Samples_Transmittance" #PUT FILE PATH WHERE DATA IS STORED HERE
my.functions <- "~/Dropbox/2_PhD/1_Project/PhD_project/0_code/0_Outlier_2_stages/2_JASA_Outlier_detector_model_FUNCTIONS.R" #PUT FILE PATH WHERE FUNCTIONS ARE SAVED HERE
my.results <- "~/Desktop/Test" #PUT FILE PATH WHERE YOU WOULD LIKE TO SAVE YOUR DATA HERE

##### Source functions
source(my.functions)

############################# DO NOT CHANGE ANYHING *IN THIS SECTION* BELOW THIS LINE #############################

############################## 2. Prepare data

##### Read in files
my.files <- list.files(path=my.path, pattern=".csv", ignore.case=TRUE)
spectra.ls <- list()
for(i in 1:length(my.files)){
   spectra.ls[[i]] <- read.csv(file.path(my.path,my.files[i]), header=FALSE, sep=",", dec=".")
   # remove zeros
   spectra.ls[[i]][spectra.ls[[i]]==0] <- min(apply(spectra.ls[[i]], 2, function(x) min(x[x!=0]))[-1])
}

############### NOTE
##### All data is in transmittance mode:
##### To get absorbance from transmittance: A = -log10(T);
##### To get transmittance from absorbance: T = 10^(-A);
##########################################################

##### Convert the spectra from transmittance to absorbance
spectra.abs.ls <- list()
for(i in 1:length(spectra.ls)){
   spectra.abs.ls[[i]] <- cbind(spectra.ls[[i]][,1],-1*log10(spectra.ls[[i]][,2:dim(spectra.ls[[i]])[2]]/100))
}

##### Separate x-axis (wavenumber) data from y-axis (absorbance) data
### Create an object to store x-axis (wavenumber)
dots <- list()
dots$wav.num <- spectra.abs.ls[[1]][,1]
### Get rid of x-axis for all elements in list so that we just have spectra in columns
spectra.abs.ls <- lapply(spectra.abs.ls, function(x){x[,2:ncol(x)]})

##### Perform baseline correction on spectra
spectra.abs.ls <- lapply(spectra.abs.ls, function(x){t(baseline.modpolyfit(t(x))$corrected)})

##### Convert lists into matrices
spectra.abs.real.mat <- matrix(unlist(spectra.abs.ls),dim(spectra.abs.ls[[1]])[1],length(spectra.abs.ls)*dim(spectra.abs.ls[[1]])[2])

##### Define some variables
### How many control objects do we want to consider? (N+M=7)
N <- 4
### How many trace objects do we want to consider? (N+M=7)
M <- 3
### How many basis functions do we want to use to define the spectra?
p.basis <- 300
### How many points do we want to use to evaluate the pseudo-spectra?
p.eval <- 1000

### At which specific points do we want to evaluate the spectra?
dots$ix <- round(seq(1,length(dots$wav.num),length=p.eval))
dots$xs.eval <- dots$wav.num[dots$ix]

##### Down-sample the real data
spectra.abs.real.mat <- spectra.abs.real.mat[dots$ix,]

############################## 3. Obtain "Gram vectors" of scores

##### How many sources are we considering?
N.sources <- 166
which.source <- 1:N.sources

##### Define items for "dots" that haven't already been defined
dots$max.lag <- 10

##### Calculate the Gram matrix for real data using a cluster
### Define a grid of source combinations
comb.grid.real <- combn(1:(N.sources*(M+N)), 2, simplify=TRUE)
### Defne the number of nodes we will use
n.nodes <- detectCores()-1
### Create the workers in each cluster
cl <- makeCluster(n.nodes)
### Souce the functions in each cluster
clusterCall(cl, fun=function(my.functions){source(my.functions)}, my.functions)
### Obtain "Gram vector" for data
Gram.real.vect <- parApply(cl, comb.grid.real, 2, function(grid, dat, kern.fun, dots){kern.fun(dat[,grid[1]], dat[,grid[2]], dots)}, dat=spectra.abs.real.mat, kern.fun=kern.CEDRIC.fun, dots=dots)
### Terminate the cluster
stopCluster(cl)

############################## 4. Observe distributions of scores - Corresponds to Figure A.1 in paper

##### Look at triplets of scores - real data
### Locate triplets of scores in pre-existing data - we will have two disjoint triplets per source
pred.mat <- diag(N.sources)%x%rbind(cbind(diag(2)%x% matrix(1,3,3),0),0)
ix.triplets <- pred.mat[lower.tri(pred.mat)]
triplet.scores.real <- Gram.real.vect[ix.triplets==1]
dim(triplet.scores.real) <- c(3, N.sources*2)

### Observe marginal plots of triplets- Real Data
get.comb <- split(t(combn(3,2)), rep(1:3))
pl <- lapply(get.comb, function(x,y){qplot(y[,x[1]], y[,x[2]], xlab=paste("Dimension", x[1]), ylab=paste("Dimension", x[2]),main="Original Space") + theme(plot.title = element_text(size=17), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15))},y=(t(triplet.scores.real)))
spectra.pca <- princomp(t(triplet.scores.real))
pl.pca <- lapply(get.comb, function(x,y){qplot(y[,x[1]], y[,x[2]], xlab=paste("Eigendimension", x[1]), ylab=paste("Eigendimension", x[2]), main="Eigenspace") + theme(plot.title = element_text(size=17), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15))}, y=spectra.pca$scores)
marrangeGrob(c(pl,pl.pca), nrow=2, ncol=3, top=NULL, layout_matrix=matrix(c(1,2,3,4,5,6),2,3,byrow=TRUE))
###################################################################################################################



###################################################################################################################
### 2. Observe Distribution of Test Statistic
###################################################################################################################

############################## 1. Set up workspace

##### Start from fresh workspace
rm(list=ls())

##### Define paths
my.path <- "~/Dropbox/2_PhD/1_Project/PhD_project/2_data/0_Cyril_Sources_168_Samples_Transmittance" #PUT FILE PATH WHERE DATA IS STORED HERE
my.functions <- "~/Dropbox/2_PhD/1_Project/PhD_project/0_code/0_Outlier_2_stages/2_JASA_Outlier_detector_model_FUNCTIONS.R" #PUT FILE PATH WHERE FUNCTIONS ARE SAVED HERE
my.results <- "~/Desktop/Test" #PUT FILE PATH WHERE YOU WOULD LIKE TO SAVE YOUR DATA HERE

##### Source functions
source(my.functions)

############################# DO NOT CHANGE ANYHING *IN THIS SECTION* BELOW THIS LINE #############################

############################## 2. Prepare real data

##### Read in files
my.files <- list.files(path=my.path, pattern=".csv", ignore.case=TRUE)
spectra.ls <- list()
for(i in 1:length(my.files)){
   spectra.ls[[i]] <- read.csv(file.path(my.path,my.files[i]), header=FALSE, sep=",", dec=".")
   # remove zeros
   spectra.ls[[i]][spectra.ls[[i]]==0] <- min(apply(spectra.ls[[i]], 2, function(x) min(x[x!=0]))[-1])
}

############### NOTE
#### All data is in transmittance mode:
##### To get absorbance from transmittance: A = -log10(T);
##### To get transmittance from absorbance: T = 10^(-A);

##### Convert the spectra from transmittance to absorbance
spectra.abs.ls <- list()
for(i in 1:length(spectra.ls)){
   spectra.abs.ls[[i]] <- cbind(spectra.ls[[i]][,1],-1*log10(spectra.ls[[i]][,2:dim(spectra.ls[[i]])[2]]/100))
}

##### Create an object to store x-axis (wavenumber)
dots <- list()
dots$wav.num <- spectra.abs.ls[[1]][,1]

##### Get rid of x-axis for all elements in list so that we just have spectra in columns
spectra.abs.ls <- lapply(spectra.abs.ls, function(x){x[,2:ncol(x)]})

##### Perform baseline correction on spectra
spectra.abs.ls <- lapply(spectra.abs.ls, function(x){t(baseline.modpolyfit(t(x))$corrected)})

############################## 4. Run a simulaiton to obtain test statistics

##### Assign values to variables for simulating spectra
### How many control objects do we want to consider?
n <- c(5,10,15)
### How many trace objects do we want to consider?
M <- 3
### How many basis functions do we want to use to define the spectra?
p.basis <- 300
### How many points do we want to use to evaluate the pseudo-spectra?
p.eval <- 1000
### At which specific points do we want to evaluate the spectra?
dots$ix <- round(seq(1,length(dots$wav.num),length=p.eval))
dots$xs.eval <- dots$wav.num[dots$ix]

##### Define items for "dots" that haven't already been defined
dots$max.lag <- 10

##### Sample sources for the simulations
N.sources <- 2000
which.source <- sample(1:length(spectra.abs.ls), N.sources, replace=TRUE)

##### Begin simulation
for(N in n){
   ##### Create a list of sampled spectra using a cluster
   ### Define the number of nodes we will use
   n.nodes <- detectCores()-1
   ### Create the workers in the cluster
   cl <- makeCluster(n.nodes)
   ### Source functions in each worker
   clusterCall(cl, fun=function(my.functions){source(my.functions)},my.functions)
   ### Create the spectra
   spectra.sims.ls <- parLapplyLB(cl, which.source, fun = function(x, p.basis, p.eval, N, M, spectra.abs.ls, wav.num, xs.eval){create.FTIR.objects.from.real.spectra.outlier.detection.fun(p.basis=p.basis, p.eval=p.eval, N=N, M=M, spectra=cbind(wav.num, spectra.abs.ls[[x]]), xs.eval=xs.eval)[,-1]},p.basis=p.basis, p.eval=p.eval, N=N, M=M, spectra.abs.ls=spectra.abs.ls, wav.num=dots$wav.num, xs.eval=dots$xs.eval)
   ### Terminate the cluster
   stopCluster(cl)

   ##### Define the P matrices
   P.M <- P.fun(N+M)
   P.N <- P.fun(N)

   ##### Pre-define the parts for the two inverse covariance matrices
   sig.N.inv.parts <- get.partial.sigma.inv.outlier.detection.fun(N)
   sig.NM.inv.parts <- get.partial.sigma.inv.outlier.detection.fun(N+M)

   ##### Begin H.val simulations
   ### How many samples do we want to consider?
   N.samples <- 10000

   ### Define the number of nodes we will use
   n.nodes <- detectCores()-1
   ### Create the workers in the cluster
   cl <- makeCluster(n.nodes)
   ### Source functions in each worker
   clusterCall(cl, fun=function(my.functions){source(my.functions)},my.functions)
   ### Generate p-values from considering different sources
   HVal.Cluster.Results <- parLapplyLB(cl, spectra.sims.ls, fun=HvalCluster.unconditional.c.alpha.FTIR.fun, M=M, N=N, N.samples=N.samples, P.M=P.M, P.N=P.N, sig.N.inv.parts=sig.N.inv.parts, sig.NM.inv.parts=sig.NM.inv.parts, kern.fun=kern.CEDRIC.fun, FTIR.dist = FTIR.dist.SCALAR.fun, my.functions=my.functions, wav.num=dots$wav.num, ix=dots$ix, xs.eval=dots$xs.eval, max.lag=dots$max.lag)
   ### Terminate the cluster
   stopCluster(cl)

   ### Save results
   save(HVal.Cluster.Results, file=paste(my.results, "/c_alpha_N", N, "_M3.rdata", sep=""))
}

############################## 5. Observe distributions of test statistics - Corresponds to Figure 3

##### Make Plot for N=5
load(paste(my.results, "/c_alpha_N", 5, "_M3.rdata", sep=""))
tmp5 <- unlist(HVal.Cluster.Results)
c.alpha5 <- quantile(tmp5,0.05)

df <- data.frame(x=tmp5)
p5 <- ggplot(df, aes(x))+stat_ecdf()+geom_abline(slope=1, intercept=0, linetype="dashed")+xlim(0,1)+ylim(0,1)+labs(x="Sample Quantile", y="ECDF", title="N=5, M=3") + theme(plot.title = element_text(size=17), axis.title.x=element_text(size=17), axis.title.y=element_text(size=17))

##### Make plot for N=10
load(paste(my.results, "/c_alpha_N", 10, "_M3.rdata", sep=""))
tmp10 <- unlist(HVal.Cluster.Results)
c.alpha10 <- quantile(tmp10,0.05)

df <- data.frame(x=tmp10)
p10 <- ggplot(df, aes(x))+stat_ecdf()+geom_abline(slope=1, intercept=0, linetype="dashed")+xlim(0,1)+ylim(0,1)+labs(x="Sample Quantile", y="ECDF", title="N=10, M=3") + theme(plot.title = element_text(size=17), axis.title.x=element_text(size=17), axis.title.y=element_text(size=17))

### Make plot for N=15
load(paste(my.results, "/c_alpha_N", 15, "_M3.rdata", sep=""))
tmp15 <- unlist(HVal.Cluster.Results)
c.alpha15 <- quantile(tmp15,0.05)

df <- data.frame(x=tmp15)
p15 <- ggplot(df, aes(x))+stat_ecdf()+geom_abline(slope=1, intercept=0, linetype="dashed")+xlim(0,1)+ylim(0,1)+labs(x="Sample Quantile", y="ECDF", title="N=15, M=3") + theme(plot.title = element_text(size=17), axis.title.x=element_text(size=17), axis.title.y=element_text(size=17))

##### Plot the three plots in a single frame
marrangeGrob(list(p5, p10, p15), nrow=1, ncol=3, top=FALSE)

############################## 6. Observe quantiles of test statistics - Corresponds to Table 1 in paper

##### Get values for table
quantile(tmp5, c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95))
quantile(tmp10, c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95))
quantile(tmp15, c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95))
###################################################################################################################



###################################################################################################################
### 3. Power Analysis
###################################################################################################################

############################## 1. Set up workspace

##### Start from fresh workspace
rm(list=ls())

##### Define paths
my.path <- "~/Dropbox/2_PhD/1_Project/PhD_project/2_data/0_Cyril_Sources_168_Samples_Transmittance" #PUT FILE PATH WHERE DATA IS STORED HERE
my.functions <- "~/Dropbox/2_PhD/1_Project/PhD_project/0_code/0_Outlier_2_stages/2_JASA_Outlier_detector_model_FUNCTIONS.R" #PUT FILE PATH WHERE FUNCTIONS ARE SAVED HERE
my.results <- "~/Desktop/Test" #PUT FILE PATH WHERE YOU WOULD LIKE TO SAVE YOUR DATA HERE

##### Source functions
source(my.functions)

############################# DO NOT CHANGE ANYHING *IN THIS SECTION* BELOW THIS LINE #############################

############################## 2. Prepare data

##### Read in files
my.files <- list.files(path=my.path, pattern=".csv", ignore.case=TRUE)
spectra.ls <- list()
for(i in 1:length(my.files)){
   spectra.ls[[i]] <- read.csv(file.path(my.path,my.files[i]), header=FALSE, sep=",", dec=".")
   # remove zeros
   spectra.ls[[i]][spectra.ls[[i]]==0] <- min(apply(spectra.ls[[i]], 2, function(x) min(x[x!=0]))[-1])
}

############### NOTE
#### All data is in transmittance mode:
##### To get absorbance from transmittance: A = -log10(T);
##### To get transmittance from absorbance: T = 10^(-A);

##### Convert the spectra from transmittance to absorbance
spectra.abs.ls <- list()
for(i in 1:length(spectra.ls)){
   spectra.abs.ls[[i]] <- cbind(spectra.ls[[i]][,1],-1*log10(spectra.ls[[i]][,2:dim(spectra.ls[[i]])[2]]/100))
}

##### Create an object to store x-axis (wavenumber)
dots <- list()
dots$wav.num <- spectra.abs.ls[[1]][,1]

##### Get rid of x-axis for all elements in list so that we just have spectra in columns
spectra.abs.ls <- lapply(spectra.abs.ls, function(x){x[,2:ncol(x)]})

##### Perform baseline correction on spectra
spectra.abs.ls <- lapply(spectra.abs.ls, function(x){t(baseline.modpolyfit(t(x))$corrected)})

############################## 3. Power study

##### Assign values to variables for simulating spectra
### How many control objects do we want to consider?
n <- c(5,10,15)
### How many trace objects do we want to consider?
M <- 3
### How many basis functions do we want to use to define the spectra?
p.basis <- 300
### How many points do we want to use to evaluate the pseudo-spectra?
p.eval <- 1000
### At which specific points do we want to evaluate the spectra?
dots$ix <- round(seq(1,length(dots$wav.num),length=p.eval))
dots$xs.eval <- dots$wav.num[dots$ix]

##### Define items for "dots" that hasn't already been defined
dots$max.lag <- 10

##### Sample sources for the simulations
### How many simulations do we want to run?
N.sims <- 1000000
### How many samples will we obtain from the distributions of our paramters?
N.samples <- 10000

for(N in n){
   ##### Define the P matrices
   P.M <- P.fun(N+M)
   P.N <- P.fun(N)
   PPt.M <- P.M%*%t(P.M)

   ##### Pre-define the parts for the two inverse covariance matrices
   sig.N.inv.parts <- get.partial.sigma.inv.outlier.detection.fun(N)
   sig.NM.inv.parts <- get.partial.sigma.inv.outlier.detection.fun(N+M)

   ##### Determine combinations of sources to study the power
   ix.grid <- matrix(sample(1:length(spectra.abs.ls), N.sims*2,  replace=TRUE), ncol=2, nrow=N.sims)

   ##### Obtain power of source defined above
   ### Define the number of nodes we will use
   n.nodes <- detectCores()-1
   ### Create the workers in the cluster
   cl <- makeCluster(n.nodes)
   ### Ensure seeds do not overlap
   clusterSetRNGStream(cl, 7675)
   ### Source functions in each worker
   clusterCall(cl, fun=function(my.functions){source(my.functions)},my.functions)
   ### Obtain power of determined trace.source
   Power.Cluster.Results <- parApply(cl=cl, X=ix.grid, MARGIN=1, FUN=POWER.overlaid.fun, M=M, N=N, N.samples=N.samples, P.M=P.M, P.N=P.N, sig.N.inv.parts=sig.N.inv.parts, sig.NM.inv.parts=sig.NM.inv.parts, p.basis=p.basis, p.eval=p.eval, spectra.abs.ls=spectra.abs.ls, kern.fun=kern.CEDRIC.fun, FTIR.dist=FTIR.dist.SCALAR.fun, my.functions=my.functions, dots=dots)
   ### Terminate the cluster
   stopCluster(cl)

   save(Power.Cluster.Results, file=paste(my.results, "/power_N", N, "_M3.rData", sep=""))
}

############################## 4. Observe power curve - Corresponds to Figure 4 in paper

##### Define values of c.alpha in Table 1 of paper (using values in paper - may not exactly correspond to c.alpha values
##### observed above due to sampling)
c.alpha5 <- 0.002
c.alpha10 <- 0.006
c.alpha15 <- 0.007

##### Make plot for N=5
### Load data
load(paste(my.results, "/power_N", 5, "_M3.rData", sep=""))
# re-structure data
h.vals5 <- t(Power.Cluster.Results)
# replace -Inf by 0
h.vals5[!is.finite(h.vals5[,2]),2] <- min(h.vals5[is.finite(h.vals5[,2]),2])-exp(11)
# order h.vals
ix5 <- order(h.vals5[,2])

### Bin observations
my.bins5 <- seq(range(h.vals5[ix5,2]/100000)[1]-1, range(h.vals5[ix5,2]/100000)[2]+1,4)
bin.means5 <- (my.bins5+2)[-length(my.bins5)]
# get the observations that fall in each bin
obs.in.bins5 <- lapply(1:(length(my.bins5)-1), function(x, my.bins, my.h.vals){which(my.h.vals >= my.bins[x] & my.h.vals < my.bins[x+1])}, my.bins5, h.vals5[,2]/100000)
# get the power for those observations in those bins
power.vals5 <- lapply(obs.in.bins5, function(x, c.alpha, my.h.vals){sum(my.h.vals[unlist(x)]<c.alpha)/length(unlist(x))}, c.alpha=c.alpha5, h.vals5[,1])
power.vals5 <- unlist(power.vals5)
# get rid of NaNs
bin.means5 <- bin.means5[which(!is.nan(power.vals5))]
power.vals5 <- power.vals5[which(!is.nan(power.vals5))]

### Create plot
df.p5 <- data.frame(X1=bin.means5, X2=power.vals5)
pc.OL <- ggplot(df.p5, aes(X1,X2))+geom_point()+geom_line(linetype="dotted")+labs(x="Distance from Mean Spectra Value", y="Power", title="Power Curve Considering All Sources (N=5,10,15, M=3)") + ylim(0,1) + xlim(80,126) + theme(plot.title = element_text(size=17), axis.title.x=element_text(size=17), axis.title.y=element_text(size=17))

##### Make plot for N=10
### Load data
load(paste(my.results, "/power_N", 10, "_M3.rData", sep=""))
# re-structure data
h.vals10 <- t(Power.Cluster.Results)
# replace -Inf by 0
h.vals10[!is.finite(h.vals10[,2]),2] <- min(h.vals10[is.finite(h.vals10[,2]),2])-exp(11)
# order h.vals
ix10 <- order(h.vals10[,2])

### Bin observations
my.bins10 <- seq(range(h.vals10[ix10,2]/100000)[1]-1, range(h.vals10[ix10,2]/100000)[2]+1,4)
bin.means10 <- (my.bins10+2)[-length(my.bins10)]
# get the observations that fall in each bin
obs.in.bins10 <- lapply(1:(length(my.bins10)-1), function(x, my.bins, my.h.vals){which(my.h.vals >= my.bins[x] & my.h.vals < my.bins[x+1])}, my.bins10, h.vals10[,2]/100000)
# get the power for those observations in those bins
power.vals10 <- lapply(obs.in.bins10, function(x, c.alpha, my.h.vals){sum(my.h.vals[unlist(x)]<c.alpha)/length(unlist(x))}, c.alpha=c.alpha10, h.vals10[,1])
power.vals10 <- unlist(power.vals10)
# get rid of NaNs
bin.means10 <- bin.means10[which(!is.nan(power.vals10))]
power.vals10 <- power.vals10[which(!is.nan(power.vals10))]

### Create plot
pc.OL <- pc.OL + geom_point(aes(bin.means10,power.vals10))+geom_line(aes(bin.means10,power.vals10), linetype="dashed")+labs(x="Distance from Mean Spectra Value", y="Power")

##### Make plot for N=15
### Load data
load(paste(my.results, "/power_N", 15, "_M3.rData", sep=""))
# re-structure data
h.vals15 <- t(Power.Cluster.Results)
# replace -Inf by 0
h.vals15[!is.finite(h.vals15[,2]),2] <- min(h.vals15[is.finite(h.vals15[,2]),2])-exp(11)
# order h.vals
ix15 <- order(h.vals15[,2])

### Bin observations
my.bins15 <- seq(range(h.vals15[ix15,2]/100000)[1]-1, range(h.vals15[ix15,2]/100000)[2]+1,4)
bin.means15 <- (my.bins15+2)[-length(my.bins15)]
# get the observations that fall in each bin
obs.in.bins15 <- lapply(1:(length(my.bins15)-1), function(x, my.bins, my.h.vals){which(my.h.vals >= my.bins[x] & my.h.vals < my.bins[x+1])}, my.bins15, h.vals15[,2]/100000)
# get the power for those observations in those bins
power.vals15 <- lapply(obs.in.bins15, function(x, c.alpha, my.h.vals){sum(my.h.vals[unlist(x)]<c.alpha)/length(unlist(x))}, c.alpha=c.alpha15, h.vals15[,1])
power.vals15 <- unlist(power.vals15)
# get rid of NaNs
bin.means15 <- bin.means15[which(!is.nan(power.vals15))]
power.vals15 <- power.vals15[which(!is.nan(power.vals15))]

### Create plot
pc.OL <- pc.OL + geom_point(aes(bin.means15,power.vals15))+geom_line(aes(bin.means15,power.vals15))+labs(x="Distance from Mean Spectra Value", y="Power")

##### Plot results in a single frame
pc.OL
###################################################################################################################



###################################################################################################################
### 4. RMP Analysis
###################################################################################################################

############################## 1. Set up workspace

##### Start from fresh workspace
rm(list=ls())

##### Define paths
my.path <- "~/Dropbox/2_PhD/1_Project/PhD_project/2_data/0_Cyril_Sources_168_Samples_Transmittance" #PUT FILE PATH WHERE DATA IS STORED HERE
my.functions <- "~/Dropbox/2_PhD/1_Project/PhD_project/0_code/0_Outlier_2_stages/2_JASA_Outlier_detector_model_FUNCTIONS.R" #PUT FILE PATH WHERE FUNCTIONS ARE SAVED HERE
my.results <- "~/Desktop/Test" #PUT FILE PATH WHERE YOU WOULD LIKE TO SAVE YOUR DATA HERE

##### Source functions
source(my.functions)

############################# DO NOT CHANGE ANYHING *IN THIS SECTION* BELOW THIS LINE #############################

############################## 2. Prepare data

##### Read in files
my.files <- list.files(path=my.path, pattern=".csv", ignore.case=TRUE)
spectra.ls <- list()
for(i in 1:length(my.files)){
   spectra.ls[[i]] <- read.csv(file.path(my.path,my.files[i]), header=FALSE, sep=",", dec=".")
   # remove zeros
   spectra.ls[[i]][spectra.ls[[i]]==0] <- min(apply(spectra.ls[[i]], 2, function(x) min(x[x!=0]))[-1])
}

############### NOTE
#### All data is in transmittance mode:
##### To get absorbance from transmittance: A = -log10(T);
##### To get transmittance from absorbance: T = 10^(-A);

##### Convert the spectra from transmittance to absorbance
spectra.abs.ls <- list()
for(i in 1:length(spectra.ls)){
   spectra.abs.ls[[i]] <- cbind(spectra.ls[[i]][,1],-1*log10(spectra.ls[[i]][,2:dim(spectra.ls[[i]])[2]]/100))
}

##### Create an object to store x-axis (wavenumber)
dots <- list()
dots$wav.num <- spectra.abs.ls[[1]][,1]

##### Get rid of x-axis for all elements in list so that we just have spectra in columns
spectra.abs.ls <- lapply(spectra.abs.ls, function(x){x[,2:ncol(x)]})

##### Perform baseline correction on spectra
spectra.abs.ls <- lapply(spectra.abs.ls, function(x){t(baseline.modpolyfit(t(x))$corrected)})

############################## 3. RMP study

##### Assign values to variables for simulating spectra
### How many control objects do we want to consider?
n <- c(5,10,15)
### How many trace objects do we want to consider?
M <- 3
# select where we want to evaluate the pseudo-spectra from the functional data
### How many basis functions do we want to use to define the spectra?
p.basis <- 300
### How many points do we want to use to evaluate the pseudo-spectra?
p.eval <- 1000
### At which specific points do we want to evaluate the spectra?
dots$ix <- round(seq(1,length(dots$wav.num),length=p.eval))
dots$xs.eval <- dots$wav.num[dots$ix]

##### Define items for "dots" that hasn't already been defined
dots$max.lag <- 10

##### Define a simulation to study the RMP
### How many samples will we obtain from the distributions of our paramters?
N.samples <- 5000 #10000

### Define the number of repetitions of the RMP for each set of trace objects
n.sims <- 20

### Define a vector of spectra
spectra.ix <- 1:length(spectra.abs.ls)

### Create the trace objects for the fixed trace RMP analysis
trace.objects.ls <- lapply(spectra.ix, FUN = function(x, p.basis, p.eval, N, M, spectra.abs.ls, wav.num, xs.eval){create.FTIR.objects.from.real.spectra.outlier.detection.fun(p.basis=p.basis, p.eval=p.eval, N=N, M=M, spectra=cbind(wav.num, spectra.abs.ls[[x]]), xs.eval, same=TRUE)[,-1]}, p.basis = p.basis, p.eval=p.eval, N=1, M=M-1, spectra.abs.ls=spectra.abs.ls, wav.num = dots$wav.num, xs.eval=dots$xs.eval)


##### Simulation for fixed traces
for (N in n){
   print(N)
   # get P matrices
   P.M <- P.fun(N+M)
   P.N <- P.fun(N)
   PPt.M <- P.M%*%t(P.M)

   # get the parts for the two inverse covariance matrices
   sig.N.inv.parts <- get.partial.sigma.inv.outlier.detection.fun(N)
   sig.NM.inv.parts <- get.partial.sigma.inv.outlier.detection.fun(N+M)

   # define a cluster
   n.nodes <- 5
   cl <- makeCluster(n.nodes)
   clusterCall(cl, fun=function(my.functions){source(my.functions)},my.functions)
   # obtain the RMPs
   RMP.Cluster.Results <- parLapplyLB(cl, spectra.ix, RMP.fun, spectra.abs.ls=spectra.abs.ls, trace.objects=trace.objects.ls, p.basis=p.basis, p.eval=p.eval, N=N, M=M, P.N=P.N, P.M=P.M, sig.N.inv.parts=sig.N.inv.parts, sig.NM.inv.parts=sig.NM.inv.parts, N.samples=N.samples, my.functions=my.functions,  kern.fun=kern.CEDRIC.fun, FTIR.dist=FTIR.dist.SCALAR.fun, n.sims=n.sims, dots=dots)
   stopCluster(cl)

   save(RMP.Cluster.Results, file=paste(my.results, "/RMP_fixed_trace_N", N, "_M3.rData",sep=""))
}


##### Simulation for random traces
for (N in n){
   print(N)

   # get P matrices
   P.M <- P.fun(N+M)
   P.N <- P.fun(N)
   PPt.M <- P.M%*%t(P.M)

   # get the parts for the two inverse covariance matrices
   sig.N.inv.parts <- get.partial.sigma.inv.outlier.detection.fun(N)
   sig.NM.inv.parts <- get.partial.sigma.inv.outlier.detection.fun(N+M)

   # define a cluster
   n.nodes <- 5
   cl <- makeCluster(n.nodes)
   clusterCall(cl, fun=function(my.functions){source(my.functions)},my.functions)
   # obtain the RMPs
   RMP.Cluster.Results <- parLapplyLB(cl, spectra.ix, RMP.fun, spectra.abs.ls=spectra.abs.ls, trace.objects=NULL, p.basis=p.basis, p.eval=p.eval, N=N, M=M, P.N=P.N, P.M=P.M, sig.N.inv.parts=sig.N.inv.parts, sig.NM.inv.parts=sig.NM.inv.parts, N.samples=N.samples, my.functions=my.functions,  kern.fun=kern.CEDRIC.fun, FTIR.dist=FTIR.dist.SCALAR.fun, n.sims, dots=dots)
   stopCluster(cl)

   save(RMP.Cluster.Results, file=paste(my.results, "/RMP_random_trace_N", N, "_M3.rData", sep=""))
}

############################## 4. Observe plots of RMPs - Corresponds to Figures 5 and 6 in paper

##### Define c.alpha for N=5, 10, 15
c.alpha5 <- 0.002
c.alpha10 <- 0.006
c.alpha15 <- 0.007

##### RMP Plots for Fixed Trace
### RMP Plot : N5 - FIXED TRACE
# Load data
load(paste(my.results, "/RMP_fixed_trace_N", 5, "_M3.rData", sep=""))
# Restructure RMP data
RMP <- sapply(RMP.Cluster.Results,FUN=function(x,c.alpha){return(rowSums(x> c.alpha)/165)},c.alpha=c.alpha5)
RMP <- as.vector(RMP)
RMP <- cbind((rep(1:166, each=20)), RMP)
colnames(RMP) <- c("Spectra", "HVal")
RMP <- data.frame(RMP)
# Create RMP plot
RMP5.F <- ggplot(RMP, aes(y=HVal, x=as.factor(Spectra)))+ylim(c(0,0.12)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x="Spectra", y="RMP", title="RMP: Fixed Trace (N=5, M=3)") + theme(plot.title = element_text(size=17), axis.title.x=element_text(size=17), axis.title.y=element_text(size=17))

### RMP Plot : N10 - FIXED TRACE
# Load data
load(paste(my.results, "/RMP_fixed_trace_N", 10, "_M3.rData", sep=""))
# Restructure RMP data
RMP <- sapply(RMP.Cluster.Results,FUN=function(x,c.alpha){return(rowSums(x> c.alpha)/165)},c.alpha=c.alpha10)
RMP <- as.vector(RMP)
RMP <- cbind((rep(1:166, each=20)), RMP)
colnames(RMP) <- c("Spectra", "HVal")
RMP <- data.frame(RMP)
# Create RMP plot
RMP10.F <- ggplot(RMP, aes(y=HVal, x=as.factor(Spectra))) + ylim(c(0,0.12)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x="Spectra", y="RMP", title="RMP: Fixed Trace (N=5, M=3)") + theme(plot.title = element_text(size=17), axis.title.x=element_text(size=17), axis.title.y=element_text(size=17))

### RMP Plot : N15 - FIXED TRACE
# Load data
load(paste(my.results, "/RMP_fixed_trace_N", 15, "_M3.rData", sep=""))
# Restructure RMP data
RMP <- sapply(RMP.Cluster.Results,FUN=function(x,c.alpha){return(rowSums(x> c.alpha)/165)},c.alpha=c.alpha15)
RMP <- as.vector(RMP)
RMP <- cbind((rep(1:166, each=20)), RMP)
colnames(RMP) <- c("Spectra", "HVal")
RMP <- data.frame(RMP)
# Create RMP plot
RMP15.F <- ggplot(RMP, aes(y=HVal, x=as.factor(Spectra))) + ylim(c(0,0.12)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x="Spectra", y="RMP", title="RMP: Fixed Trace (N=5, M=3)") + theme(plot.title = element_text(size=17), axis.title.x=element_text(size=17), axis.title.y=element_text(size=17))

### Plot fixed RMP results in a single frame
marrangeGrob(list(RMP5.F, RMP10.F, RMP15.F), nrow=3, ncol=1, top=FALSE)

##### RMP Plots for Random Trace
### RMP Plot : N5 - RANDOM TRACE
# Load data
load(paste(my.results, "/RMP_random_trace_N", 5, "_M3.rData", sep=""))
# Restructure RMP data
RMP <- sapply(RMP.Cluster.Results,FUN=function(x,c.alpha){return(rowSums(x> c.alpha)/165)},c.alpha=c.alpha5)
RMP <- as.vector(RMP)
RMP <- cbind(rep(1:166, each=20), RMP)
colnames(RMP) <- c("Spectra", "HVal")
RMP <- data.frame(RMP)
# Create RMP plot
RMP5.R <- ggplot(RMP, aes(y=HVal, x=as.factor(Spectra)))+ylim(c(0,0.12))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(x="Spectra", y="RMP", title="RMP: Fixed Trace (N=5, M=3)")+theme(plot.title = element_text(size=17), axis.title.x=element_text(size=17), axis.title.y=element_text(size=17))

### RMP Plot : N10 - RANDOM TRACE
# Load data
load(paste(my.results, "/RMP_random_trace_N", 10, "_M3.rData", sep=""))
# Restructure RMP data
RMP <- sapply(RMP.Cluster.Results,FUN=function(x,c.alpha){return(rowSums(x> c.alpha)/165)},c.alpha=c.alpha10)
RMP <- as.vector(RMP)
RMP <- cbind((rep(1:166, each=20)), RMP)
colnames(RMP) <- c("Spectra", "HVal")
RMP <- data.frame(RMP)
# Create RMP plot
RMP10.R <- ggplot(RMP, aes(y=HVal, x=as.factor(Spectra)))+ylim(c(0,0.12))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(x="Spectra", y="RMP", title="RMP: Fixed Trace (N=5, M=3)")+theme(plot.title = element_text(size=17), axis.title.x=element_text(size=17), axis.title.y=element_text(size=17))

### RMP Plot : N15 - RANDOM TRACE
# Load data
load(paste(my.results, "/RMP_random_trace_N", 15, "_M3.rData", sep=""))
# Restructure data
RMP <- sapply(RMP.Cluster.Results,FUN=function(x,c.alpha){return(rowSums(x> c.alpha)/165)},c.alpha=c.alpha15)
RMP <- as.vector(RMP)
RMP <- cbind((rep(1:166, each=20)), RMP)
colnames(RMP) <- c("Spectra", "HVal")
RMP <- data.frame(RMP)
# Create RMP plot
RMP15.R <- ggplot(RMP, aes(y=HVal, x=as.factor(Spectra)))+ylim(c(0,0.12))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(x="Spectra", y="RMP", title="RMP: Fixed Trace (N=5, M=3)")+theme(plot.title = element_text(size=17), axis.title.x=element_text(size=17), axis.title.y=element_text(size=17))

### Plot random RMP results in a single frame
marrangeGrob(list(RMP5.R, RMP10.R, RMP15.R), nrow=3, ncol=1, top=FALSE)
###################################################################################################################



###################################################################################################################
### 5. Additional Plots in Paper
###################################################################################################################

##### Start from fresh workspace
rm(list=ls())

##### Define paths
my.path <- "~/Dropbox/2_PhD/1_Project/PhD_project/2_data/0_Cyril_Sources_168_Samples_Transmittance" #PUT FILE PATH WHERE DATA IS STORED HERE
my.functions <- "~/Dropbox/2_PhD/1_Project/PhD_project/0_code/0_Outlier_2_stages/2_JASA_Outlier_detector_model_FUNCTIONS.R" #PUT FILE PATH WHERE FUNCTIONS ARE SAVED HERE
my.results <- "~/Desktop/Test" #PUT FILE PATH WHERE YOU WOULD LIKE TO SAVE YOUR DATA HERE

##### Source functions
source(my.functions)

############################# DO NOT CHANGE ANYHING *IN THIS SECTION* BELOW THIS LINE #############################

############################## 2. Prepare real data

##### Read in files
my.files <- list.files(path=my.path, pattern=".csv", ignore.case=TRUE)
spectra.ls <- list()
for(i in 1:length(my.files)){
   spectra.ls[[i]] <- read.csv(file.path(my.path,my.files[i]), header=FALSE, sep=",", dec=".")
   # remove zeros
   spectra.ls[[i]][spectra.ls[[i]]==0] <- min(apply(spectra.ls[[i]], 2, function(x) min(x[x!=0]))[-1])
}

############### NOTE
#### All data is in transmittance mode:
##### To get absorbance from transmittance: A = -log10(T);
##### To get transmittance from absorbance: T = 10^(-A);


##### Convert the spectra from transmittance to absorbance
spectra.abs.ls <- list()
for(i in 1:length(spectra.ls)){
   spectra.abs.ls[[i]] <- cbind(spectra.ls[[i]][,1],-1*log10(spectra.ls[[i]][,2:dim(spectra.ls[[i]])[2]]/100))
}

##### Create an object to store x-axis (wavenumber)
dots <- list()
dots$wav.num <- spectra.abs.ls[[1]][,1]

##### Get rid of x-axis for all elements in list so that we just have spectra in columns
spectra.abs.ls <- lapply(spectra.abs.ls, function(x){x[,2:ncol(x)]})

##### Perform baseline correction on spectra
spectra.abs.ls <- lapply(spectra.abs.ls, function(x){t(baseline.modpolyfit(t(x))$corrected)})

############################## 3. Prepare pseudo data

##### Assign values to variables for simulating spectra
### How many control objects do we want to consider?
N <- 4
### How many trace objects do we want to consider?
M <- 3
### How many basis functions do we want to use to define the spectra?
p.basis <- 300
### How many points do we want to use to evaluate the pseudo-spectra?
p.eval <- 1000
### At which specific points do we want to evaluate the spectra?
dots$ix <- round(seq(1,length(dots$wav.num),length=p.eval))
dots$xs.eval <- dots$wav.num[dots$ix]

##### Define items for "dots" that hasn't already been defined
dots$max.lag <- 10

##### Sample sources for the simulations
N.sources <- 166
which.source <- 1:N.sources

##### Create a list of sampled spectra using a cluster
### Define the number of nodes we will use
n.nodes <- detectCores()-1
### Create the workers in the cluster
cl <- makeCluster(n.nodes)
### Source functions in each worker
clusterCall(cl, fun=function(my.functions){source(my.functions)},my.functions)
### Create the spectra
spectra.sims.ls <- parLapplyLB(cl, which.source, fun = function(x, p.basis, p.eval, N, M, spectra.abs.ls, wav.num, xs.eval){create.FTIR.objects.from.real.spectra.outlier.detection.fun(p.basis=p.basis, p.eval=p.eval, N=N, M=M, spectra=cbind(wav.num, spectra.abs.ls[[x]]), xs.eval=xs.eval)[,-1]},p.basis=p.basis, p.eval=p.eval, N=N, M=M, spectra.abs.ls=spectra.abs.ls, wav.num=dots$wav.num, xs.eval=dots$xs.eval)
### Terminate the cluster
stopCluster(cl)

############################## 4. Create plots

#####  Figure 1 - Sampled Spectra versus True Spectra
### Define which source we would like to plot
which.true.source <- 5
which.H1.source <- which.true.source

### Create plot of real spectra overlaid with pseudo spectra
# Real spectra
pl.spec <- ggplot() + geom_line(aes(dots$xs.eval, spectra.abs.ls[[which.true.source]][dots$ix,1]),size=0.25, colour="navy") + labs(x="Wave Number", y="Absorbance", title="Overlayed Spectra") + theme(plot.title = element_text(size=17), axis.title.x=element_text(size=17), axis.title.y=element_text(size=17)) + geom_line(aes(dots$xs.eval, spectra.abs.ls[[which.true.source]][dots$ix,2]), size=0.25, colour="navy") + geom_line(aes(dots$xs.eval, spectra.abs.ls[[which.true.source]][dots$ix,3]), size=0.25, colour="navy") + geom_line(aes(dots$xs.eval, spectra.abs.ls[[which.true.source]][dots$ix,4]), size=0.25, colour="navy") + geom_line(aes(dots$xs.eval, spectra.abs.ls[[which.true.source]][dots$ix,5]), size=0.25, colour="navy") + geom_line(aes(dots$xs.eval, spectra.abs.ls[[which.true.source]][dots$ix,6]), size=0.25, colour="navy") + geom_line(aes(dots$xs.eval, spectra.abs.ls[[which.true.source]][dots$ix,7]), size=0.25, colour="navy")

# Pseudo spectra
pl.spec <- pl.spec + geom_line(aes(dots$xs.eval, spectra.sims.ls[[which.H1.source]][,1]), colour="grey55", size=0.2, linetype="twodash") + geom_line(aes(dots$xs.eval, spectra.sims.ls[[which.H1.source]][,2]), colour="grey55", size=0.25, linetype="twodash") + geom_line(aes(dots$xs.eval, spectra.sims.ls[[which.H1.source]][,3]), colour="grey55", size=0.25, linetype="twodash") + geom_line(aes(dots$xs.eval, spectra.sims.ls[[which.H1.source]][,4]), colour="grey55", size=0.25, linetype="twodash") + geom_line(aes(dots$xs.eval, spectra.sims.ls[[which.H1.source]][,5]), colour="grey55", size=0.25, linetype="twodash") + geom_line(aes(dots$xs.eval, spectra.sims.ls[[which.H1.source]][,6]), colour="grey55", size=0.25, linetype="twodash") + geom_line(aes(dots$xs.eval, spectra.sims.ls[[which.H1.source]][,7]), colour="grey55", size=0.25, linetype="twodash")

### Create plot
pl.spec

# ### Get plot of just real spectra
# pl.spec.real <- ggplot() + geom_line(aes(dots$xs.eval, spectra.abs.ls[[which.true.source]][dots$ix,1]),size=0.25, colour="navy") + labs(x="Wave Number", y="Absorbance", title="Real Spectra") + theme(plot.title = element_text(size=17), axis.title.x=element_text(size=17), axis.title.y=element_text(size=17)) + geom_line(aes(dots$xs.eval, spectra.abs.ls[[which.true.source]][dots$ix,2]), size=0.25, colour="navy") + geom_line(aes(dots$xs.eval, spectra.abs.ls[[which.true.source]][dots$ix,3]), size=0.25, colour="navy") + geom_line(aes(dots$xs.eval, spectra.abs.ls[[which.true.source]][dots$ix,4]), size=0.25, colour="navy") + geom_line(aes(dots$xs.eval, spectra.abs.ls[[which.true.source]][dots$ix,5]), size=0.25, colour="navy") + geom_line(aes(dots$xs.eval, spectra.abs.ls[[which.true.source]][dots$ix,6]), size=0.25, colour="navy") + geom_line(aes(dots$xs.eval, spectra.abs.ls[[which.true.source]][dots$ix,7]), size=0.25, colour="navy")
#
# ### Get plot of just pseudo spectra
# pl.spec.pseudo <- ggplot() + geom_line(aes(dots$xs.eval, spectra.sims.ls[[which.H1.source]][,1]), colour="grey15", size=0.2, linetype="twodash") + labs(x="Wave Number", y="Absorbance", title="Pseudo Spectra") + theme(plot.title = element_text(size=17), axis.title.x=element_text(size=17), axis.title.y=element_text(size=17)) + geom_line(aes(dots$xs.eval, spectra.sims.ls[[which.H1.source]][,2]), colour="grey15", size=0.25, linetype="twodash") + geom_line(aes(dots$xs.eval, spectra.sims.ls[[which.H1.source]][,3]), colour="grey15", size=0.25, linetype="twodash") + geom_line(aes(dots$xs.eval, spectra.sims.ls[[which.H1.source]][,4]), colour="grey15", size=0.25, linetype="twodash") + geom_line(aes(dots$xs.eval, spectra.sims.ls[[which.H1.source]][,5]), colour="grey15", size=0.25, linetype="twodash") + geom_line(aes(dots$xs.eval, spectra.sims.ls[[which.H1.source]][,6]), colour="grey15", size=0.25, linetype="twodash") + geom_line(aes(dots$xs.eval, spectra.sims.ls[[which.H1.source]][,7]), colour="grey15", size=0.25, linetype="twodash")
#
# ### Combine plots in a single frame
# marrangeGrob(list(pl.spec.real, pl.spec.pseudo, pl.spec), nrow=3, ncol=1, top=FALSE)


#####  Figure 2 - True spectra versus reduced spectra of two easily confused sources (19 and 34)
### Define which sources we would like to plot
which.true.source <- 19
which.H1.source <- 34

### Get average spectra
avg.spectra.H0 <- rowMeans(spectra.abs.ls[[which.true.source]][dots$ix, ])
avg.spectra.H1 <- rowMeans(spectra.abs.ls[[which.H1.source]][dots$ix, ])

### Get indices of x-axis values that are maintained for the average spectra of H0 source and H1 source
ix.pairs <- index.pair.fun(avg.spectra.H0, avg.spectra.H1, wav.num=dots$xs.eval)

### Start plot
pl.spec <- ggplot()

### Add reduced true source spectra
pl.spec <- pl.spec + labs(x="Wave Number", y="Absorbabce", title="Filtered Spectra (Downsampled to 1000 points; Union of points above some threshold)") + theme(plot.title = element_text(size=17), axis.title.x=element_text(size=17), axis.title.y=element_text(size=17)) + geom_line(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.true.source]][dots$ix[ix.pairs], 1]), size=0.1, colour="grey75") + geom_point(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.true.source]][dots$ix[ix.pairs],1]), size=0.7, colour="navy") + geom_line(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.true.source]][dots$ix[ix.pairs], 2]), size=0.1, colour="grey75") + geom_point(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.true.source]][dots$ix[ix.pairs], 2]), size=0.7, colour="navy") + geom_line(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.true.source]][dots$ix[ix.pairs], 3]), size=0.1, colour="grey75") + geom_point(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.true.source]][dots$ix[ix.pairs], 3]), size=0.7, colour="navy") + geom_line(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.true.source]][dots$ix[ix.pairs], 4]), size=0.1, colour="grey75") + geom_point(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.true.source]][dots$ix[ix.pairs], 4]), size=0.7, colour="navy") + geom_line(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.true.source]][dots$ix[ix.pairs], 5]), size=0.1, colour="grey75") + geom_point(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.true.source]][dots$ix[ix.pairs], 5]), size=0.7, colour="navy") + geom_line(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.true.source]][dots$ix[ix.pairs], 6]), size=0.1, colour="grey75") + geom_point(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.true.source]][dots$ix[ix.pairs], 6]), size=0.7, colour="navy") + geom_line(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.true.source]][dots$ix[ix.pairs], 7]), size=0.1, colour="grey75") + geom_point(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.true.source]][dots$ix[ix.pairs], 7]), size=0.7, colour="navy")

### Add pseudo spectra
pl.spec <- pl.spec + geom_line(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.H1.source]][dots$ix[ix.pairs], 1]), size=0.1, colour="grey75") + geom_point(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.H1.source]][dots$ix[ix.pairs],1]), size=0.7, colour="slategray3") + geom_line(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.H1.source]][dots$ix[ix.pairs], 2]), size=0.1, colour="grey75") + geom_point(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.H1.source]][dots$ix[ix.pairs], 2]), size=0.7, colour="slategray3") + geom_line(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.H1.source]][dots$ix[ix.pairs], 3]), size=0.1, colour="grey75") + geom_point(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.H1.source]][dots$ix[ix.pairs], 3]), size=0.7, colour="slategray3") + geom_line(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.H1.source]][dots$ix[ix.pairs], 4]), size=0.1, colour="grey75") + geom_point(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.H1.source]][dots$ix[ix.pairs], 4]), size=0.7, colour="slategray3") + geom_line(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.H1.source]][dots$ix[ix.pairs], 5]), size=0.1, colour="grey75") + geom_point(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.H1.source]][dots$ix[ix.pairs], 5]), size=0.7, colour="slategray3") + geom_line(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.H1.source]][dots$ix[ix.pairs], 6]), size=0.1, colour="grey75") + geom_point(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.H1.source]][dots$ix[ix.pairs], 6]), size=0.7, colour="slategray3") + geom_line(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.H1.source]][dots$ix[ix.pairs], 7]), size=0.1, colour="grey75") + geom_point(aes(dots$wav.num[dots$ix[ix.pairs]], spectra.abs.ls[[which.H1.source]][dots$ix[ix.pairs], 7]), size=0.7, colour="slategray3")

### Create plot
pl.spec
###################################################################################################################

###################################################################################################################
###                                                End Document                                                 ###
###################################################################################################################
