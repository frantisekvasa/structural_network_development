# The code below reproduces most analyses and figures (with the exception of certain supplementary analyses)
# conducted in the manuscript "Adolescent tuning of association cortex in human structural brain networks" by Vasa et al.
# For details regarding the motivation behind analyses and the interpretation of results, see the manuscript.
#
# written by Frantisek Vasa (fv247@cam.ac.uk), June 2015 - April 2017
 
library(vows)         # spline fitting
library(fields)       # matrix plots (image.plot)
library(matrixStats)  # rowMedians
library(seqinr)       # col2alpha

# dependencies
# rp.main - function for formatting of r-squared and p-value within titles (separate due to bugs in interaction of rstudio and the b-quote function)
source('rp.main.R') 

# set path
dir.path = '~/R/str_net_dev/'

load(paste(dir.path,'mri.data.RData',sep='')) 
# ct            cortical thickness [mm] (297 participants x 308 regions of interest)
# mt            cortical magnetization transfer [PU = percentage units] - a measure of myelination *at 70% cortical depth* (297 participants x 308 regions of interest)
# mt.depth      cortical magnetization transfer [PU = percentage units] - a measure of myelination *at 13 cortical/white matter depths* (13 depths x 297 participants x 308 regions of interest)
# age           participant age [years]
# male          participant gender [0 = female, 1 = male]
# roi.nm        region of interest name (subparcellation of desikan-killiany atlas into 308 regions with approximately equal surface area)
# roi.coord     region coordinates (x,y,z)
# nspn.age.bin  age stratum (bin), from 1-5, of subjects the NSPN study (1: 14-15 y, 2: 16-17 y, 3: 18-19 y, 4: 20-21 y, 5: 22-14 y - inclusive)
# mod.id        node affiliation to 7 modules of the age-invariant network (obtaine using the Louvain algorithm to minimise node versatility)
# ve.id         node affiliation to 7 cytoarchitectonic classes of the von Economo atlas (von Economo and Koskinas, 1925)
# yeo.id        node affiliation to 7 resting-state fMRI networks (Yeo, Krienen et al. J. Neurophysiol. 2011)

ns = length(age)    # number of subjects
nroi = dim(ct)[2]   # number of regions of interest
triup = upper.tri(matrix(nrow=nroi,ncol=nroi))  # upper triangular mask (for matrix averaging)
nboot = 1000        # number of bootstraps

# age-invariant structural correlation (SC) network ---------------------------------------
if (!dir.exists(paste(dir.path,'age_inv',sep=''))) dir.create(paste(dir.path,'age_inv',sep=''),showWarning=T, recursive=F) # directory for age-invariant results

# empirical age-invariant SC
sc = cor(ct) 

# bootstrap age-invariant SC
sc.boot = array(NA,dim=c(nroi,nroi,nboot)) 
for (i in 1:nboot) {
  if (i%%100 == 0) print(i)                 # keep track of progress
  boot.id = sample(1:ns,size=ns,replace=T)  # sample participants w/ replacement
  sc.boot[,,i] = cor(ct[boot.id,])          # bootstrapped SC
}
boot.p = 1-sign(sc)*abs(apply(sign(sc.boot),c(1,2),sum)/nboot) # bootstrap p-value - what proportion of bootstrap signs are inconsistent with empirical value?
boot.p.fdr = matrix(p.adjust(boot.p, method = "fdr"),nrow=nroi,byrow=T) # FDR correct p-values
sc.thr = array(NA,dim=c(nroi,nroi)) # initialise thresholded matrix
sc.thr[boot.p.fdr<0.01] = sc[boot.p.fdr<0.01] # set alpha; here alphaFDR = 0.01
sc.thr[as.logical(diag(nroi))] = NA; # set diagonal to NA (self-correlations to be ignored)

d = sum(!is.na(sc.thr[triup]))/((nroi*(nroi-1))/2) # edge density of age-invariant network

# matrix properties
str = rowMeans(sc)              # node strength (in unthresholded network)
deg = rowSums(!is.na(sc.thr))   # degree (number of retained edges)
wdeg = rowSums(sc.thr,na.rm=T)  # weighted degree (sum of weights of retained edges)

# histograms of strength, degree, weighted degree
# strength
pdf(paste(dir.path,'age_inv/str_dist.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
hist(str,50, xlab='strength', ylab='frequency',main='',col='grey90')
dev.off()

# degree
pdf(paste(dir.path,'age_inv/deg_dist.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
hist(deg,50, xlab='degree', ylab='frequency',main='',col='grey90')
dev.off()

# weighted degree
pdf(paste(dir.path,'age_inv/wdeg_dist.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
hist(wdeg,50, xlab='weighted degree', ylab='frequency',main='',col='grey90')
dev.off()

# export text files for surface plotting with pysurfer
if (!dir.exists(paste(dir.path,'age_inv/surf_txt',sep=''))) dir.create(paste(dir.path,'age_inv/surf_txt',sep=''),showWarning=T, recursive=F)
write(str, file = paste(dir.path,'age_inv/surf_txt/str.txt',sep=''), ncolumns = 1)
write(deg, file = paste(dir.path,'age_inv/surf_txt/deg.txt',sep=''), ncolumns = 1)
write(wdeg, file = paste(dir.path,'age_inv/surf_txt/wdeg.txt',sep=''), ncolumns = 1)

# age-dependent (varying) structural correlation (SC) network -----------------------
if (!dir.exists(paste(dir.path,'age_var',sep=''))) dir.create(paste(dir.path,'age_var',sep=''),showWarning=T, recursive=F) # directory for age-varying results

# set-up age-windows (bins)

# 1) half-interpolation of 5 NSPN-specific age-bins into nine bins
# interpolate five NSPN age-bins to create nine
nbin = 9
bin.id.all = list() # store ids of subjects in each bin
for (b in 1:nbin) {
  if (!(b%%2==0)) { # odd numbers correspond to existing NSPN bins
    bin.id = which(nspn.age.bin==((b+1)/2))
    bin.id.all[[b]] = bin.id
  } else {          # even numbers correspond to interpolated bins (parameters 29,30 ensure 60 subjects in interpolated bin)
    l.b = which(nspn.age.bin==(b/2))[(length(bin.id)-29):length(bin.id)]
    u.b = which(nspn.age.bin==((b+2)/2))[1:30]
    bin.id = c(l.b,u.b)
    bin.id.all[[b]] = bin.id
  }
}
ns.bin = unlist(lapply(bin.id.all,length)) # number of subjects in each bin

# # 2) more general sliding window analysis
# # parameters:
# ww = 40 # window width - number of participants per window
# sts = 5 # step-size - number of participants per window "step"
#     
# id.str = paste('w',toString(ww),'_s',toString(sts),sep='') # id-string for window parameters
#     
# # calculate number of age windows (bins) based on parameters ww and sts
# if ((ns-ww)%%sts > 0) { 
#   nbin = ceiling((ns-ww)/sts); 
# } else {
#   nbin = ceiling((ns-ww)/sts)+1; # unless the parameters align, the last window will contain fewer participants
# }
# 
# # obtain participant id's for general case of any combination of window width and step size
# bin.id.all = list()
# for (b in 1:nbin) {
#   bin.id = (1+sts*(b-1)):(ww+sts*(b-1))
#   bin.id.all[[b]] = bin.id
# }
# ns.bin = unlist(lapply(bin.id.all,length)) # number of subjects in each bin

# estimate empirical and bootstrapped SC within each age bin

# initialise variables
age.bin = vector(length=nbin)             # median age of participants in each bin
age.bin.b = array(NA,dim=c(nbin,nboot))   # bootstrap estimate of median age 
sc.bin = sc.bin.thr = array(NA,dim=c(nroi,nroi,nbin))   # "raw" and thresholded correlation matrices
sc.bin.b = array(NA,dim=c(nroi,nroi,nbin,nboot))        # bootstrapped correlation matrices (can be very large - ~6Gb file for 9 bins and 1000 bootstraps)
sc.bin.m = vector(length = nbin)          # mean of structural correlation in each bin
sc.bin.b.m = array(NA,dim=c(nbin,nboot))  # mean of bootstraps per bin

boot.id.all = list() # store ID's of subjects included in each bootstrap

for (b in 1:nbin) {
  print(paste('age bin',b,'of',nbin)) # keep track of age-bin progress
  
  bin.id = bin.id.all[[b]]
  age.bin[b] = median(age[bin.id])
  
  # matrices
  sc.bin[,,b] = cor(ct[bin.id,])    # empirical SC
  
  # means (of upper triangular parts) of empirical matrices
  sc.bin.m[b] = mean(sc.bin[,,b][triup])
  
  # bootstrap (re-estimation of correlation upon sampling of subjects with replacement)
  temp.id = array(NA,dim=c(nboot,ns.bin[b])) # temporarily store ID's of resampled participants
  for (i in 1:nboot) {
    if (i%%100 == 0) print(paste('bootstrap',i,'of',nboot)) # keep track of bootstrap progress
    
    boot.id = sample(1:ns.bin[b],size=ns.bin[b],replace=T) # current bootstrap sample
    temp.id[i,] = boot.id # store all bootstraps of current window
    
    age.bin.b[b,i] = median(age[bin.id[boot.id]]) # median age of bootstrapped participants
    
    sc.bin.b[,,b,i] = cor(ct[bin.id[boot.id],])     # bootstrapped SC
    
    # means (of upper triangular parts) of bootstrapped matrices
    sc.bin.b.m[b,i] = mean(sc.bin.b[,,b,i][triup])
  }
  boot.id.all[[b]] = temp.id # store bootstrapped IDs for all subjects in window
  
  # use bootstraps to threshold correlation matrices
  boot.p = 1-abs(apply(sign(sc.bin.b[,,b,]),c(1,2),sum)/nboot) # bootstrap p-value (how consistent is the sign of correlation across bootstraps?)
  boot.p.fdr = matrix(p.adjust(boot.p, method = "fdr"),nrow=nroi,byrow=T) # FDR correct p-values
  temp.thr = array(NA,dim=c(nroi,nroi))                     # initialise window matrix
  temp.thr[boot.p.fdr<0.01] = sc.bin[,,b][boot.p.fdr<0.01]  # set alpha; here alphaFDR = 0.01
  temp.thr[as.logical(diag(nroi))] = NA                     # set diagonal to NA
  sc.bin.thr[,,b] = temp.thr                                # assign thresholded matrix
  
} 
rm(temp.id,boot.p,boot.p.fdr,temp.thr) # remove temporary variables

### calculate nodal network measures
# nodal strength and degree
str.bin = apply(sc.bin,c(2,3),mean) 
deg.bin = apply(!is.na(sc.bin.thr),c(2,3),sum) 

# construct matrix of euclidean distances between nodes using coordinates of node centroids
dist.mat = array(NA,dim=c(nroi,nroi))
for (i in 1:nroi) { for (j in 1:nroi) { dist.mat[i,j] = sqrt(sum((roi.coord[i,]-roi.coord[j,])^2)) } }

# nodal distances in bootstrap thresholded network
dist.bin = array(NA,dim=c(nroi,nroi,nbin))  # distance matrix per window (bin)
dist.nod.bin = array(NA,dim=c(nroi,nbin))   # nodal average distance per window (bin)
for (b in 1:nbin) {
  temp.dist = array(NA,dim=c(nroi,nroi))          # temporary distance matrix
  temp.dist[!is.na(sc.bin.thr[,,b])] = dist.mat[!is.na(sc.bin.thr[,,b])] # select only retained edges
  dist.bin[,,b] = temp.dist                       # distance matrix 
  dist.nod.bin[,b] = rowMeans(temp.dist,na.rm=T)  # nodal distance
}
rm(temp.dist) # clear temporary variables

# calculate graph measures of organisation for the thresholded networks
dens = dist.thr = vector(length=nbin)
for (b in 1:nbin) {
  dens[b] = sum(!is.na(sc.bin.thr[,,b][triup]))/((nroi*(nroi-1))/2)       # edge density
  temp.dist = array(NA,dim=c(nroi,nroi))
  temp.dist[!is.na(sc.bin.thr[,,b])] = dist.mat[!is.na(sc.bin.thr[,,b])]
  dist.thr[b] = mean(temp.dist[triup],na.rm=T)                            # average euclidean distance (spanned by retained edges)
}

save('nbin','bin.id.all','ns.bin','age.bin','age.bin.b','sc.bin','sc.bin.thr','sc.bin.m','sc.bin.b.m','boot.id.all',
     'str.bin','deg.bin','dist.mat','dist.bin','dist.nod.bin','dens','dist.thr',file=paste(dir.path,'age_var/sc.age.var.RData',sep=''))
#'cov.bin','var.bin','cov.bin.b.m','var.bin.b.m','cov.bin.m','var.bin.m', # optional covariance and product of variances

# number of male and female participants per bin
nm.bin.all = nf.bin.all = vector(length=nbin)
for (b in 1:nbin) {
  nm.bin.all[b] = sum(male[bin.id.all[[b]]])  # number male
  nf.bin.all[b] = sum(!male[bin.id.all[[b]]]) # number female
}

# do any negative entries remain in the thresholded matrix?
# sum(sc.bin.thr<0)/2 # divide by two as every edge is counted twice (within both upper and lower triangular parts)
# sc.bin.thr[sc.bin.thr<0] = 0 # optionally, negative edges can be set to zero 

# plot correlation distributions ------------------------------------------
if (!dir.exists(paste(dir.path,'age_var/dist',sep=''))) dir.create(paste(dir.path,'age_var/dist',sep=''),showWarning=T, recursive=F)
dens.step = 1:nbin # steps at which densities are plotted (in this case for all age-bins)
rbPal <- colorRampPalette(c('blue','red')); colBar = rbPal(100) # colorbar set-up

# map age to colobar
col.l = floor(min(age.bin))   # lower limit of colorbar
col.u = ceiling(max(age.bin)) # upper limit of colorbar
age.temp = ceiling(100*((age.bin-col.l)/(col.u-col.l))) # mapping of age to colorbar (on scale from 1 to 100)
dens.col = colBar[age.temp[dens.step]] # color per age-bin
dens.col.alph = sapply(dens.col, col2alpha, 0.1) # transparent colors (tranparency parameter alpha = 0.1)

dens.n = 200 # number of points at which density is estimated

# set up plot parameters
# dens.i / dens.f = initial / final "cut-off" values (between which values of distribution are estimated)
# dens.thr.i / dens.thr.f = as above, for bootstrap-thresholded distributions
# xlab = x-axis label
# dist.sc / dist.sc.thr = parameter to scale size of distibutions (controls overlap across age-bins)
# abl / abl.thr = positions of "guiding" vertical dashed lines on x-axis

# structural correlation
dens.i = -0.5; dens.f = 1 # sc
dens.thr.i = 0.25; dens.thr.f = 0.75 # sc.thr
xlab = expression(paste('correlation'))
dist.sc = 1; dist.sc.thr = 0.3 # cor
abl = c(0,0.5); abl.thr = c(0.4,0.6) # cor

# general
dens.x = seq(from=dens.i,to=dens.f,length.out = dens.n)
dens.y = matrix(NA, nrow=length(dens.step), ncol=dens.n) 
dens.thr.x = seq(from=dens.thr.i,to=dens.thr.f,length.out = dens.n)
dens.thr.y = matrix(NA, nrow=length(dens.step), ncol=dens.n) 
dens.boot.y = array(NA,dim=c(length(dens.step),nboot,dens.n))

# changes in weight (+ distance) distribution
for (b in 1:length(dens.step)) {
  print(paste('age bin',b,'of',nbin)) # keep track of age-bin progress
  dens.y[b,] = density(sc.bin[,,dens.step[b]][triup],n=dens.n,from=dens.i,to=dens.f)$y                  # full (unthresholded)
  for (i in 1:nboot) { # bootstrap 
    if (i%%100 == 0) print(paste('bootstrap',i,'of',nboot)) # keep track of bootstrap progress
    dens.boot.y[b,i,] = density(sc.bin.b[,,dens.step[b],i][triup],n=dens.n,from=dens.i,to=dens.f)$y 
    }
  dens.thr.y[b,] = density(sc.bin.thr[,,dens.step[b]][triup],n=dens.n,from=dens.thr.i,to=dens.thr.f,na.rm=T)$y  # thresholded
}

# plot empirical distributions
pdf(paste(dir.path,'age_var/dist/sc_dens_line.pdf',sep=''),width=5,height=7)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
plot(0, yaxt='n', type='n', ylab='median age (years)', xlab=xlab, ylim=rev(c(14-2,24)), xlim=c(dens.i,dens.f))
abline(v=abl,col='grey',lty=2)
axis(2, at=seq(14,24,by=2))
for (b in 1:nbin) {
  abline(a=age.bin[b],b=0,col='grey')
  lines(dens.x,age.bin[b]-dist.sc*dens.y[b,],col=dens.col[b],lwd=2)
}
dev.off()

# plot bootstrapped distributions
pdf(paste(dir.path,'age_var/dist/sc_dens_boot.pdf',sep=''),width=5,height=7)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
plot(0, yaxt='n', type='n', ylab='median age (years)', xlab=xlab, ylim=rev(c(14-2,24)), xlim=c(dens.i,dens.f))
abline(v=abl,col='grey',lty=2)
axis(2, at=seq(14,24,by=2))
for (b in 1:nbin) {
  for (i in 1:nboot) {
    lines(dens.x,age.bin[b]-dist.sc*dens.boot.y[b,i,],col=dens.col.alph[b])
  }
  lines(dens.x,age.bin[b]-dist.sc*colMeans(dens.boot.y[b,,]),col='white',lwd=2,lty=1)
  abline(a=age.bin[b],b=0,col='grey')
}
dev.off()

# plot thresholded distributions
pdf(paste(dir.path,'age_var/dist/sc_thr_dens_line.pdf',sep=''),width=5,height=7)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
plot(0, yaxt='n', type='n', ylab='median age (years)', xlab=xlab, ylim=rev(c(14-2,24)), xlim=c(dens.thr.i,dens.thr.f))
abline(v=abl.thr,col='grey',lty=2)
axis(2, at=seq(14,24,by=2))
for (b in 1:nbin) {
  abline(a=age.bin[b],b=0,col='grey')
  lines(dens.thr.x,age.bin[b]-dist.sc.thr*dens.thr.y[b,],col=dens.col[b],lwd=2)
}
dev.off()

# plot colorbar
color.bar <- function(lut, alpha, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  lut = sapply(lut, col2alpha, alpha)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, at=ticks, labels=rev(ticks), las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
pdf(paste(dir.path,'age_var/dist/bin_wei_dens_colorbar.pdf',sep=''),width=1.5,height=6)
color.bar(colorRampPalette(c('red','blue'))(100), alpha = 0.4, min = 15, max =  23, nticks = 5, ticks = rev(seq(15, 23, len=5)))
dev.off()

# fit linear models to global and nodal measures --------------------------
if (!dir.exists(paste(dir.path,'age_var/traj',sep=''))) dir.create(paste(dir.path,'age_var/traj',sep=''),showWarning=T, recursive=F) # folder for storing trajectories

### empirical results
# correlation
l = lm(sc.bin.m~age.bin)
sc.bin.m.l.aic = AIC(l) # AIC
summary(l) # statistics - r-squared and p-value

# edge density
l = lm(dens~age.bin)
dens.l.aic = AIC(l) # AIC
summary(l) # statistics - r-squared and p-value

# distance
l = lm(dist.thr~age.bin)
dist.glob.l.aic = AIC(l) # AIC
summary(l) # statistics - r-squared and p-value

### calculate bootstrap medians and confidence intervals for global measures
# bm = bootstrap mean
# bqrt = bootstrap inter-quartile range (corresponding to [0.25,0.75] quantile)

# bootstrap - global trajectories
sc.m.bm = rowMedians(sc.bin.b.m); 
sc.m.bqrt = apply(sc.bin.b.m, 1, quantile, probs = c(0.25,0.75), na.rm = T)
# bootstrap - age
age.bm = rowMedians(age.bin.b); 
age.bqrt = apply(age.bin.b, 1, quantile, probs = c(0.25,0.75), na.rm = T)

# fit linear models for each bootstrap (using bootstrap-specific median ages)
sc.b.m.sl = sc.b.m.int = sc.b.m.l.aic = vector(length = nboot) # correlation
for (i in 1:nboot) {
  if (i%%100 == 0) print(paste('bootstrap',i,'of',nboot)) # keep track of bootstrap progress
  l = lm(sc.bin.b.m[,i] ~ age.bin.b[,i]); sc.b.m.int[i] = l$coefficients[1]; sc.b.m.sl[i] = l$coefficients[2]; sc.b.m.l.aic[i] = AIC(l)       # correlations
}

# fit linear models to nodal: strength, degree, distance (for all values except sl and int, l indicates linear model)
# sl = slope
# int = intercept
# p = p-value
# rsq = r-squared (variance explained by fit)
# rsq.adj = adjusted r-squared
# aic = akaike's information criterion
# bic = bayesian information criterion
str.sl = str.int = str.l.p = str.l.rsq = str.l.rsq.adj = str.l.aic = str.l.bic = vector(length = nroi)        # strength
deg.sl = deg.int = deg.l.p = deg.l.rsq = deg.l.rsq.adj = deg.l.aic = deg.l.bic = vector(length = nroi)        # degree
dist.sl = dist.int = dist.l.p = dist.l.rsq = dist.l.rsq.adj = dist.l.aic = dist.l.bic = vector(length = nroi) # nodal distance
for (n in 1:nroi) {
  # str
  l = lm(str.bin[n,] ~ age.bin)   # fit linear model
  str.int[n] = l$coefficients[1]  # intercept
  str.sl[n] = l$coefficients[2]   # slope
  f = summary(l)$fstatistic; str.l.p[n] = pf(f[1],f[2],f[3], lower.tail = FALSE)    # p-value
  str.l.rsq[n] = summary(l)$r.squared; str.l.rsq.adj[n] = summary(l)$adj.r.squared  # r-squared (raw and adjusted)
  str.l.aic[n] = AIC(l) # aic
  str.l.bic[n] = BIC(l) # bic
  # deg
  l = lm(deg.bin[n,] ~ age.bin)   # fit linear model
  deg.int[n] = l$coefficients[1]  # intercept
  deg.sl[n] = l$coefficients[2]   # slope
  f = summary(l)$fstatistic; deg.l.p[n] = pf(f[1],f[2],f[3], lower.tail = FALSE)    # p-value
  deg.l.rsq[n] = summary(l)$r.squared; deg.l.rsq.adj[n] = summary(l)$adj.r.squared  # r-squared (raw and adjusted)
  deg.l.aic[n] = AIC(l) # aic
  deg.l.bic[n] = BIC(l) # bic
  # dist
  l = lm(dist.nod.bin[n,] ~ age.bin)  # fit linear model
  dist.int[n] = l$coefficients[1]     # intercept
  dist.sl[n] = l$coefficients[2]      # slope
  f = summary(l)$fstatistic; dist.l.p[n] = pf(f[1],f[2],f[3], lower.tail = FALSE)    # p-value 
  dist.l.rsq[n] = summary(l)$r.squared; dist.l.rsq.adj[n] = summary(l)$adj.r.squared # r-squared (raw and adjusted)
  dist.l.aic[n] = AIC(l) # aic
  dist.l.bic[n] = BIC(l) # bic
}
# # optionally write out files for plotting of slopes
# if (!dir.exists(paste(dir.path,'age_var/surf_txt',sep=''))) dir.create(paste(dir.path,'age_var/surf_txt',sep=''),showWarning=T, recursive=F)
# write(str.sl, file = paste(dir.path,'age_var/surf_txt/str_sl.txt',sep=''), ncolumns = 1)
# write(deg.sl, file = paste(dir.path,'age_var/surf_txt/deg_sl.txt',sep=''), ncolumns = 1)
# write(dist.sl, file = paste(dir.path,'age_var/surf_txt/dist_sl.txt',sep=''), ncolumns = 1)

save('str.bin','str.sl','str.int','str.l.p','str.l.rsq','str.l.rsq.adj','str.l.aic',
     'deg.bin','deg.sl','deg.int','deg.l.p','deg.l.rsq','deg.l.rsq.adj','deg.l.aic',
     'dist.nod.bin','dist.sl','dist.int','dist.l.p','dist.l.rsq','dist.l.rsq.adj','dist.l.aic',
     'sc.bin.m.l.aic','sc.m.bm','sc.m.bqrt','age.bm','sc.b.m.sl','sc.b.m.int','sc.b.m.l.aic',
     file=paste(dir.path,'age_var/sc.age.var.stats.lin.RData',sep=''))

# fit smoothing splines to global and nodal measures ---------------------------------------------------

# spline parameters
kspl = 6        # number of basis functions
min.sp = 0.065  # minimum log smoothing parameter (for 9 windows, correspond to max ~3.5 degrees of freedom)
npred = 100     # number of predicted time-points (for plots of smoothing spline)
age.pred = seq(from=14,to=25,length.out=100) # predicted ages (for plots of smoothing spline)
age.min.pred = seq(from=14,to=25,length.out=1000) # ages for calculation of age at minimum (greater precision than for plotting)

### empirical results

# correlation
g = gam(sc.bin.m ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML')
sc.bin.m.g.aic = AIC(g) # AIC
summary(g)              # statistics - r-squared and p-value
# age at minimum correlation
pred = predict(g,newdata=data.frame('x'=age.min.pred),se.fit=T)
age.min.pred[which(pred$fit==min(pred$fit))]

# edge density
g = gam(dens ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML')
dens.g.aic = AIC(g)     # AIC
summary(g)              # statistics - r-squared and p-value
# age at minimum density
pred = predict(g,newdata=data.frame('x'=age.min.pred),se.fit=T)
age.min.pred[which(pred$fit==min(pred$fit))] # age at minimum

# euclidean distance
g = gam(dist.thr ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML')
dist.glob.g.aic = AIC(g)  # AIC
summary(g)                # statistics - r-squared and p-value
# age at minimum distance
pred = predict(g,newdata=data.frame('x'=age.min.pred),se.fit=T)
age.min.pred[which(pred$fit==min(pred$fit))] # age at minimum

# fit splines to global measures across bootstraps (using bootstrap-specific median ages)
# df = degrees of freedom
# aic = akaike's information criterion
# sc.b.m.g = store predicted smoothing spline curve for later plotting
sc.b.m.df = sc.b.m.g.aic = vector(length=nboot); 
sc.b.m.g = array(NA,dim=c(nboot,npred)) 
for (i in 1:nboot) {
  if (i%%100 == 0) print(paste('bootstrap',i,'of',nboot)) # keep track of bootstrap progress
  x = age.bin.b[,i] # bootstrap age (needs to be assigned to "x" variable for ease during prediction)
  # correlation
  g = gam(sc.bin.b.m[,i] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML') # fit smoothing spline
  sc.b.m.g[i,] = predict(g,newdata=data.frame('x'=age.pred),se.fit=F)         # predicted values (for plotting of trajectory)
  sc.b.m.df[i] = sum(g$edf); sc.b.m.g.aic[i] = AIC(g)                         # average degrees of freedom and AIC
}

# sc - mean
pdf(paste(dir.path,'age_var/traj/sc_bin_mean_spline.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
plot(0, type='n', xlab='median age (years)', ylab=expression(paste(mu, ' correlation')), xlim=c(14,25), ylim=c(0.1,0.5)) # set-up plot
for (i in 1:nboot) lines(age.pred,sc.b.m.g[i,],col='grey',lwd=1)  # plot bootstrap estimates of trajectories
lines(age.pred,colMeans(sc.b.m.g),col='white',lty=2,lwd=5)        # plot means of bootsrapped estimates
points(age.bin,sc.bin.m,col='black',pch=19,cex=0.75)              # plot empirical data 
for (b in 1:nbin) { lines(c(age.bm[b],age.bm[b]),sc.m.bqrt[,b],col='black'); lines(age.bqrt[,b],c(sc.m.bm[b],sc.m.bm[b]),col='black') } # plot inter-quartile ranges (for data (y-axis) and ages(x-axis))
x = age.bin; g = gam(sc.bin.m ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML')                    # re-fit smoothing spline for plotting
pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=F); lines(age.pred,pred,col='black',lwd=1.5) # empirical smoothing spline
points(rep(age.pred[which(pred==min(pred))],2),rep(min(pred),2),col=c('black','white'),pch=c(19,19),cex=c(1,0.5)) # marker at the age at minimum
dev.off()

# edge density 
pdf(paste(dir.path,'age_var/traj/sc_thr_dens.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
plot(0, type='n', xlab='median age (years)', ylab='edge density (%)', xlim=c(14,25), ylim=c(0,50))
points(age.bin,100*dens,col='black',pch=19,cex=0.75)
x = age.bin; g = gam(100*dens ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML')
pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=T)
polygon(c(rev(age.pred), age.pred), c(rev(pred$fit+2*pred$se.fit), pred$fit-2*pred$se.fit), col = col2alpha('grey',alpha=0.5), border = NA)
lines(age.pred,pred$fit,col='black',lwd=3); 
points(rep(age.pred[which(pred$fit==min(pred$fit))],2),rep(min(pred$fit),2),col=c('black','white'),pch=c(19,19),cex=c(1,0.5))
dev.off()

# distance
pdf(paste(dir.path,'age_var/traj/sc_thr_dist.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
plot(0, type='n', xlab='median age (years)', ylab=expression(paste(mu, ' distance (mm)')), xlim=c(14,25), ylim=c(55,90))
points(age.bin,dist.thr,col='black',pch=19,cex=0.75)
x = age.bin; g = gam(dist.thr ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML')
pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=T)
polygon(c(rev(age.pred), age.pred), c(rev(pred$fit+2*pred$se.fit), pred$fit-2*pred$se.fit), col = col2alpha('grey',alpha=0.5), border = NA)
lines(age.pred,pred$fit,col='black',lwd=3); 
points(rep(age.pred[which(pred$fit==min(pred$fit))],2),rep(min(pred$fit),2),col=c('black','white'),pch=c(19,19),cex=c(1,0.5))
dev.off()

### statistics for nodal spline fits (prefix g signifies spline fit (g for generalised additive model))
# df = average degrees of freedom
# rsq.adj = adjusted r-squared
# p = p-value
# aic = akaike's information criterion
# bic = bayesian information criterion
# min.age = age at minimum
# d.max = maximum change (difference between maximum and minimum, with direction coded in sign: + = increase, - = decrease)
# pred = predicted values (for curve plotting)

# initialise predicted variables
str.df = str.g.rsq.adj = str.g.p = str.g.aic = str.g.bic = str.min.age = str.d.max = vector(length = nroi)        # strength
deg.df = deg.g.rsq.adj = deg.g.p = deg.g.aic = deg.g.bic = deg.min.age = deg.d.max = vector(length = nroi)        # degree
dist.df = dist.g.rsq.adj = dist.g.p = dist.g.aic = dist.g.bic = dist.min.age = dist.d.max = vector(length = nroi) # distance
for (n in 1:nroi) {
  if (n%%30 == 0) print(paste('region',n,'of',nroi)) # keep track of bootstrap progress
  x = age.bin # for ease of prediction
  # strength
  g = gam(str.bin[n,] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML') # fit spline
  str.df[n] = sum(g$edf)              # average degrees of freedom
  str.g.rsq.adj[n] = summary(g)$r.sq  # adjusted r-squared
  str.g.p[n] = summary(g)$s.pv        # p-value
  str.g.aic[n] = AIC(g)               # akaike's information criterion
  str.g.bic[n] = BIC(g)               # bayesian information criterion
  pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=F) # predict spline trajectory
  str.min.age[n] = age.pred[which(pred==min(pred))]   # age at minimum
  str.d.max[n] = (max(pred)-min(pred))/sign(age.pred[which(pred==max(pred))]-age.pred[which(pred==min(pred))]) # maximum change
  # degree
  g = gam(deg.bin[n,] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML') # fit spline
  deg.df[n] = sum(g$edf)              # average degrees of freedom
  deg.g.rsq.adj[n] = summary(g)$r.sq  # adjusted r-squared
  deg.g.p[n] = summary(g)$s.pv        # p-value
  deg.g.aic[n] = AIC(g)               # akaike's information criterion
  deg.g.bic[n] = BIC(g)               # bayesian information criterion
  pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=F) # predict spline trajectory
  deg.min.age[n] = age.pred[which(pred==min(pred))]   # age at minimum
  deg.d.max[n] = (max(pred)-min(pred))/sign(age.pred[which(pred==max(pred))]-age.pred[which(pred==min(pred))]) # maximum change
  # distance
  g = gam(dist.nod.bin[n,] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML') # fit spline
  dist.df[n] = sum(g$edf)              # average degrees of freedom
  dist.g.rsq.adj[n] = summary(g)$r.sq  # adjusted r-squared
  dist.g.p[n] = summary(g)$s.pv        # p-value
  dist.g.aic[n] = AIC(g)               # akaike's information criterion
  dist.g.bic[n] = BIC(g)               # bayesian information criterion
  pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=F) # predict spline trajectory
  dist.min.age[n] = age.pred[which(pred==min(pred))]   # age at minimum
  dist.d.max[n] =(max(pred)-min(pred))/sign(age.pred[which(pred==max(pred))]-age.pred[which(pred==min(pred))]) # maximum change
}

# # optionally write out files for plotting of summary statistics
# # strength
# write(str.min.age, file = paste(dir.path,'age_var/surf_txt/str_min_age.txt',sep=''), ncolumns = 1)
# write(str.d.max, file = paste(dir.path,'age_var/surf_txt/str_d_max.txt',sep=''), ncolumns = 1)
# # degree
# write(deg.min.age, file = paste(dir.path,'age_var/surf_txt/deg_min_age.txt',sep=''), ncolumns = 1)
# write(deg.d.max, file = paste(dir.path,'age_var/surf_txt/deg_d_max.txt',sep=''), ncolumns = 1)
# # nodal distance
# write(dist.min.age, file = paste(dir.path,'age_var/surf_txt/dist_min_age.txt',sep=''), ncolumns = 1)
# write(dist.d.max, file = paste(dir.path,'age_var/surf_txt/dist_d_max.txt',sep=''), ncolumns = 1)

save('kspl','min.sp','npred','age.pred','sc.bin.m.g.aic','sc.b.m.df','sc.b.m.g.aic','sc.b.m.g',
     'str.df','str.g.rsq.adj','str.g.p','str.g.aic','str.g.bic','str.min.age','str.d.max',
     'deg.df','deg.g.rsq.adj','deg.g.p','deg.g.aic','deg.g.bic','deg.min.age','deg.d.max',
     'dist.df','dist.g.rsq.adj','dist.g.p','dist.g.aic','dist.g.bic','dist.min.age','dist.d.max',
     file=paste(dir.path,'age_var/sc.age.var.stats.spline.RData',sep=''))
#'cov.bin.m.g.aic','var.bin.m.g.aic','cov.b.m.df','cov.b.m.g.aic','cov.b.m.g','var.b.m.df','var.b.m.g.aic','var.b.m.g'

### plot splines at nodes
# strength
pdf(paste(dir.path,'age_var/traj/str_bins_spline.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
par(mfrow=c(1,1))
for (n in 1:nroi) {
  if (n%%30 == 0) print(paste('region',n,'of',nroi)) 
  plot(0, type='n', xlab='median age (years)',ylab=expression(paste('strength')),font.main=1,main=paste('ROI ',toString(n),' = ',roi.nm[n],sep=''),xlim=c(14,25),ylim=c(0,0.8))
  points(age.bin,str.bin[n,],col='black',pch=19,cex=0.75);
  x = age.bin; g = gam(str.bin[n,] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML');
  pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=T)
  polygon(c(rev(age.pred), age.pred), c(rev(pred$fit+2*pred$se.fit), pred$fit-2*pred$se.fit), col = col2alpha('grey',alpha=0.5), border = NA)
  lines(age.pred,pred$fit,col='black',lwd=1.5)
}
dev.off()

# degree
pdf(paste(dir.path,'age_var/traj/deg_bins_spline.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
par(mfrow=c(1,1))
for (n in 1:nroi) {
  if (n%%30 == 0) print(paste('region',n,'of',nroi))
  plot(0, type='n', xlab='median age (years)',ylab=expression(paste('degree')),font.main=1,main=paste('ROI ',toString(n),' = ',roi.nm[n],sep=''),xlim=c(14,25),ylim=c(-10,230))
  points(age.bin,deg.bin[n,],col='black',pch=19,cex=0.75);
  x = age.bin; g = gam(deg.bin[n,] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML');
  pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=T)
  polygon(c(rev(age.pred), age.pred), c(rev(pred$fit+2*pred$se.fit), pred$fit-2*pred$se.fit), col = col2alpha('grey',alpha=0.5), border = NA)
  lines(age.pred,pred$fit,col='black',lwd=1.5)
}
dev.off()

# compare linear and spline fits ------------------------------------------

# 1) globally

# correlation
if (sc.bin.m.g.aic<sc.bin.m.l.aic) print('smoothing spline provides a better fit to mean correlation than linear model')
if (sc.bin.m.g.aic>sc.bin.m.l.aic) print('linear model provides a better fit to mean correlation than smoothing spline')
# edge density
if (dens.g.aic<dens.l.aic) print('smoothing spline provides a better fit to edge density than linear model')
if (dens.g.aic>dens.l.aic) print('linear model provides a better fit to edge density than smoothing spline')
# euclidean distance
if (dist.glob.g.aic<dist.glob.l.aic) print('smoothing spline provides a better fit to average euclidean distance than linear model')
if (dist.glob.g.aic>dist.glob.l.aic) print('linear model provides a better fit to average euclidean distance than smoothing spline')

# 2) locally

### strength
# correct model p-values for multiple comparisons
str.l.p.fdr = p.adjust(str.l.p, method = 'fdr') # linear
str.g.p.fdr = p.adjust(str.g.p, method = 'fdr') # spline

str.lin = which(str.g.aic > str.l.aic)                        # linear trajectories
length(intersect(str.lin,which(str.l.p.fdr<0.05)))            # number of significant linear trajectories
sum(str.d.max[str.lin]<0)                                     # number of linear decreases
sum(str.d.max[str.lin]>0)                                     # number of linear increases
sum(str.d.max[intersect(str.lin,which(str.l.p.fdr<0.05))]<0)  # number of significant linear decreases
sum(str.d.max[intersect(str.lin,which(str.l.p.fdr<0.05))]>0)  # number of significant linear increases

str.nonlin = which(str.g.aic < str.l.aic)                       # non-linear trajectories
length(intersect(str.nonlin,which(str.g.p.fdr<0.05)))           # number of significant nonlinear trajectories
sum(str.d.max[str.nonlin]<0)                                    # number of nonlinear decreases
sum(str.d.max[str.nonlin]>0)                                    # number of nonlinear increases
sum(str.d.max[intersect(str.nonlin,which(str.g.p.fdr<0.05))]<0) # number of significant nonlinear decreases
sum(str.d.max[intersect(str.nonlin,which(str.g.p.fdr<0.05))]>0) # number of significant nonlinear increases

# plot trajectories by four categories (colours) - lin up / down // nonlin up / down (all + sig)
traj.type.str.d.max = traj.type.str.d.max.sig = array(NA,dim=c(nroi))
# all
traj.type.str.d.max[intersect(str.lin,which(str.d.max<0))] = 1    # lin decrease
traj.type.str.d.max[intersect(str.lin,which(str.d.max>0))] = 2    # lin increase
traj.type.str.d.max[intersect(str.nonlin,which(str.d.max<0))] = 3 # nonlin decrease
traj.type.str.d.max[intersect(str.nonlin,which(str.d.max>0))] = 4 # nonlin increase
# significant (FDR-adjusted threshold alphaFDR = 0.05)
traj.type.str.d.max.sig[intersect(which(str.d.max<0),intersect(str.lin,which(str.l.p.fdr<0.05)))] = 1     # significant lin decrease
traj.type.str.d.max.sig[intersect(which(str.d.max>0),intersect(str.lin,which(str.l.p.fdr<0.05)))] = 2     # significant lin increase
traj.type.str.d.max.sig[intersect(which(str.d.max<0),intersect(str.nonlin,which(str.g.p.fdr<0.05)))] = 3  # significant nonlin decrease
traj.type.str.d.max.sig[intersect(which(str.d.max>0),intersect(str.nonlin,which(str.g.p.fdr<0.05)))] = 4  # significant nonlin increase

# # export trajectory types for plotting
# write(traj.type.str.d.max, file = paste(dir.path,'age_var/surf_txt/traj_type_str_d_max.txt',sep=''), ncolumns = 1)
# write(traj.type.str.d.max.sig, file = paste(dir.path,'age_var/surf_txt/traj_type_str_d_max_sig.txt',sep=''), ncolumns = 1)

# # plot str.d.max and str.age.min for subsets of regions with significant strength change
# str.sig = which(!is.na(traj.type.str.sig))
# str.d.max.sig = str.min.age.sig = array(NA,dim=c(nroi))
# str.d.max.sig[str.sig] = str.d.max[str.sig]
# str.min.age.sig[str.sig] = str.min.age[str.sig]
# write(str.d.max.sig, file = paste(dir.path,'age_var/surf_txt/str_d_max_sig.txt',sep=''), ncolumns = 1)
# write(str.min.age.sig, file = paste(dir.path,'age_var/surf_txt/str_min_age_sig.txt',sep=''), ncolumns = 1)

# plot trajectories, with averages, for each of four classes
col.traj.type = c('powderblue','lightcoral','steelblue','red3')

# all
spl.all = matrix(nrow = nroi, ncol = length(age.pred))
pdf(paste(dir.path,'age_var/traj/traj_type_str.pdf',sep=''),width=5.2,height=5)
par(mfrow=c(2,2),oma=c(4,4,2,2),mar=c(2,2,1,1)+0.1,cex.lab=2,cex.axis=1.5,cex.main=2,font.main=1,bg='white')
for (i in 1:4) {
  print(i)
  npred = length(age.pred)
  plot(0, type='n', xaxt='n', yaxt='n', xlim=c(min(age)-0.5,max(age)+0.5), ylim=c(0,0.6))
  axis(1,at=seq(15,25,by=2.5),labels=c(15,'',20,'',25))
  axis(2,at=seq(0,0.6,by=0.15),labels=c(0,'',0.3,'',0.6))
  for (m in which(traj.type.str.d.max==i)) {
    x = age.bin; g = gam(str.bin[m,] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML');
    pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=F); spl.all[m,] = pred
    lines(age.pred,pred,col=col.traj.type[i],lwd=1)
  }
  if (length(which(traj.type.str.d.max==i)) > 1) lines(age.pred,colMeans(spl.all[which(traj.type.str.d.max==i),]), type = 'l', lty = 1, lwd = 3, col = 'grey10') # if cluster contains more than one traj, plot mean
}
mtext('median age (years)',outer=T,side=1,line=1.5,cex=2)
mtext('strength',outer=T,side=2,line=1.5,cex=2)
mtext('decrease',outer=T,side=3,line=-0.2,cex=1.75,adj=0.18); mtext('increase',outer=T,side=3,line=-0.2,cex=1.75,adj=0.88)
mtext('nonlinear',outer=T,side=4,line=0.2,cex=1.75,adj=0.18); mtext('linear',outer=T,side=4,line=0.2,cex=1.75,adj=0.85)
dev.off()

# significant regions only
spl.all = matrix(nrow = nroi, ncol = length(age.pred))
pdf(paste(dir.path,'age_var/traj/traj_type_str_sig.pdf',sep=''),width=5.2,height=5)
par(mfrow=c(2,2),oma=c(4,4,2,2),mar=c(2,2,1,1)+0.1,cex.lab=2,cex.axis=1.5,cex.main=2,font.main=1,bg='white')
for (i in 1:4) {
  print(i)
  npred = length(age.pred)
  plot(0, type='n', xaxt='n', yaxt='n', xlim=c(min(age)-0.5,max(age)+0.5), ylim=c(0,0.6))
  axis(1,at=seq(15,25,by=2.5),labels=c(15,'',20,'',25))
  axis(2,at=seq(0,0.6,by=0.15),labels=c(0,'',0.3,'',0.6))
  for (m in which(traj.type.str.d.max.sig==i)) {
    x = age.bin; g = gam(str.bin[m,] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML');
    pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=F); spl.all[m,] = pred
    lines(age.pred,pred,col=col.traj.type[i],lwd=1)
  }
  if (length(which(traj.type.str.d.max.sig==i)) > 1) lines(age.pred,colMeans(spl.all[which(traj.type.str.d.max.sig==i),]), type = 'l', lty = 1, lwd = 3, col = 'grey10') # if cluster contains more than one traj, plot mean
}
mtext('median age (years)',outer=T,side=1,line=1.5,cex=2)
mtext('strength',outer=T,side=2,line=1.5,cex=2)
mtext('decrease',outer=T,side=3,line=-0.2,cex=1.75,adj=0.18); mtext('increase',outer=T,side=3,line=-0.2,cex=1.75,adj=0.88)
mtext('nonlinear',outer=T,side=4,line=0.2,cex=1.75,adj=0.18); mtext('linear',outer=T,side=4,line=0.2,cex=1.75,adj=0.85)
dev.off()

### degree
# correct model p-values for multiple comparisons
deg.l.p.fdr = p.adjust(deg.l.p, method = 'fdr') # linear
deg.g.p.fdr = p.adjust(deg.g.p, method = 'fdr') # spline

deg.lin = which(deg.g.aic > deg.l.aic)                        # linear trajectories
length(intersect(deg.lin,which(deg.l.p.fdr<0.05)))            # number of significant linear trajectories
sum(deg.d.max[deg.lin]<0)                                     # number of linear decreases
sum(deg.d.max[deg.lin]>0)                                     # number of linear increases
sum(deg.d.max[intersect(deg.lin,which(deg.l.p.fdr<0.05))]<0)  # number of significant linear decreases
sum(deg.d.max[intersect(deg.lin,which(deg.l.p.fdr<0.05))]>0)  # number of significant linear increases

deg.nonlin = which(deg.g.aic < deg.l.aic)                       # non-linear trajectories
length(intersect(deg.nonlin,which(deg.g.p.fdr<0.05)))           # number of significant nonlinear trajectories
sum(deg.d.max[deg.nonlin]<0)                                    # number of nonlinear decreases
sum(deg.d.max[deg.nonlin]>0)                                    # number of nonlinear increases
sum(deg.d.max[intersect(deg.nonlin,which(deg.g.p.fdr<0.05))]<0) # number of significant nonlinear decreases
sum(deg.d.max[intersect(deg.nonlin,which(deg.g.p.fdr<0.05))]>0) # number of significant nonlinear increases

# plot trajectories by four categories (colours) - lin up / down // nonlin up / down (all + sig)
traj.type.deg = traj.type.deg.sig = array(NA,dim=c(nroi))
# all
traj.type.deg[intersect(deg.lin,which(deg.d.max<0))] = 1    # lin decrease
traj.type.deg[intersect(deg.lin,which(deg.d.max>0))] = 2    # lin increase
traj.type.deg[intersect(deg.nonlin,which(deg.d.max<0))] = 3 # nonlin decrease
traj.type.deg[intersect(deg.nonlin,which(deg.d.max>0))] = 4 # nonlin increase
# significant (FDR-adjusted threshold alphaFDR = 0.05)
traj.type.deg.sig[intersect(which(deg.d.max<0),intersect(deg.lin,which(deg.l.p.fdr<0.05)))] = 1     # lin decrease
traj.type.deg.sig[intersect(which(deg.d.max>0),intersect(deg.lin,which(deg.l.p.fdr<0.05)))] = 2     # lin increase
traj.type.deg.sig[intersect(which(deg.d.max<0),intersect(deg.nonlin,which(deg.g.p.fdr<0.05)))] = 3  # nonlin decrease
traj.type.deg.sig[intersect(which(deg.d.max>0),intersect(deg.nonlin,which(deg.g.p.fdr<0.05)))] = 4  # nonlin increase

# # export trajectory types for plotting
# write(traj.type.deg, file = paste(dir.path,'age_var/surf_txt/traj_type_deg.txt',sep=''), ncolumns = 1)
# write(traj.type.deg.sig, file = paste(dir.path,'age_var/surf_txt/traj_type_deg_sig.txt',sep=''), ncolumns = 1)

# # plot deg.d.max and deg.age.min for subsets of regions with significant degree change
# deg.sig = which(!is.na(traj.type.deg.sig))
# deg.d.max.sig = deg.min.age.sig = array(NA,dim=c(nroi))
# deg.d.max.sig[deg.sig] = deg.d.max[deg.sig]
# deg.min.age.sig[deg.sig] = deg.min.age[deg.sig]
# write(deg.d.max.sig, file = paste(dir.path,'age_var/surf_txt/deg_d_max_sig.txt',sep=''), ncolumns = 1)
# write(deg.min.age.sig, file = paste(dir.path,'age_var/surf_txt/deg_min_age_sig.txt',sep=''), ncolumns = 1)

# plot trajectories, with averages, for each of four classes
col.traj.type = c('powderblue','lightcoral','steelblue','red3')

# all
spl.all = matrix(nrow = nroi, ncol = length(age.pred))
pdf(paste(dir.path,'age_var/traj/traj_type_deg.pdf',sep=''),width=5.2,height=5)
par(mfrow=c(2,2),oma=c(4,4,2,2),mar=c(2,2,1,1)+0.1,cex.lab=2,cex.axis=1.5,cex.main=2,font.main=1,bg='white')
for (i in 1:4) {
  print(i)
  npred = length(age.pred)
  plot(0, type='n', xaxt='n', xlim=c(min(age)-0.5,max(age)+0.5), ylim=c(-10,250))
  axis(1,at=seq(15,25,by=2.5),labels=c(15,'',20,'',25))
  axis(2,at=seq(0,250,by=50),labels=c(0,'',100,'',200,''))
  for (m in which(traj.type.deg==i)) {
    x = age.bin; g = gam(deg.bin[m,] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML');
    pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=F); spl.all[m,] = pred
    lines(age.pred,pred,col=col.traj.type[i],lwd=1)
  }
  if (length(which(traj.type.deg==i)) > 1) lines(age.pred,colMeans(spl.all[which(traj.type.deg==i),]), type = 'l', lty = 1, lwd = 3, col = 'grey10') # if cluster contains more than one traj, plot mean
}
mtext('median age (years)',outer=T,side=1,line=1.5,cex=2)
mtext('degree',outer=T,side=2,line=1.5,cex=2)
mtext('decrease',outer=T,side=3,line=-0.2,cex=1.75,adj=0.18); mtext('increase',outer=T,side=3,line=-0.2,cex=1.75,adj=0.88)
mtext('nonlinear',outer=T,side=4,line=0.2,cex=1.75,adj=0.18); mtext('linear',outer=T,side=4,line=0.2,cex=1.75,adj=0.85)
dev.off()

# significant
spl.all = matrix(nrow = nroi, ncol = length(age.pred))
pdf(paste(dir.path,'age_var/traj/traj_type_deg_sig.pdf',sep=''),width=5.2,height=5)
par(mfrow=c(2,2),oma=c(4,4,2,2),mar=c(2,2,1,1)+0.1,cex.lab=2,cex.axis=1.5,cex.main=2,font.main=1,bg='white')
for (i in 1:4) {
  print(i)
  npred = length(age.pred)
  plot(0, type='n', xaxt='n', xlim=c(min(age)-0.5,max(age)+0.5), ylim=c(-10,250))
  axis(1,at=seq(15,25,by=2.5),labels=c(15,'',20,'',25))
  axis(2,at=seq(0,250,by=50),labels=c(0,'',100,'',200,''))
  for (m in which(traj.type.deg.sig==i)) {
    x = age.bin; g = gam(deg.bin[m,] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML');
    pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=F); spl.all[m,] = pred
    lines(age.pred,pred,col=col.traj.type[i],lwd=1)
  }
  if (length(which(traj.type.deg.sig==i)) > 1) lines(age.pred,colMeans(spl.all[which(traj.type.deg.sig==i),]), type = 'l', lty = 1, lwd = 3, col = 'grey10') # if cluster contains more than one traj, plot mean
}
mtext('median age (years)',outer=T,side=1,line=1.5,cex=2)
mtext('degree',outer=T,side=2,line=1.5,cex=2)
mtext('decrease',outer=T,side=3,line=-0.2,cex=1.75,adj=0.18); mtext('increase',outer=T,side=3,line=-0.2,cex=1.75,adj=0.88)
mtext('nonlinear',outer=T,side=4,line=0.2,cex=1.75,adj=0.18); mtext('linear',outer=T,side=4,line=0.2,cex=1.75,adj=0.85)
dev.off()

### distance
# correct model p-values for multiple comparisons
dist.l.p.fdr = p.adjust(dist.l.p, method = 'fdr') # linear
dist.g.p.fdr = p.adjust(dist.g.p, method = 'fdr') # spline

dist.lin = which(dist.g.aic > dist.l.aic)                       # linear trajectories
length(intersect(dist.lin,which(dist.l.p.fdr<0.05)))            # number of significant linear trajectories
sum(dist.d.max[dist.lin]<0)                                     # number of linear decreases
sum(dist.d.max[dist.lin]>0)                                     # number of linear increases
sum(dist.d.max[intersect(dist.lin,which(dist.l.p.fdr<0.05))]<0) # number of significant linear decreases
sum(dist.d.max[intersect(dist.lin,which(dist.l.p.fdr<0.05))]>0) # number of significant linear increases

dist.nonlin = which(dist.g.aic < dist.l.aic)                       # non-linear trajectories
length(intersect(dist.nonlin,which(dist.g.p.fdr<0.05)))            # number of significant nonlinear trajectories
sum(dist.d.max[dist.nonlin]<0)                                     # number of nonlinear decreases
sum(dist.d.max[dist.nonlin]>0)                                     # number of nonlinear increases
sum(dist.d.max[intersect(dist.nonlin,which(dist.g.p.fdr<0.05))]<0) # number of significant nonlinear decreases
sum(dist.d.max[intersect(dist.nonlin,which(dist.g.p.fdr<0.05))]>0) # number of significant nonlinear increases

# plot trajectories by four categories (colours) - lin up / down // nonlin up / down (all + sig)
traj.type.dist = traj.type.dist.sig = array(NA,dim=c(nroi))
# all
traj.type.dist[intersect(dist.lin,which(dist.d.max<0))] = 1    # lin decrease
traj.type.dist[intersect(dist.lin,which(dist.d.max>0))] = 2    # lin increase
traj.type.dist[intersect(dist.nonlin,which(dist.d.max<0))] = 3 # nonlin decrease
traj.type.dist[intersect(dist.nonlin,which(dist.d.max>0))] = 4 # nonlin increase
# significant (FDR-adjusted threshold alphaFDR = 0.05)
traj.type.dist.sig[intersect(which(dist.d.max<0),intersect(dist.lin,which(dist.l.p.fdr<0.05)))] = 1     # lin decrease
traj.type.dist.sig[intersect(which(dist.d.max>0),intersect(dist.lin,which(dist.l.p.fdr<0.05)))] = 2     # lin increase
traj.type.dist.sig[intersect(which(dist.d.max<0),intersect(dist.nonlin,which(dist.g.p.fdr<0.05)))] = 3  # nonlin decrease
traj.type.dist.sig[intersect(which(dist.d.max>0),intersect(dist.nonlin,which(dist.g.p.fdr<0.05)))] = 4  # nonlin increase

# # # export trajectory types for plotting
# write(traj.type.dist, file = paste(dir.path,'age_var/surf_txt/traj_type_dist.txt',sep=''), ncolumns = 1)
# write(traj.type.dist.sig, file = paste(dir.path,'age_var/surf_txt/traj_type_dist_sig.txt',sep=''), ncolumns = 1)

# # plot dist.d.max and dist.min.age for regions with significant distance change
# dist.sig = which(!is.na(traj.type.dist.sig))
# dist.d.max.sig = dist.min.age.sig = array(NA,dim=c(nroi))
# dist.d.max.sig[dist.sig] = dist.d.max[dist.sig]
# dist.min.age.sig[dist.sig] = dist.min.age[dist.sig]
# write(dist.d.max.sig, file = paste(dir.path,'age_var/surf_txt/dist_d_max_sig.txt',sep=''), ncolumns = 1)
# write(dist.min.age.sig, file = paste(dir.path,'age_var/surf_txt/dist_min_age_sig.txt',sep=''), ncolumns = 1)

# plot trajectories, with averages, for each of four classes
col.traj.type = c('powderblue','lightcoral','steelblue','red3')

# all
spl.all = matrix(nrow = nroi, ncol = length(age.pred))
pdf(paste(dir.path,'age_var/traj/traj_type_dist.pdf',sep=''),width=5.2,height=5)
par(mfrow=c(2,2),oma=c(4,4,2,2),mar=c(2,2,1,1)+0.1,cex.lab=2,cex.axis=1.5,cex.main=2,font.main=1,bg='white')
for (i in 1:4) {
  print(i)
  npred = length(age.pred)
  plot(0, type='n', xaxt='n', xlim=c(min(age)-0.5,max(age)+0.5), ylim=c(0,130))
  axis(1,at=seq(15,25,by=2.5),labels=c(15,'',20,'',25))
  axis(2,at=seq(0,80,by=20),labels=c(0,'',40,'',80))
  for (m in which(traj.type.dist==i)) {
    x = age.bin; g = gam(dist.nod.bin[m,] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML');
    pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=F); spl.all[m,] = pred
    lines(age.pred,pred,col=col.traj.type[i],lwd=1)
  }
  if (length(which(traj.type.dist==i)) > 1) lines(age.pred,colMeans(spl.all[which(traj.type.dist==i),]), type = 'l', lty = 1, lwd = 3, col = 'grey10') # if cluster contains more than one traj, plot mean
}
mtext('median age (years)',outer=T,side=1,line=1.5,cex=2)
mtext('distance (mm)',outer=T,side=2,line=1.5,cex=2)
mtext('decrease',outer=T,side=3,line=-0.2,cex=1.75,adj=0.18); mtext('increase',outer=T,side=3,line=-0.2,cex=1.75,adj=0.88)
mtext('nonlinear',outer=T,side=4,line=0.2,cex=1.75,adj=0.18); mtext('linear',outer=T,side=4,line=0.2,cex=1.75,adj=0.85)
dev.off()

# significant
col.traj.type = c('powderblue','lightcoral','steelblue','red3')
spl.all = matrix(nrow = nroi, ncol = length(age.pred))
pdf(paste(dir.path,'age_var/traj/traj_type_dist_sig.pdf',sep=''),width=5.2,height=5)
par(mfrow=c(2,2),oma=c(4,4,2,2),mar=c(2,2,1,1)+0.1,cex.lab=2,cex.axis=1.5,cex.main=2,font.main=1,bg='white')
for (i in 1:4) {
  print(i)
  npred = length(age.pred)
  plot(0, type='n', xaxt='n', xlim=c(min(age)-0.5,max(age)+0.5), ylim=c(0,130))
  axis(1,at=seq(15,25,by=2.5),labels=c(15,'',20,'',25))
  axis(2,at=seq(0,80,by=20),labels=c(0,'',40,'',80))
  for (m in which(traj.type.dist.sig==i)) {
    x = age.bin; g = gam(dist.nod.bin[m,] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML');
    pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=F); spl.all[m,] = pred
    lines(age.pred,pred,col=col.traj.type[i],lwd=1)
  }
  if (length(which(traj.type.dist.sig==i)) > 1) lines(age.pred,colMeans(spl.all[which(traj.type.dist.sig==i),]), type = 'l', lty = 1, lwd = 3, col = 'grey10') # if cluster contains more than one traj, plot mean
}
mtext('median age (years)',outer=T,side=1,line=1.5,cex=2)
mtext('distance (mm)',outer=T,side=2,line=1.5,cex=2)
mtext('decrease',outer=T,side=3,line=-0.2,cex=1.75,adj=0.18); mtext('increase',outer=T,side=3,line=-0.2,cex=1.75,adj=0.88)
mtext('nonlinear',outer=T,side=4,line=0.2,cex=1.75,adj=0.18); mtext('linear',outer=T,side=4,line=0.2,cex=1.75,adj=0.85)
dev.off()

save('str.l.p.fdr','str.g.p.fdr','str.lin','str.nonlin','traj.type.str.d.max','traj.type.str.d.max.sig',
     'deg.l.p.fdr','deg.g.p.fdr','deg.lin','deg.nonlin','traj.type.deg.d.max','traj.type.deg.d.max.sig',
     'dist.l.p.fdr','dist.g.p.fdr','dist.lin','dist.nonlin','traj.type.dist.d.max','traj.type.dist.d.max.sig',
     file=paste(dir.path,'age_var/sc.age.var.traj.type.RData',sep=''))

# relationships ("relns") between regional measures ----------------------------------
if (!dir.exists(paste(dir.path,'age_var/relns',sep=''))) dir.create(paste(dir.path,'age_var/relns',sep=''),showWarning=T, recursive=F)

### relationships between measures of change and age at minimum

# "predicted" values (for regression line plotting)
str.d.pred = seq(from=-0.6,to=0.45,by=0.05) # max d strength
deg.d.pred = seq(from=-250,to=180,by=10)    # max d degree
dist.d.pred = seq(from=-100,to=80,by=10)    # max d distance

# strength (str.d.max vs str.age.min - combined)
id.nolim = intersect(which(str.min.age!=min(age.pred)),which(str.min.age!=max(age.pred))) # id's of regions whose minimum does not occur at a limit of the age-range
pdf(paste(dir.path,'age_var/relns/str_d_max_vs_str_age_min.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
# all regions, including those whose minima occur at limits (grey markers, dashed grey regression line)
x = str.d.max; l = lm(str.min.age~x); 
spear = cor.test(str.d.max,str.min.age,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(x, str.min.age, xlab=expression(paste(Delta,k[max],sep='')), ylab=expression(paste('age(',k[min],') (years)',sep='')), pch = 19, col='grey75', xlim=c(min(str.d.pred),max(str.d.pred)), ylim=c(13,26))
pred = predict(l,newdata=data.frame('x'=str.d.pred),interval='confidence')
polygon(c(rev(str.d.pred), str.d.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey85',alpha=0.5), border = NA)
lines(str.d.pred,pred[,1],col='grey50',lwd=2,lty=2); 
# only regions whose minima do not occur at limits (black markers, solid black regression line)
x = str.d.max[id.nolim]; l = lm(str.min.age[id.nolim]~x); 
spear = cor.test(str.d.max[id.nolim],str.min.age[id.nolim],method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
points(x, str.min.age[id.nolim],pch = 19, col='black')
pred = predict(l,newdata=data.frame('x'=str.d.pred),interval='confidence')
polygon(c(rev(str.d.pred), str.d.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(str.d.pred,pred[,1],col='black',lwd=3); 
dev.off()

# degree (deg.d.max vs deg.age.min - combined)
id.nolim = intersect(which(deg.min.age!=min(age.pred)),which(deg.min.age!=max(age.pred)))
pdf(paste(dir.path,'age_var/relns/deg_d_max_vs_deg_age_min.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
# all regions, including those whose minima occur at limits (grey markers, dashed grey regression line)
x = deg.d.max; l = lm(deg.min.age~x); 
spear = cor.test(deg.d.max,deg.min.age,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(x, deg.min.age, xlab=expression(paste(Delta,k[max],sep='')), ylab=expression(paste('age(',k[min],') (years)',sep='')), pch = 19, col='grey75', xlim=c(min(deg.d.pred),max(deg.d.pred)), ylim=c(13,26)) 
pred = predict(l,newdata=data.frame('x'=deg.d.pred),interval='confidence')
polygon(c(rev(deg.d.pred), deg.d.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey85',alpha=0.5), border = NA)
lines(deg.d.pred,pred[,1],col='grey50',lwd=2,lty=2); 
# only regions whose minima do not occur at limits (black markers, solid black regression line)
x = deg.d.max[id.nolim]; l = lm(deg.min.age[id.nolim]~x); 
spear = cor.test(deg.d.max[id.nolim],deg.min.age[id.nolim],method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
points(x, deg.min.age[id.nolim],pch = 19, col='black')
pred = predict(l,newdata=data.frame('x'=deg.d.pred),interval='confidence')
polygon(c(rev(deg.d.pred), deg.d.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(deg.d.pred,pred[,1],col='black',lwd=3); 
dev.off()

# distance (dist.d.max vs dist.age.min - combined)
id.nolim = intersect(which(dist.min.age!=min(age.pred)),which(dist.min.age!=max(age.pred)))
pdf(paste(dir.path,'age_var/relns/ddist_abs_vs_dist_age_min_comb.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
# all regions, including those whose minima occur at limits (grey markers, dashed grey regression line)
x = dist.d.max; l = lm(dist.min.age~x); 
spear = cor.test(dist.d.max,dist.min.age,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(x, dist.min.age, xlab=expression(paste(Delta,d[max],sep='')), ylab=expression(paste('age(',d[min],') (years)',sep='')), pch = 19, col='grey75', xlim=c(min(dist.d.pred),max(dist.d.pred)), ylim=c(13,26))
pred = predict(l,newdata=data.frame('x'=dist.d.pred),interval='confidence')
polygon(c(rev(dist.d.pred), dist.d.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey85',alpha=0.5), border = NA)
lines(dist.d.pred,pred[,1],col='grey50',lwd=2,lty=2); 
# only regions whose minima do not occur at limits (black markers, solid black regression line)
x = dist.d.max[id.nolim]; l = lm(dist.min.age[id.nolim]~x); 
spear = cor.test(dist.d.max[id.nolim],dist.min.age[id.nolim],method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
points(x, dist.min.age[id.nolim],pch = 19, col='black')
pred = predict(l,newdata=data.frame('x'=dist.d.pred),interval='confidence')
polygon(c(rev(dist.d.pred), dist.d.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(dist.d.pred,pred[,1],col='black',lwd=3); 
dev.off()

# deg.d.max vs dist.d.max
pdf(paste(dir.path,'age_var/relns/deg_d_max_vs_dist_d_max.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
x = deg.d.max; l = lm(dist.d.max~x);
spear = cor.test(deg.d.max,dist.d.max,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(x, dist.d.max, xlab=expression(paste(Delta,k[max],sep='')), ylab=expression(paste(Delta,d[max],' (mm)',sep='')), main = rp.main.sp(sp.rho,sp.p,2), pch = 19, col='black', xlim=c(min(deg.d.pred),max(deg.d.pred)), ylim=c(-100,80))
pred = predict(l,newdata=data.frame('x'=deg.d.pred),interval='confidence')
polygon(c(rev(deg.d.pred), deg.d.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(deg.d.pred,pred[,1],col='black',lwd=3); 
dev.off()

### relationships of (changes in) structural network to changes in morphology

# loop over regions, estimate slopes and interecepts for CT and MT 
# (main MT values correspond to 70% cortical depth, all depths used below)
ct.sl = ct.p = vector(length = nroi)
mt.sl = mt.p = vector(length = nroi)
for (n in 1:nroi) {
  # CT
  l = lm(ct[,n] ~ age + male)
  ct.sl[n] = l$coefficients[2]
  ct.p[n] = summary(l)$coefficients[2,4]
  # MT
  l = lm(mt[,n] ~ age + male)
  mt.sl[n] = l$coefficients[2]
  mt.p[n] = summary(l)$coefficients[2,4]
}

ct.p.fdr = p.adjust(ct.p, method = "fdr")
mt.p.fdr = p.adjust(mt.p, method = "fdr")

ct.sl.sig = mt.sl.sig = array(NA,dim=c(nroi,1))
ct.sl.sig[ct.p.fdr<0.05] = ct.sl[ct.p.fdr<0.05]
mt.sl.sig[mt.p.fdr<0.05] = mt.sl[mt.p.fdr<0.05]

# # optionally write out files for plotting
# if (!dir.exists(paste(dir.path,'age_var/relns/surf_txt',sep=''))) dir.create(paste(dir.path,'age_var/relns/surf_txt',sep=''),showWarning=T, recursive=F)
# write(ct.sl.sig, file = paste(dir.path,'age_var/relns/surf_txt/ct_sl_sig.txt',sep=''), ncolumns = 1)
# write(mt.sl.sig, file = paste(dir.path,'age_var/relns/surf_txt/mt_sl_sig.txt',sep=''), ncolumns = 1)

# "predicted" values (for regression line plotting)
str.pred = seq(from=10,to=125,by=5)           # strength (age-invariant)
deg.pred = seq(from=0,to=330,by=10)           # degree (age-invariant)
dct.pred = seq(from=-0.04,to=0.015,by=0.001)  # delta CT
dmt.pred = seq(from=-0.002,to=0.01,by=0.002)  # delta MT

# d.ct vs str.d.max
pdf(paste(dir.path,'age_var/relns/dct_vs_str_d_max.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
x = ct.sl; l = lm(str.d.max~x);
spear = cor.test(ct.sl,str.d.max,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(ct.sl,str.d.max,xaxt='n',yaxt='n',xlab=expression(paste(Delta,'CT (mm/year)',sep='')),ylab=expression(paste(Delta,s[max])), main = rp.main.sp(sp.rho,sp.p,2), pch = 19, col='black', xlim=c(min(dct.pred),max(dct.pred)), ylim=c(-0.6,0.4));
axis(1,at=seq(-0.04,0.01,by=0.01),labels=c(-0.04,'',-0.02,'',0,''))
axis(2,at=seq(-0.6,0.3,by=0.15),labels=c(-0.6,'',-0.3,'',0,'','0.3'))
pred = predict(l,newdata=data.frame('x'=dct.pred),interval='confidence')
polygon(c(rev(dct.pred), dct.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(dct.pred,pred[,1],col='black',lwd=3);
dev.off()

# d.mt vs str.d.max
pdf(paste(dir.path,'age_var/relns/dmt_vs_str_d_max.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
x = mt.sl; l = lm(str.d.max~x);
spear = cor.test(mt.sl,str.d.max,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(mt.sl,str.d.max,yaxt='n',xlab=expression(paste(Delta,'MT (PU/year)',sep='')),ylab=expression(paste(Delta,s[max])), main = rp.main.sp(sp.rho,sp.p,2), pch = 19, col='black', xlim=c(min(dmt.pred),max(dmt.pred)), ylim=c(-0.6,0.4));
axis(2,at=seq(-0.6,0.3,by=0.15),labels=c(-0.6,'',-0.3,'',0,'','0.3'))
pred = predict(l,newdata=data.frame('x'=dmt.pred),interval='confidence')
polygon(c(rev(dmt.pred), dmt.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(dmt.pred,pred[,1],col='black',lwd=3);
dev.off()

# d.ct vs deg.d.max
pdf(paste(dir.path,'age_var/relns/dct_vs_deg_d_max.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
x = ct.sl; l = lm(deg.d.max~x);
spear = cor.test(ct.sl,deg.d.max,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(ct.sl,deg.d.max,xaxt='n',yaxt='n',xlab=expression(paste(Delta,'CT (mm/year)',sep='')),ylab=expression(paste(Delta,k[max])), main = rp.main.sp(sp.rho,sp.p,2), pch = 19, col='black', xlim=c(min(dct.pred),max(dct.pred)), ylim=c(-300,220));
axis(1,at=seq(-0.04,0.01,by=0.01),labels=c(-0.04,'',-0.02,'',0,''))
axis(2,at=c(-300,0)); axis(2,at=c(-150,150))
pred = predict(l,newdata=data.frame('x'=dct.pred),interval='confidence')
polygon(c(rev(dct.pred), dct.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(dct.pred,pred[,1],col='black',lwd=3);
dev.off()

# d.ct vs deg.d.max - w/ and w/out outliers (regions with d.ct>0)
# excluding outliers
ct.decr.ind = which(ct.sl<0) # outlier indices
pdf(paste(dir.path,'age_var/relns/dct_vs_deg_d_max_outliers.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
# w/ outliers
x = ct.sl; l = lm(deg.d.max~x);
spear = cor.test(ct.sl,deg.d.max,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(ct.sl,deg.d.max,xaxt='n',yaxt='n',xlab=expression(paste(Delta,'CT (mm/year)',sep='')),ylab=expression(paste(Delta,k[max])), main = rp.main.sp(sp.rho,sp.p,2), pch = 19, col='grey75', xlim=c(min(dct.pred),max(dct.pred)), ylim=c(-300,220));
axis(1,at=seq(-0.04,0.01,by=0.01),labels=c(-0.04,'',-0.02,'',0,''))
axis(2,at=c(-300,0)); axis(2,at=c(-150,150))
pred = predict(l,newdata=data.frame('x'=dct.pred),interval='confidence')
polygon(c(rev(dct.pred), dct.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey85',alpha=0.5), border = NA)
lines(dct.pred,pred[,1],col='grey50',lwd=2,lty=2);
# w/out outliers
x = ct.sl[ct.decr.ind]; l = lm(deg.d.max[ct.decr.ind]~x);
spear = cor.test(ct.sl[ct.decr.ind],deg.d.max[ct.decr.ind],method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
points(ct.sl[ct.decr.ind],deg.d.max[ct.decr.ind],pch = 19, col='black');
pred = predict(l,newdata=data.frame('x'=dct.pred),interval='confidence')
polygon(c(rev(dct.pred), dct.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(dct.pred,pred[,1],col='black',lwd=3);
dev.off()

# d.mt vs deg.d.max
pdf(paste(dir.path,'age_var/relns/dmt_vs_deg_d_max.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
x = mt.sl; l = lm(deg.d.max~x);
spear = cor.test(mt.sl,deg.d.max,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(mt.sl,deg.d.max,yaxt='n',xaxt='n',xlab=expression(paste(Delta,'MT (PU/year)',sep='')),ylab=expression(paste(Delta,k[max])), main = rp.main.sp(sp.rho,sp.p,2), pch = 19, col='black', xlim=c(min(dmt.pred),max(dmt.pred)), ylim=c(-300,220));
axis(1,at=seq(-0.002,0.01,by=0.002),labels=c('',0,'',0.004,'',0.008,''))
axis(2,at=c(-300,0)); axis(2,at=c(-150,150))
pred = predict(l,newdata=data.frame('x'=dmt.pred),interval='confidence')
polygon(c(rev(dmt.pred), dmt.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(dmt.pred,pred[,1],col='black',lwd=3);
dev.off()

# d.ct vs str.age.min
id.nolim = intersect(which(str.min.age!=min(age.pred)),which(str.min.age!=max(age.pred))) # regions whose minima do not occur at limits of the age range
# all regions including those whose minimum occurs at one of the limits (grey markers, dashed grey regression line)
pdf(paste(dir.path,'age_var/relns/dct_vs_str_age_min_comb.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
x = ct.sl; l = lm(str.min.age~x);
spear = cor.test(ct.sl,str.min.age,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(ct.sl,str.min.age,xaxt='n',xlab=expression(paste(Delta,'CT (mm/year)',sep='')),ylab=expression(paste('age(',s[min],') (years)',sep='')), pch = 19, col='grey75', xlim=c(min(dct.pred),max(dct.pred)), ylim=c(min(age)-1,max(age)+1));
axis(1,at=seq(-0.04,0.01,b=0.01),labels=c(-0.04,'',-0.02,'',0,''))
pred = predict(l,newdata=data.frame('x'=dct.pred),interval='confidence')
polygon(c(rev(dct.pred), dct.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey85',alpha=0.5), border = NA)
lines(dct.pred,pred[,1],col='grey50',lwd=2,lty=2);
# only regions including whose minimum does not occur at one of the limits (black markers, solid black regression line)
x = ct.sl[id.nolim]; l = lm(str.min.age[id.nolim]~x);
spear = cor.test(ct.sl[id.nolim],str.min.age[id.nolim],method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
points(ct.sl[id.nolim],str.min.age[id.nolim],pch = 19, col='black')
axis(1,at=seq(-0.04,0.01,b=0.01),labels=c(-0.04,'',-0.02,'',0,''))
pred = predict(l,newdata=data.frame('x'=dct.pred),interval='confidence')
polygon(c(rev(dct.pred), dct.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(dct.pred,pred[,1],col='black',lwd=3);
dev.off()

# d.mt vs str.age.min
id.nolim = intersect(which(str.min.age!=min(age.pred)),which(str.min.age!=max(age.pred))) # regions whose minima do not occur at limits of the age range
# all regions including those whose minimum occurs at one of the limits (grey markers, dashed grey regression line)
pdf(paste(dir.path,'age_var/relns/dmt_vs_str_age_min_comb.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
x = mt.sl; l = lm(str.min.age~x);
spear = cor.test(mt.sl,str.min.age,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(mt.sl,str.min.age,xaxt='n',xlab=expression(paste(Delta,'MT (PU/year)',sep='')),ylab=expression(paste('age(',k[min],') (years)',sep='')), pch = 19, col='grey75', xlim=c(min(dmt.pred),max(dmt.pred)), ylim=c(min(age)-1,max(age)+1));
axis(1,at=seq(-0.002,0.01,by=0.002),labels=c('',0,'',0.004,'',0.008,''))
pred = predict(l,newdata=data.frame('x'=dmt.pred),interval='confidence')
polygon(c(rev(dmt.pred), dmt.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey85',alpha=0.5), border = NA)
lines(dmt.pred,pred[,1],col='grey50',lwd=2,lty=2);
# only regions including whose minimum does not occur at one of the limits (black markers, solid black regression line)
x = mt.sl[id.nolim]; l = lm(str.min.age[id.nolim]~x); 
spear = cor.test(mt.sl[id.nolim],str.min.age[id.nolim],method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
points(mt.sl[id.nolim],str.min.age[id.nolim],pch = 19, col='black')
pred = predict(l,newdata=data.frame('x'=dmt.pred),interval='confidence')
polygon(c(rev(dmt.pred), dmt.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(dmt.pred,pred[,1],col='black',lwd=3);
dev.off()

# d.ct VS deg.age.min
id.nolim = intersect(which(deg.min.age!=min(age.pred)),which(deg.min.age!=max(age.pred))) # regions whose minima do not occur at limits of the age range
# all regions including those whose minimum occurs at one of the limits (grey markers, dashed grey regression line)
pdf(paste(dir.path,'age_var/relns/dct_vs_deg_age_min_comb.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
x = ct.sl; l = lm(deg.min.age~x);
spear = cor.test(ct.sl,deg.min.age,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(ct.sl,deg.min.age,xaxt='n',xlab=expression(paste(Delta,'CT (mm/year)',sep='')),ylab=expression(paste('age(',k[min],') (years)',sep='')), pch = 19, col='grey75', xlim=c(-0.04,0.015), ylim=c(min(age)-1,max(age)+1));
axis(1,at=seq(-0.04,0.01,b=0.01),labels=c(-0.04,'',-0.02,'',0,''))
pred = predict(l,newdata=data.frame('x'=dct.pred),interval='confidence')
polygon(c(rev(dct.pred), dct.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey85',alpha=0.5), border = NA)
lines(dct.pred,pred[,1],col='grey50',lwd=2,lty=2);
# only regions including whose minimum does not occur at one of the limits (black markers, solid black regression line)
x = ct.sl[id.nolim]; l = lm(deg.min.age[id.nolim]~x);
spear = cor.test(ct.sl[id.nolim],deg.min.age[id.nolim],method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
points(ct.sl[id.nolim],deg.min.age[id.nolim],pch = 19, col='black')
axis(1,at=seq(-0.04,0.01,b=0.01),labels=c(-0.04,'',-0.02,'',0,''))
pred = predict(l,newdata=data.frame('x'=dct.pred),interval='confidence')
polygon(c(rev(dct.pred), dct.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(dct.pred,pred[,1],col='black',lwd=3);
dev.off()

# d.mt vs deg.age.min
id.nolim = intersect(which(deg.min.age!=min(age.pred)),which(deg.min.age!=max(age.pred))) # regions whose minima do not occur at limits of the age range
# all regions including those whose minimum occurs at one of the limits (grey markers, dashed grey regression line)
pdf(paste(dir.path,'age_var/relns/dmt_vs_deg_age_min_comb.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
x = mt.sl; l = lm(deg.min.age~x);
spear = cor.test(mt.sl,deg.min.age,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(mt.sl,deg.min.age,xaxt='n',xlab=expression(paste(Delta,'MT (PU/year)',sep='')),ylab=expression(paste('age(',k[min],') (years)',sep='')), pch = 19, col='grey75', xlim=c(min(dmt.pred),max(dmt.pred)), ylim=c(min(age)-1,max(age)+1));
axis(1,at=seq(-0.002,0.01,by=0.002),labels=c('',0,'',0.004,'',0.008,''))
pred = predict(l,newdata=data.frame('x'=dmt.pred),interval='confidence')
polygon(c(rev(dmt.pred), dmt.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey85',alpha=0.5), border = NA)
lines(dmt.pred,pred[,1],col='grey50',lwd=2,lty=2);
# only regions including whose minimum does not occur at one of the limits (black markers, solid black regression line)
x = mt.sl[id.nolim]; l = lm(deg.min.age[id.nolim]~x); 
spear = cor.test(mt.sl[id.nolim],deg.min.age[id.nolim],method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
points(mt.sl[id.nolim],deg.min.age[id.nolim],pch = 19, col='black')
pred = predict(l,newdata=data.frame('x'=dmt.pred),interval='confidence')
polygon(c(rev(dmt.pred), dmt.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(dmt.pred,pred[,1],col='black',lwd=3);
dev.off()

### relationship of change in myelination to changes in network properties as a function of cortical depth

# mt.depth contains the magnetization transfer (a measure of myelination) as a function of cortical depth
# it increases from the pial surface to the grey/white matter boundary and into white matter
ndepth = dim(mt.depth)[1]                                 # number of depths
depth.lab = c('pial','10%','20%','30%','40%','50%','60%','70%','80%','90%','GM/WM','0.4mm','0.8mm') # depth labels

# to convince yourself of its ordering, plot the values across regions and depths for an example subject; eg:
image.plot(mt.depth[,1,],xaxt='n',yaxt='n',ylab='region ID',xlab='cortical depth')
title('MT as a function of depth for example subject')
axis(1,at=seq(from=0,to=1,length.out=ndepth),labels=depth.lab,las=3)

# change in myelination and its relationship to str.d.max and deg.d.max *as a function of cortical depth*
dmt.depth.pred = seq(from=-0.007,to=0.01,by=0.002)            # predicted values (for regression curve)
dmt.depth = array(NA,dim=c(nroi,ndepth))                      # initialise delta MT as a function of depth
dmt.vs.dstr.depth = dmt.vs.ddeg.depth = vector(length=ndepth) # initialise relationships between d.mt and str, deg
for (i in 1:ndepth) {
  print(i)
  # loop over rois, estimate betas for MT
  for (n in 1:nroi) {
    l = lm(mt.depth[i,,n] ~ age + male)
    dmt.depth[n,i] = l$coefficients[2]
  }
  # d.mt VS str.d.max
  spear = cor.test(str.d.max,dmt.depth[,i],method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
  dmt.vs.dstr.depth[i] = sp.rho #rsq
  # d.mt VS deg.d.max
  spear = cor.test(deg.d.max,dmt.depth[,i],method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
  dmt.vs.ddeg.depth[i] = sp.rho #rsq
}

# dmt as a function of depth
library(vioplot)
pdf(paste(dir.path,'age_var/relns/dmt_all.pdf',sep=''),width=4,height=8)
par(mar=c(6, 6.5, 1, 2.5) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white', las=1)
plot(0, type='n', yaxt='n', xlim=c(-0.008,0.01), ylim=c(0.75,13.25),xlab='',ylab='')
abline(v=c(0),lty=2,col='grey')
for (i in 1:ndepth) vioplot(dmt.depth[,i],at=ndepth-(i-1),col='grey90',add=T,wex=0.75,horizontal=T) # violin plot
axis(side=1,at=seq(0,0.1,by=0.025),labels=c(0,'',0.05,'',0.1))
axis(side=2,at=c(1:ndepth),labels=rev(depth.lab))
abline(h=c(3,13),lty=2,col='grey')
par(las=0); mtext(side=1,bquote(Delta*'MT'~'(PU/year)'),line=4,cex=2)
dev.off()

# dmt VS d.str.max as a function of depth
pdf(paste(dir.path,'age_var/relns/dmt_vs_str_d_max_depth.pdf',sep=''),width=4,height=8)
par(mar=c(6, 6.5, 1, 2.5) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white', las=1)
plot(rev(dmt.vs.dstr.depth),c(1:ndepth),yaxt='n',pch=19,cex=0.75,xlab='',ylab='',xlim=c(-0.4,0),ylim=c(0.75,13.25))
lines(rev(dmt.vs.dstr.depth),c(1:ndepth))
axis(side=2,at=c(1:ndepth),labels=rev(depth.lab))
abline(h=c(3,13),lty=2,col='grey')
par(las=0); mtext(side=1,bquote(Delta*'MT VS'~Delta*s[max]~'('*rho*')'),line=4,cex=2)
dev.off()

# dmt VS d.deg.max as a function of depth
pdf(paste(dir.path,'age_var/relns/dmt_vs_deg_d_max_depth.pdf',sep=''),width=4,height=8)
par(mar=c(6, 6.5, 1, 2.5) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white', las=1)
plot(rev(dmt.vs.ddeg.depth),c(1:ndepth),yaxt='n',pch=19,cex=0.75,xlab='',ylab='',xlim=c(-0.4,0),ylim=c(0.75,13.25))
lines(rev(dmt.vs.ddeg.depth),c(1:ndepth))
axis(side=2,at=c(1:ndepth),labels=rev(depth.lab))
abline(h=c(3,13),lty=2,col='grey')
par(las=0); mtext(side=1,bquote(Delta*'MT VS'~Delta*k[max]~'('*rho*')'),line=4,cex=2)
dev.off()

### linear relationships
# d.ct VS d.str (linear slope) 
pdf(paste(dir.path,'age_var/relns/dct_vs_dstr.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
x = ct.sl; l = lm(str.sl~x); 
spear = cor.test(ct.sl,deg.sl,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(ct.sl,str.sl,xaxt='n',xlab=expression(paste(Delta,' CT (mm/year)',sep='')),ylab=expression(paste(Delta,s[lin],' (year'^'-1',')')), main = rp.main.sp(sp.rho,sp.p,2), pch = 19, col='black', xlim=c(-0.04,0.015), ylim=c(-0.04,0.03));
axis(1,at=seq(-0.04,0.01,b=0.01),labels=c(-0.04,'',-0.02,'',0,''))
pred = predict(l,newdata=data.frame('x'=dct.pred),interval='confidence')
polygon(c(rev(dct.pred), dct.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(dct.pred,pred[,1],col='black',lwd=3); 
dev.off()

# d.mt VS d.str (linear slope) 
pdf(paste(dir.path,'age_var/relns/dmt_vs_dstr.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
x = mt.sl; l = lm(str.sl~x); 
spear = cor.test(ct.sl,deg.sl,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(mt.sl,str.sl,xlab=expression(paste(Delta,'MT (PU/year)',sep='')),ylab=expression(paste(Delta,s[lin],' (year'^'-1',')')), main = rp.main.sp(sp.rho,sp.p,2), pch = 19, col='black', xlim=c(min(dmt.pred),max(dmt.pred)), ylim=c(-0.04,0.03));
pred = predict(l,newdata=data.frame('x'=dmt.pred),interval='confidence')
polygon(c(rev(dmt.pred), dmt.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(dmt.pred,pred[,1],col='black',lwd=3); 
dev.off()

# d.ct VS d.deg (linear slope)
pdf(paste(dir.path,'age_var/relns/dct_vs_ddeg.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
x = ct.sl; l = lm(deg.sl~x);
spear = cor.test(ct.sl,deg.sl,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(ct.sl,deg.sl,xaxt='n',xlab=expression(paste(Delta,'CT (mm/year)',sep='')),ylab=expression(paste(Delta,k[lin],' (year'^'-1',')')), main = rp.main.sp(sp.rho,sp.p,2), pch = 19, col='black', xlim=c(-0.04,0.015), ylim=c(-25,15));
axis(1,at=seq(-0.04,0.01,b=0.01),labels=c(-0.04,'',-0.02,'',0,''))
pred = predict(l,newdata=data.frame('x'=dct.pred),interval='confidence')
polygon(c(rev(dct.pred), dct.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(dct.pred,pred[,1],col='black',lwd=3); 
dev.off()

# d.mt VS d.deg (linear slope)
pdf(paste(dir.path,'age_var/relns/dmt_vs_ddeg.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
x = mt.sl; l = lm(deg.sl~x);
spear = cor.test(mt.sl,deg.sl,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(mt.sl,deg.sl,xaxt='n',xlab=expression(paste(Delta,'MT (PU/year)',sep='')),ylab=expression(paste(Delta,k[lin],' (year'^'-1',')')), main = rp.main.sp(sp.rho,sp.p,2), pch = 19, col='black', xlim=c(min(dmt.pred),max(dmt.pred)), ylim=c(-25,10));
axis(1,at=seq(-0.002,0.01,by=0.002),labels=c('',0,'',0.004,'',0.008,''))
pred = predict(l,newdata=data.frame('x'=dmt.pred),interval='confidence')
polygon(c(rev(dmt.pred), dmt.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(dmt.pred,pred[,1],col='black',lwd=3);
dev.off()

### compare properties of age-invariant network to age-varying network changes
str.pred = seq(from=0,to=0.45,by=0.01)  # strength
deg.pred = seq(from=0,to=330,by=10)     # degree
wdeg.pred = seq(0,120,by=2);            # weighted degree

# strength

# str (age-inv) VS str.d.max
pdf(paste(dir.path,'age_var/relns/str_vs_dstr_max.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
x = str; l = lm(str.d.max~x); f = summary(l)$fstatistic; p = pf(f[1],f[2],f[3], lower.tail = FALSE); rsq = summary(l)$r.squared
plot(str,str.d.max,xlab='s',ylab=expression(paste(Delta,s[max],sep='')), main = rp.main(rsq,p,2), pch = 19, col='black', xlim=c(min(str.pred),max(str.pred)), ylim=c(-0.6,0.45));
pred = predict(l,newdata=data.frame('x'=str.pred),interval='confidence')
polygon(c(rev(str.pred), str.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(str.pred,pred[,1],col='black',lwd=3); 
dev.off()

# str (age-inv) vs str.min.age
id.nolim = intersect(which(str.min.age!=min(age.pred)),which(str.min.age!=max(age.pred))) # regions whose minimum does not occur at limits of the age range
pdf(paste(dir.path,'age_var/relns/str_vs_str_age_min_comb.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
# all regions including those whose minimum occurs at one of the limits (grey markers, dashed grey regression line)
x = str; l = lm(str.min.age~x); 
spear = cor.test(str,str.min.age,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(x, str.min.age, xlab=expression(paste(s,sep='')), ylab=expression(paste('age(',s[min],') (years)',sep='')), pch = 19, col='grey75', xlim=c(min(str.pred),max(str.pred)), ylim=c(13,26)) 
pred = predict(l,newdata=data.frame('x'=str.pred),interval='confidence')
polygon(c(rev(str.pred), str.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey85',alpha=0.5), border = NA)
lines(str.pred,pred[,1],col='grey50',lwd=2,lty=2); 
# only regions including whose minimum does not occur at one of the limits (black markers, solid black regression line)
x = str[id.nolim]; l = lm(str.min.age[id.nolim]~x); 
spear = cor.test(str[id.nolim],str.min.age[id.nolim],method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
points(x, str.min.age[id.nolim],pch = 19, col='black')
pred = predict(l,newdata=data.frame('x'=str.pred),interval='confidence')
polygon(c(rev(str.pred), str.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(str.pred,pred[,1],col='black',lwd=3); 
dev.off()

# degree

# deg (age-inv) VS deg.d.max
pdf(paste(dir.path,'age_var/relns/deg_vs_ddeg_max.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2.25, font.main = 1, bg='white')
x = deg; l = lm(deg.d.max~x);
spear = cor.test(deg,deg.d.max,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(deg,deg.d.max,xlab='k',ylab=expression(paste(Delta,k[max],sep='')), main = rp.main.sp(sp.rho,sp.p,2), pch = 19, col='black', xlim=c(min(deg.pred),max(deg.pred)), ylim=c(-300,200));
pred = predict(l,newdata=data.frame('x'=deg.pred),interval='confidence')
polygon(c(rev(deg.pred), deg.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(deg.pred,pred[,1],col='black',lwd=3); 
dev.off()

# deg (age-inv) vs deg.min.age
id.nolim = intersect(which(deg.min.age!=min(age.pred)),which(deg.min.age!=max(age.pred))) # regions whose minimum does not occur at limits of the age range
pdf(paste(dir.path,'age_var/relns/deg_vs_deg_age_min_comb.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
# all regions including those whose minimum occurs at one of the limits (grey markers, dashed grey regression line)
x = deg; l = lm(deg.min.age~x); 
spear = cor.test(deg,deg.min.age,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(x, deg.min.age, xlab=expression(paste(k,sep='')), ylab=expression(paste('age(',k[min],') (years)',sep='')), pch = 19, col='grey75', xlim=c(min(deg.pred),max(deg.pred)), ylim=c(13,26)) 
pred = predict(l,newdata=data.frame('x'=deg.pred),interval='confidence')
polygon(c(rev(deg.pred), deg.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey85',alpha=0.5), border = NA)
lines(deg.pred,pred[,1],col='grey50',lwd=2,lty=2); 
# only regions including whose minimum does not occur at one of the limits (black markers, solid black regression line)
x = deg[id.nolim]; l = lm(deg.min.age[id.nolim]~x); 
spear = cor.test(deg[id.nolim],deg.min.age[id.nolim],method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
points(x, deg.min.age[id.nolim],pch = 19, col='black')
pred = predict(l,newdata=data.frame('x'=deg.pred),interval='confidence')
polygon(c(rev(deg.pred), deg.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(deg.pred,pred[,1],col='black',lwd=3); 
dev.off()

# weighted degree

# wdeg (age-inv) VS deg.d.max
pdf(paste(dir.path,'age_var/relns/wdeg_vs_ddeg_max.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
x = wdeg; l = lm(deg.d.max~x);
spear = cor.test(wdeg,deg.d.max,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(wdeg,deg.d.max,xlab=expression(paste(k[w],sep='')),ylab=expression(paste(Delta,k[max],sep='')), main = rp.main.sp(sp.rho,sp.p,2), pch = 19, col='black', xlim=c(min(wdeg.pred),max(wdeg.pred)), ylim=c(-300,200));
pred = predict(l,newdata=data.frame('x'=wdeg.pred),interval='confidence')
polygon(c(rev(wdeg.pred), wdeg.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(wdeg.pred,pred[,1],col='black',lwd=3); 
dev.off()

# wdeg (age-inv) vs deg.min.age
id.nolim = intersect(which(deg.min.age!=min(age.pred)),which(deg.min.age!=max(age.pred)))
pdf(paste(dir.path,'age_var/relns/wdeg_vs_deg_age_min_comb.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2.25, cex.axis = 1.75, cex.main = 2.25, font.main = 1, bg='white')
# all regions including those whose minimum occurs at one of the limits (grey markers, dashed grey regression line)
x = wdeg; l = lm(deg.min.age~x);
spear = cor.test(wdeg,deg.min.age,method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
plot(x, deg.min.age, xlab=expression(paste(k[w],sep='')), ylab=expression(paste('age(',k[min],') (years)',sep='')), pch = 19, col='grey75', xlim=c(min(wdeg.pred),max(wdeg.pred)), ylim=c(13,26)) 
pred = predict(l,newdata=data.frame('x'=wdeg.pred),interval='confidence')
polygon(c(rev(wdeg.pred), wdeg.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey85',alpha=0.5), border = NA)
lines(wdeg.pred,pred[,1],col='grey50',lwd=2,lty=2); 
# only regions including whose minimum does not occur at one of the limits (black markers, solid black regression line)
x = wdeg[id.nolim]; l = lm(deg.min.age[id.nolim]~x); 
spear = cor.test(wdeg[id.nolim],deg.min.age[id.nolim],method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
points(x, deg.min.age[id.nolim],pch = 19, col='black')
pred = predict(l,newdata=data.frame('x'=wdeg.pred),interval='confidence')
polygon(c(rev(wdeg.pred), wdeg.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
lines(wdeg.pred,pred[,1],col='black',lwd=3); 
dev.off()

# changes within and between communities of nodes -------------------------
if (!dir.exists(paste(dir.path,'age_var/communities/',sep=''))) dir.create(paste(dir.path,'age_var/communities/',sep=''),showWarning=T, recursive=F)

# modules of the age-invariant network (obtained with Louvain and versatility codes in Matlab)
nmod = max(mod.id)
col.mod <- c('red3','steelblue3','darkseagreen4','lightcoral','gold2','darkviolet','darkorange')

# cytoarchitectonic classes of the von Economo atlas (von Economo and Koskinas, 1925)
nve = max(ve.id) # number of classes
col.ve = c('darkmagenta','blue1','forestgreen','orange1','gold1','cyan','magenta') # class colours

# resting-state networks (Yeo, Krienen et al., J. Neurophysiol. 2011)
nyeo = max(yeo.id)
col.yeo <- c('darkmagenta','cadetblue','forestgreen','hotpink2','cornsilk1','orange','salmon')

# blue-red colobar used for plotting maximum change
x=matrix(rnorm(100),nrow=10)*100; xmin=0; xmax=100; x[x<xmin]=xmin; x[x>xmax]=xmax;
collist <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0","#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")
ColorRamp<-colorRampPalette(collist)(10000); ColorLevels<-seq(from=xmin, to=xmax, length=10000)
bluered.col <- ColorRamp[round(1+(min(x)-xmin)*10000/(xmax-xmin)) : round( (max(x)-xmin)*10000/(xmax-xmin) )]
rm(collist,ColorRamp,ColorLevels)

# cool colorbar used for plotting age at minimum
cool = c("#00FFFF", "#0DF2FF", "#19E6FF", "#26D9FF","#33CCFF", "#3FBFFF", "#4CB3FF", "#59A6FF", "#6699FF", "#738CFF", 
         "#7F7FFF", "#8C73FF", "#9966FF", "#A659FF", "#B24DFF", "#BF3FFF", "#CC33FF", "#D926FF", "#E619FF", "#F20DFF")

# relationships between communities
library(NMI)
NMI(data.frame(c(1:nroi),mod.id),data.frame(c(1:nroi),ve.id))   # 0.2894
NMI(data.frame(c(1:nroi),mod.id),data.frame(c(1:nroi),yeo.id))  # 0.3848
NMI(data.frame(c(1:nroi),ve.id),data.frame(c(1:nroi),yeo.id))   # 0.2875

### modules

# average strength (str) and partial density (pdens) by module
# ("partial" density here is relative to the total number of edges than can exist within a community, or between two communities)
mod.str = mod.pdens = array(NA,dim=c(nmod,nmod,nbin))
for (b in 1:nbin) {
  print(b)
  for (i in 1:nmod) {
    for (j in 1:i) {
      if (i == j) { # within-module (upper triangular only)
        mod.str[i,j,b] = mean(sc.bin[mod.id==i,mod.id==j,b][triup[mod.id==i,mod.id==j]])
        mod.pdens[i,j,b] = sum(!is.na(sc.bin.thr[mod.id==i,mod.id==j,b]))/(sum(mod.id==i)*(sum(mod.id==j)-1))
      } else { # between-module
        mod.str[i,j,b] = mean(sc.bin[mod.id==i,mod.id==j,b])
        mod.pdens[i,j,b] = sum(!is.na(sc.bin.thr[mod.id==i,mod.id==j,b]))/(sum(mod.id==i)*sum(mod.id==j))
      }
    }
  }
}

# linear models
mod.str.l.sl = mod.str.l.aic = mod.str.l.p = array(NA,dim=c(nmod,nmod))
mod.str.b.l.sl = mod.str.b.l.aic = array(NA,dim=c(nmod,nmod,nboot))
for (i in 1:nmod) {
  print(i)
  for (j in 1:i) {
    print(j)
    l = lm(mod.str[i,j,] ~ age.bin); 
    mod.str.l.sl[i,j] = l$coefficients[2]; 
    mod.str.l.aic[i,j] = AIC(l); 
    f = summary(l)$fstatistic; 
    mod.str.l.p[i,j] = pf(f[1],f[2],f[3], lower.tail = FALSE)
  }
}

# splines
mod.str.df = mod.str.g.aic = mod.str.g.p = mod.str.min.age = mod.str.d.max = array(NA,dim=c(nmod,nmod)) # strength
mod.pdens.df = mod.pdens.g.aic = mod.pdens.g.p = mod.pdens.min.age = mod.pdens.d.max = array(NA,dim=c(nmod,nmod)) # density
for (i in 1:nmod) {
  print(i)
  for (j in 1:i) {
    print(j)
    x = age.bin
    # strength
    g = gam(mod.str[i,j,] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML'); mod.str.df[i,j] = sum(g$edf); mod.str.g.p[i,j] = summary(g)$s.pv; mod.str.g.aic[i,j] = AIC(g)
    pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=F); mod.str.min.age[i,j] = age.pred[which(pred==min(pred))]
    mod.str.d.max[i,j] = (max(pred)-min(pred))/sign(age.pred[which(pred==max(pred))]-age.pred[which(pred==min(pred))])
    # density
    g = gam(mod.pdens[i,j,] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML'); mod.pdens.df[i,j] = sum(g$edf); mod.pdens.g.p[i,j] = summary(g)$s.pv; mod.pdens.g.aic[i,j] = AIC(g)
    pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=F); mod.pdens.min.age[i,j] = age.pred[which(pred==min(pred))]
    mod.pdens.d.max[i,j] = (max(pred)-min(pred))/sign(age.pred[which(pred==max(pred))]-age.pred[which(pred==min(pred))])
  }
}

# adjust p-values for multiple comparisons (FDR)
mod.pdens.g.p.vec = vector(length=nmod*(nmod+1)/2); c = 0
for (i in 1:nve) { for (j in 1:i) { c = c+1; mod.pdens.g.p.vec[c] = mod.pdens.g.p[i,j] } }
mod.pdens.g.p.fdr = p.adjust(mod.pdens.g.p.vec, method = 'fdr')

# mod.d.max matrix
pdf(paste(dir.path,'age_var/communities/mod_dens_d_max_mat.pdf',sep=''),width=8,height=7)
par(mar=c(3,3,6,8), cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
image.plot(100*mod.pdens.d.max[1:nmod,rev(1:nmod)],add=F,axes=F,legend.mar=8,col=bluered.col,zlim=(c(-max(abs(100*mod.pdens.d.max),na.rm=T),max(abs(100*mod.pdens.d.max),na.rm=T))))
for (j in 1:nmod) {
  points((j-1)/(nmod-1), 1+(1.2/nmod), pch = 19, cex = 3.9, col = col.mod[j],xpd=T)
  points(1+(1.2/nmod), (j-1)/(nmod-1), pch = 19, cex = 3.9, col = col.mod[rev(1:nmod)[j]],xpd=T)
}
c = 0
for (i in 1:nmod) {
  for (j in 1:i) {
    c=c+1
    if (mod.pdens.g.p.fdr[c] < 0.05) points((i-1)/(nmod-1), (rev(1:nmod)[j]-1)/(nmod-1), pch = 19, cex = 1, col = 'black')
    if (mod.pdens.g.p.fdr[c] < 0.01) points((i-1)/(nmod-1), (rev(1:nmod)[j]-1)/(nmod-1), pch = 19, cex = 2, col = 'black')
  }
}
dev.off()

# mod.age.min matrix
pdf(paste(dir.path,'age_var/communities/mod_dens_min_age_mat.pdf',sep=''),width=8,height=7)
par(mar=c(3,3,6,8), cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
image.plot(mod.pdens.min.age[1:nmod,rev(1:nmod)],add=F,axes=F,legend.mar=8,col=cool,zlim=c(min(age),max(age)))
for (j in 1:nmod) {
  points((j-1)/(nmod-1), 1+(1.2/nmod), pch = 19, cex = 3.9, col = col.mod[j],xpd=T)
  points(1+(1.2/nmod), (j-1)/(nmod-1), pch = 19, cex = 3.9, col = col.mod[rev(1:nmod)[j]],xpd=T)
}
c = 0
for (i in 1:nmod) {
  for (j in 1:i) {
    c=c+1
    if (mod.pdens.g.p.fdr[c] < 0.05) points((i-1)/(nmod-1), (rev(1:nmod)[j]-1)/(nmod-1), pch = 19, cex = 1, col = 'black')
    if (mod.pdens.g.p.fdr[c] < 0.01) points((i-1)/(nmod-1), (rev(1:nmod)[j]-1)/(nmod-1), pch = 19, cex = 2, col = 'black')
  }
}
dev.off()

### von economo classes

# average strength (str) and partial density (pdens) by von economo class
# ("partial" density here is relative to the total number of edges than can exist within a community, or between two communities)
ve.str = ve.pdens = array(NA,dim=c(nve,nve,nbin))
for (b in 1:nbin) {
  print(b)
  for (i in 1:nve) {
    for (j in 1:i) {
      if (i == j) { # within-class (upper triangular only)
        ve.str[i,j,b] <- mean(sc.bin[ve.id==i,ve.id==j,b][triup[ve.id==i,ve.id==j]])
        ve.pdens[i,j,b] = sum(!is.na(sc.bin.thr[ve.id==i,ve.id==j,b]))/(sum(ve.id==i)*(sum(ve.id==j)-1))
      } else {      # between-class
        ve.str[i,j,b] <- mean(sc.bin[ve.id==i,ve.id==j,b])
        ve.pdens[i,j,b] = sum(!is.na(sc.bin.thr[ve.id==i,ve.id==j,b]))/(sum(ve.id==i)*sum(ve.id==j))
      }
    }
  }
}

# linear models
ve.str.l.sl = ve.str.l.aic = ve.str.l.p = array(NA,dim=c(nve,nve))
for (i in 1:nve) {
  print(i)
  for (j in 1:i) {
    print(j)
    l = lm(ve.str[i,j,] ~ age.bin); 
    ve.str.l.sl[i,j] = l$coefficients[2]; 
    ve.str.l.aic[i,j] = AIC(l); 
    f = summary(l)$fstatistic; ve.str.l.p[i,j] = pf(f[1],f[2],f[3], lower.tail = FALSE)
  }
}

# splines
ve.str.df = ve.str.g.aic = ve.str.g.p = ve.str.min.age = ve.str.d.max = array(NA,dim=c(nve,nve)) # strength
ve.pdens.df = ve.pdens.g.aic = ve.pdens.g.p = ve.pdens.min.age = ve.pdens.d.max = array(NA,dim=c(nve,nve)) # density
for (i in 1:nve) {
  print(i)
  for (j in 1:i) {
    print(j)
    x = age.bin
    # strength
    g = gam(ve.str[i,j,] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML'); ve.str.df[i,j] = sum(g$edf); ve.str.g.p[i,j] = summary(g)$s.pv; ve.str.g.aic[i,j] = AIC(g)
    pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=F); ve.str.min.age[i,j] = age.pred[which(pred==min(pred))]
    ve.str.d.max[i,j] = (max(pred)-min(pred))/sign(age.pred[which(pred==max(pred))]-age.pred[which(pred==min(pred))])
    # density
    g = gam(ve.pdens[i,j,] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML'); ve.pdens.df[i,j] = sum(g$edf); ve.pdens.g.p[i,j] = summary(g)$s.pv; ve.pdens.g.aic[i,j] = AIC(g)
    pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=F); ve.pdens.min.age[i,j] = age.pred[which(pred==min(pred))]
    ve.pdens.d.max[i,j] = (max(pred)-min(pred))/sign(age.pred[which(pred==max(pred))]-age.pred[which(pred==min(pred))])
  }
}

# adjust p-values for multiple comparisons (FDR; first unwrap them into a vector)
ve.pdens.g.p.vec = vector(length=nve*(nve+1)/2); c = 0
for (i in 1:nve) { for (j in 1:i) { c = c+1; ve.pdens.g.p.vec[c] = ve.pdens.g.p[i,j] } }
ve.pdens.g.p.fdr = p.adjust(ve.pdens.g.p.vec, method = 'fdr')

# ve.d.max matrix (with markers)
pdf(paste(dir.path,'age_var/communities/ve_dens_d_max_mat.pdf',sep=''),width=8,height=6.62)
par(mar=c(3,3,6,8), cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
image.plot(100*ve.pdens.d.max[1:nve,rev(1:nve)],add=F,axes=F,legend.mar=8,col=bluered.col,zlim=(c(-max(abs(100*ve.pdens.d.max),na.rm=T),max(abs(100*ve.pdens.d.max),na.rm=T))))
for (j in 1:nve) {
  points((j-1)/(nve-1), 1+(1.2/nve), pch = 19, cex = 3.9, col = col.ve[j],xpd=T)
  points(1+(1.2/nve), (j-1)/(nve-1), pch = 19, cex = 3.9, col = col.ve[rev(1:nve)[j]],xpd=T)
}
c = 0
for (i in 1:nve) {
  for (j in 1:i) {
    c=c+1
    if (ve.pdens.g.p.fdr[c] < 0.05) points((i-1)/(nve-1), (rev(1:nve)[j]-1)/(nve-1), pch = 19, cex = 1, col = 'black')
    if (ve.pdens.g.p.fdr[c] < 0.01) points((i-1)/(nve-1), (rev(1:nve)[j]-1)/(nve-1), pch = 19, cex = 2, col = 'black')
  }
}
dev.off()

# ve.age.min matrix
pdf(paste(dir.path,'age_var/communities/ve_dens_min_age_mat.pdf',sep=''),width=8,height=6.62)
par(mar=c(3,3,6,8), cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
image.plot(ve.pdens.min.age[1:nve,rev(1:nve)],add=F,axes=F,legend.mar=8,col=cool,zlim=c(min(age),max(age)))
for (j in 1:nve) {
  points((j-1)/(nve-1), 1+(1.2/nve), pch = 19, cex = 3.9, col = col.ve[j],xpd=T)
  points(1+(1.2/nve), (j-1)/(nve-1), pch = 19, cex = 3.9, col = col.ve[rev(1:nve)[j]],xpd=T)
}
c = 0
for (i in 1:nve) {
  for (j in 1:i) {
    c=c+1
    if (ve.pdens.g.p.fdr[c] < 0.05) points((i-1)/(nve-1), (rev(1:nve)[j]-1)/(nve-1), pch = 19, cex = 1, col = 'black')
    if (ve.pdens.g.p.fdr[c] < 0.01) points((i-1)/(nve-1), (rev(1:nve)[j]-1)/(nve-1), pch = 19, cex = 2, col = 'black')
  }
}
dev.off()

### resting-state networks 

# average strength (str) and partial density (pdens) by yeo network
# ("partial" density here is relative to the total number of edges than can exist within a community, or between two communities)
yeo.str = yeo.pdens = array(NA,dim=c(nyeo,nyeo,nbin))
for (b in 1:nbin) {
  print(b)
  for (i in 1:nyeo) {
    for (j in 1:i) {
      if (i == j) {  # within-network (upper triangular only)
        yeo.str[i,j,b] = mean(sc.bin[yeo.id==i,yeo.id==j,b][triup[yeo.id==i,yeo.id==j]])
        yeo.pdens[i,j,b] = sum(!is.na(sc.bin.thr[yeo.id==i,yeo.id==j,b]))/(sum(yeo.id==i)*(sum(yeo.id==j)-1))
      } else {  # between-network
        yeo.str[i,j,b] = mean(sc.bin[yeo.id==i,yeo.id==j,b])
        yeo.pdens[i,j,b] = sum(!is.na(sc.bin.thr[yeo.id==i,yeo.id==j,b]))/(sum(yeo.id==i)*sum(yeo.id==j))
      }
    }
  }
}

# linear models
yeo.str.l.sl = yeo.str.l.aic = yeo.str.l.p = array(NA,dim=c(nyeo,nyeo))
yeo.str.b.l.sl = yeo.str.b.l.aic = array(NA,dim=c(nyeo,nyeo,nboot))
for (i in 1:nyeo) {
  print(i)
  for (j in 1:i) {
    print(j)
    l = lm(yeo.str[i,j,] ~ age.bin); 
    yeo.str.l.sl[i,j] = l$coefficients[2]; 
    yeo.str.l.aic[i,j] = AIC(l); 
    f = summary(l)$fstatistic; 
    yeo.str.l.p[i,j] = pf(f[1],f[2],f[3], lower.tail = FALSE)
  }
}

# splines
# strength
yeo.str.df = yeo.str.g.aic = yeo.str.g.p = array(NA,dim=c(nyeo,nyeo))
yeo.str.min.age = yeo.str.d.max = array(NA,dim=c(nyeo,nyeo))
# density
yeo.pdens.df = yeo.pdens.g.aic = yeo.pdens.g.p = array(NA,dim=c(nyeo,nyeo))
yeo.pdens.min.age = yeo.pdens.d.max = array(NA,dim=c(nyeo,nyeo))
for (i in 1:nyeo) {
  print(i)
  for (j in 1:i) {
    print(j)
    x = age.bin
    # strength
    g = gam(yeo.str[i,j,] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML'); yeo.str.df[i,j] = sum(g$edf); yeo.str.g.p[i,j] = summary(g)$s.pv; yeo.str.g.aic[i,j] = AIC(g)
    pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=F); yeo.str.min.age[i,j] = age.pred[which(pred==min(pred))]
    yeo.str.d.max[i,j] = (max(pred)-min(pred))/sign(age.pred[which(pred==max(pred))]-age.pred[which(pred==min(pred))])
    # density
    g = gam(yeo.pdens[i,j,] ~ s(x, k=kspl, m=2), min.sp=min.sp, method = 'REML'); yeo.pdens.df[i,j] = sum(g$edf); yeo.pdens.g.p[i,j] = summary(g)$s.pv; yeo.pdens.g.aic[i,j] = AIC(g)
    pred = predict(g,newdata=data.frame('x'=age.pred),se.fit=F); yeo.pdens.min.age[i,j] = age.pred[which(pred==min(pred))]
    yeo.pdens.d.max[i,j] = (max(pred)-min(pred))/sign(age.pred[which(pred==max(pred))]-age.pred[which(pred==min(pred))])
  }
}

# adjust p-values for multiple comparisons (FDR)
yeo.pdens.g.p.vec = vector(length=nyeo*(nyeo+1)/2); c = 0
for (i in 1:nve) { for (j in 1:i) { c = c+1; yeo.pdens.g.p.vec[c] = yeo.pdens.g.p[i,j] } }
yeo.pdens.g.p.fdr = p.adjust(yeo.pdens.g.p.vec, method = 'fdr')

# yeo.d.max matrix
pdf(paste(dir.path,'age_var/communities/yeo_dens_d_max_mat.pdf',sep=''),width=8,height=6.62)
par(mar=c(3,3,6,8), cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
image.plot(100*yeo.pdens.d.max[1:nyeo,rev(1:nyeo)],add=F,axes=F,legend.mar=8,col=bluered.col,zlim=(c(-max(abs(100*yeo.pdens.d.max),na.rm=T),max(abs(100*yeo.pdens.d.max),na.rm=T))))
for (j in 1:nyeo) {
  points((j-1)/(nyeo-1), 1+(1.2/nyeo), pch = 19, cex = 3.9, col = col.yeo[j],xpd=T)
  points(1+(1.2/nyeo), (j-1)/(nyeo-1), pch = 19, cex = 3.9, col = col.yeo[rev(1:nyeo)[j]],xpd=T)
}
c = 0
for (i in 1:nyeo) {
  for (j in 1:i) {
    c=c+1
    if (yeo.pdens.g.p.fdr[c] < 0.05) points((i-1)/(nyeo-1), (rev(1:nyeo)[j]-1)/(nyeo-1), pch = 19, cex = 1, col = 'black')
    if (yeo.pdens.g.p.fdr[c] < 0.01) points((i-1)/(nyeo-1), (rev(1:nyeo)[j]-1)/(nyeo-1), pch = 19, cex = 2, col = 'black')
  }
}
dev.off()

# yeo.age.min matrix
pdf(paste(dir.path,'age_var/communities/yeo_dens_min_age_mat.pdf',sep=''),width=8,height=6.62)
par(mar=c(3,3,6,8), cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
image.plot(yeo.pdens.min.age[1:nyeo,rev(1:nyeo)],add=F,axes=F,legend.mar=8,col=cool,zlim=c(min(age),max(age)))
for (j in 1:nyeo) {
  points((j-1)/(nyeo-1), 1+(1.2/nyeo), pch = 19, cex = 3.9, col = col.yeo[j],xpd=T)
  points(1+(1.2/nyeo), (j-1)/(nyeo-1), pch = 19, cex = 3.9, col = col.yeo[rev(1:nyeo)[j]],xpd=T)
}
c = 0
for (i in 1:nyeo) {
  for (j in 1:i) {
    c=c+1
    if (yeo.pdens.g.p.fdr[c] < 0.05) points((i-1)/(nyeo-1), (rev(1:nyeo)[j]-1)/(nyeo-1), pch = 19, cex = 1, col = 'black')
    if (yeo.pdens.g.p.fdr[c] < 0.01) points((i-1)/(nyeo-1), (rev(1:nyeo)[j]-1)/(nyeo-1), pch = 19, cex = 2, col = 'black')
  }
}
dev.off()

save('nmod','col.mod','mod.str','mod.pdens','mod.str.l.sl','mod.str.l.aic','mod.str.l.p','mod.str.b.l.sl','mod.str.b.l.aic','mod.str.df','mod.str.g.aic',
     'mod.str.g.p','mod.str.min.age','mod.str.d.max','mod.pdens.df','mod.pdens.g.aic','mod.pdens.g.p','mod.pdens.min.age','mod.pdens.d.max','mod.pdens.g.p.fdr',
     'nve','col.ve','ve.str','ve.pdens','ve.str.l.sl','ve.str.l.aic','ve.str.l.p','ve.str.b.l.sl','ve.str.b.l.aic','ve.str.df','ve.str.g.aic',
     've.str.g.p','ve.str.min.age','ve.str.d.max','ve.pdens.df','ve.pdens.g.aic','ve.pdens.g.p','ve.pdens.min.age','ve.pdens.d.max','ve.pdens.g.p.fdr',
     'nyeo','col.yeo','yeo.str','yeo.pdens','yeo.str.l.sl','yeo.str.l.aic','yeo.str.l.p','yeo.str.b.l.sl','yeo.str.b.l.aic','yeo.str.df','yeo.str.g.aic',
     'yeo.str.g.p','yeo.str.min.age','yeo.str.d.max','yeo.pdens.df','yeo.pdens.g.aic','yeo.pdens.g.p','yeo.pdens.min.age','yeo.pdens.d.max','yeo.pdens.g.p.fdr',
     file=paste(dir.path,'age_var/sc.age.var.communities.RData',sep=''))
     
