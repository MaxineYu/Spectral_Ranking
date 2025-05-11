### Instructions 
## AA Matrix of the comparion graph 
## WW Matrix of winner indices: encodes observed winner information 
## Output reference guide: matrix RR
##-----------------------------------------------------------------    ## Table EC.8 2011-2015
## RR[1,] theta.hat
## RR[2,] rank based on theta hat
## RR{3,] & RR[4,] two-sided CI for rank                               ## RR[1,] to RR[6,] are based on the Vanilla spectral method
## RR[5,] left-sided CI for rank
## RR[6,] uniform left-sided CI for rank
##-----------------------------------------------------------------    ## Table 1 2011-2015
## RR[7,] theta.hat
## RR[8,] rank based on theta hat
## RR{9,] & RR[10,] two-sided CI for rank                              ## RR[7,] to RR[12,] are based on the two-stage theta estimators 
## RR[11,] left-sided CI for rank
## RR[12,] uniform left-sided CI for rank







## Refer to Journal_Citation_Data_Input.R for Ranking_AA/Ranking_Journal/Ranking_ToYear/Ranking_Year/Ranking_WW ## 





AA0 = read.table('Ranking_AA')
AA0 = as.matrix(AA0)
IdxJournal = read.table('Ranking_Journal')
IdxYearTo = read.table('Ranking_ToYear')
IdxYear = read.table('Ranking_Year')
WW0 = read.table('Ranking_WW')
WW0 = as.matrix(WW0)









##-----------------------------------------------------------------    ## generate comparison graph AA & top choice WW ##
NullIdx = c(13,14,15)
tmp.idx2 = (IdxYear$x > 2010 & (IdxYear$x - IdxYearTo$x) < 10)
AAA = AA0[tmp.idx2,]
WWW = WW0[tmp.idx2,]
tmp.I = (AAA[,NullIdx] == 1)
tmp.R = rowSums(tmp.I)
tmp.R = (tmp.R == 0)
AA = AAA[tmp.R,-NullIdx]                                               ## comparison graph ## 
WW = WWW[tmp.R,-NullIdx]                                               ## top choice ## 








##-----------------------------------------------------------------    ## Vanilla spectral method ##
B = 2000
n = ncol(AA)
L = nrow(AA)
fAvec = numeric(L)+2                                                   ## weight for vanilla spectral method ##




##-----------------------------------------------------------------    ## compute matrix P ##
dval = 2*max(colSums(AA))
P = matrix(0,n,n)
for(i in 1:n){
	for(j in 1:n){
		if(j != i){
			P[i,j] = sum(AA[,i]*AA[,j]*WW[,j]/fAvec)/dval        ## dval ##
		}
	}
	P[i,i] = 1-sum(P[i,])
}




##-----------------------------------------------------------------    ## solve theta and pi ##
tmp.P = t(t(P)-diag(n))%*%(t(P)-diag(n))
tmp.svd = svd(tmp.P)
pihat = abs(tmp.svd$v[,n])
thetahat = log(pihat)-mean(log(pihat))




##-----------------------------------------------------------------    ## output rank ##
RR = matrix(0,12,n)
colnames(RR) = IdxJournal$x[-NullIdx]
RR[1,] = thetahat
RR[2,] = n+1-rank(thetahat)




##-----------------------------------------------------------------    ## compute sample xihat ##   
Vmatrix = matrix(0,L,n)
tauhatvec = numeric(n)
tmp.pimatrix = t(AA)*pihat
tmp.pivec = colSums(tmp.pimatrix)
tmp.var = numeric(n)
for(oo in 1:n){
	tauhatvec[oo] = sum(AA[,oo]*(1-pihat[oo]/tmp.pivec)*pihat[oo]/fAvec)/dval
	tmp.var[oo] = sum(AA[,oo]*(tmp.pivec-pihat[oo])/fAvec/fAvec)*pihat[oo]/dval/dval/tauhatvec[oo]/tauhatvec[oo]
	Vmatrix[,oo] = (AA[,oo]*WW[,oo]*tmp.pivec-AA[,oo]*pihat[oo])/fAvec  
}
sigmahatmatrix = matrix(tmp.var,n,n)+t(matrix(tmp.var,n,n))




##-----------------------------------------------------------------    ## Weighted bootstrap ##
Wmatrix = matrix(rnorm(L*B),L,B)
tmp.Vtau = (t(Vmatrix)/tauhatvec)%*%Wmatrix




R.left.m = numeric(n)
R.right.m = numeric(n)
R.left.one.m = numeric(n)
for(ooo in 1:n){
	print(ooo)
	tmpGMmatrix0 = matrix(rep(tmp.Vtau[ooo,],n)-c(t(tmp.Vtau)),B,n)
	tmpGMmatrix = abs(t(t(tmpGMmatrix0)/sqrt(sigmahatmatrix[ooo,]))/dval)
	tmpGMmatrixone = t(t(tmpGMmatrix0)/sqrt(sigmahatmatrix[ooo,]))/dval
	tmp.GMvecmax = apply(tmpGMmatrix,1,max)
	tmp.GMvecmaxone = apply(tmpGMmatrixone,1,max)
    cutval = quantile(tmp.GMvecmax,0.95)
    cutvalone = quantile(tmp.GMvecmaxone,0.95)
	tmp.theta.sd = sqrt(sigmahatmatrix[ooo,])
	tmp.theta.sd = tmp.theta.sd[-ooo]
	R.left.m[ooo] = 1+sum(1*(((thetahat[-ooo]-thetahat[ooo])/tmp.theta.sd)>cutval))
	R.right.m[ooo] = n-sum(1*(((thetahat[-ooo]-thetahat[ooo])/tmp.theta.sd)<(-cutval)))
	R.left.one.m[ooo] = 1+sum(1*(((thetahat[-ooo]-thetahat[ooo])/tmp.theta.sd)>cutvalone))
}
## two-sided CI for rank ##
RR[3,] = R.left.m
RR[4,] = R.right.m

## left-sided CI for rank ##
RR[5,] = R.left.one.m










##-----------------------------------------------------------------    ## Weighted bootstrap ##
Wmatrix = matrix(rnorm(L*B),L,B)
tmp.Vtau = (t(Vmatrix)/tauhatvec)%*%Wmatrix
Mval = n
GMvecmax = numeric(B)-1
GMvecmaxone = numeric(B)-Inf
tmpTMval = -1
tmpTMvalone = -Inf





##-----------------------------------------------------------------    ##
for(ooo in 1:n){
	tmpGMmatrix0 = matrix(rep(tmp.Vtau[ooo,],n)-c(t(tmp.Vtau)),B,n)
	tmpGMmatrix = abs(t(t(tmpGMmatrix0)/sqrt(sigmahatmatrix[ooo,]))/dval)
	tmpGMmatrixone = t(t(tmpGMmatrix0)/sqrt(sigmahatmatrix[ooo,]))/dval
	tmp.GMvecmax = apply(tmpGMmatrix,1,max)
	tmp.GMvecmaxone = apply(tmpGMmatrixone,1,max)
	GMvecmax = c(GMvecmax,tmp.GMvecmax)
	GMvecmaxone = c(GMvecmaxone,tmp.GMvecmaxone)
}
GMmaxmatrixone = matrix(GMvecmaxone,B)
GMmaxone = apply(GMmaxmatrixone,1,max)
cutvalone = quantile(GMmaxone,0.95)







##-----------------------------------------------------------------    ##
R.left.one = numeric(n)
for(oooo in 1:n){
	tmp.theta.sd = sqrt(sigmahatmatrix[oooo,])
	tmp.theta.sd = tmp.theta.sd[-oooo]
	R.left.one[oooo] = 1+sum(1*(((thetahat[-oooo]-thetahat[oooo])/tmp.theta.sd)>cutvalone))
}
RR[6,] = R.left.one                                                    ## uniform left-sided CI for rank ##

















##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ## Two-stage spectal method ##
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ## 
fMLE = numeric(L)
ftwo = rowSums(AA*t(matrix(exp(thetahat),n,nrow(AA))))
fMLE = ftwo                                                            ## weight for two-stage spectral method ##




##-----------------------------------------------------------------    ## compute matrix P ##
dval = 2*max(colSums(AA))
PMLE = matrix(0,n,n)
for(i in 1:n){
	for(j in 1:n){
		if(j != i){
			PMLE[i,j] = sum(AA[,i]*AA[,j]*WW[,j]/fMLE)/dval
		}
	}
	PMLE[i,i] = 1-sum(PMLE[i,])
}




##-----------------------------------------------------------------    ## solve theta and pi ##
tmp.PMLE = t(t(PMLE)-diag(n))%*%(t(PMLE)-diag(n))
tmp.svd.MLE = svd(tmp.PMLE)
pihatMLE = abs(tmp.svd.MLE$v[,n])
thetahatMLE = log(pihatMLE)-mean(log(pihatMLE))
RR[7,] = thetahatMLE
RR[8,] = n+1-rank(thetahatMLE)




##-----------------------------------------------------------------    ## compute sample xihat ##
pihat = pihatMLE
thetahat = thetahatMLE
P = PMLE
fAvec = fMLE




##-----------------------------------------------------------------    ##
Vmatrix = matrix(0,L,n)
tauhatvec = numeric(n)
tmp.pimatrix = t(AA)*pihat
tmp.pivec = colSums(tmp.pimatrix)
tmp.var = numeric(n)
for(oo in 1:n){
	tauhatvec[oo] = sum(AA[,oo]*(1-pihat[oo]/tmp.pivec)*pihat[oo]/fAvec)/dval
	tmp.var[oo] = sum(AA[,oo]*(tmp.pivec-pihat[oo])/fAvec/fAvec)*pihat[oo]/dval/dval/tauhatvec[oo]/tauhatvec[oo]
	Vmatrix[,oo] = (AA[,oo]*WW[,oo]*tmp.pivec-AA[,oo]*pihat[oo])/fAvec  
}
sigmahatmatrix = matrix(tmp.var,n,n)+t(matrix(tmp.var,n,n))




##-----------------------------------------------------------------    ## Weighted bootstrap ## 
Wmatrix = matrix(rnorm(L*B),L,B)
tmp.Vtau = (t(Vmatrix)/tauhatvec)%*%Wmatrix




##-----------------------------------------------------------------    ##
R.left.m = numeric(n)
R.right.m = numeric(n)
R.left.one.m = numeric(n)
for(ooo in 1:n){
	print(ooo)
	tmpGMmatrix0 = matrix(rep(tmp.Vtau[ooo,],n)-c(t(tmp.Vtau)),B,n)
	tmpGMmatrix = abs(t(t(tmpGMmatrix0)/sqrt(sigmahatmatrix[ooo,]))/dval)
	tmpGMmatrixone = t(t(tmpGMmatrix0)/sqrt(sigmahatmatrix[ooo,]))/dval
	tmp.GMvecmax = apply(tmpGMmatrix,1,max)
	tmp.GMvecmaxone = apply(tmpGMmatrixone,1,max)
    cutval = quantile(tmp.GMvecmax,0.95)
    cutvalone = quantile(tmp.GMvecmaxone,0.95)
	tmp.theta.sd = sqrt(sigmahatmatrix[ooo,])
	tmp.theta.sd = tmp.theta.sd[-ooo]
	R.left.m[ooo] = 1+sum(1*(((thetahat[-ooo]-thetahat[ooo])/tmp.theta.sd)>cutval))
	R.right.m[ooo] = n-sum(1*(((thetahat[-ooo]-thetahat[ooo])/tmp.theta.sd)<(-cutval)))
	R.left.one.m[ooo] = 1+sum(1*(((thetahat[-ooo]-thetahat[ooo])/tmp.theta.sd)>cutvalone))
}
RR[9,] = R.left.m
RR[10,] = R.right.m
RR[11,] = R.left.one.m




##-----------------------------------------------------------------    ## Weighted bootstrap##
Wmatrix = matrix(rnorm(L*B),L,B)
tmp.Vtau = (t(Vmatrix)/tauhatvec)%*%Wmatrix
Mval = n
GMvecmax = numeric(B)-1
GMvecmaxone = numeric(B)-Inf
tmpTMval = -1
tmpTMvalone = -Inf




##-----------------------------------------------------------------    ##
for(ooo in 1:n){
	tmpGMmatrix0 = matrix(rep(tmp.Vtau[ooo,],n)-c(t(tmp.Vtau)),B,n)
	tmpGMmatrix = abs(t(t(tmpGMmatrix0)/sqrt(sigmahatmatrix[ooo,]))/dval)
	tmpGMmatrixone = t(t(tmpGMmatrix0)/sqrt(sigmahatmatrix[ooo,]))/dval
	tmp.GMvecmax = apply(tmpGMmatrix,1,max)
	tmp.GMvecmaxone = apply(tmpGMmatrixone,1,max)
	GMvecmax = c(GMvecmax,tmp.GMvecmax)
	GMvecmaxone = c(GMvecmaxone,tmp.GMvecmaxone)
}
GMmaxmatrixone = matrix(GMvecmaxone,B)
GMmaxone = apply(GMmaxmatrixone,1,max)
cutvalone = quantile(GMmaxone,0.95)
R.left.one = numeric(n)
for(oooo in 1:n){
	tmp.theta.sd = sqrt(sigmahatmatrix[oooo,])
	tmp.theta.sd = tmp.theta.sd[-oooo]
	R.left.one[oooo] = 1+sum(1*(((thetahat[-oooo]-thetahat[oooo])/tmp.theta.sd)>cutvalone))
}
RR[12,] = R.left.one                                                   ## uniform left-sided CI for rank ##

