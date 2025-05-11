### Instructions 
## AA Matrix of the comparion graph 
## WW Matrix of winner indices: encodes observed winner information 
## Output reference guide: matrix RR
##-----------------------------------------------------------------    ## Table EC.8 2006-2010
## RR[1,] theta.hat
## RR[2,] rank based on theta hat
## RR{3,] & RR[4,] two-sided CI for rank                               ## RR[1,] to RR[6,] are based on the Vanilla spectral method
## RR[5,] left-sided CI for rank
## RR[6,] uniform left-sided CI for rank
##-----------------------------------------------------------------    ## Table 1 2006-2010
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




##-----------------------------------------------------------------    ##
NullIdx = c(13,14,15)
tmp.idx2 = (IdxYear$x > 2005 & IdxYear$x <= 2010 & (IdxYear$x - IdxYearTo$x) < 10)
AAA = AA0[tmp.idx2,]
WWW = WW0[tmp.idx2,]



##-----------------------------------------------------------------    ## generate comparison graph AA & top choice WW ##
tmp.I = (AAA[,NullIdx] == 1)
tmp.R = rowSums(tmp.I)
tmp.R = (tmp.R == 0)
AA2 = AAA[tmp.R,-NullIdx]                                              ## comparison graph ##
WW2 = WWW[tmp.R,-NullIdx]                                              ## top choice ## 











##-----------------------------------------------------------------    ## Vanilla spectral method ##
B = 2000
n = ncol(AA2)
L2 = nrow(AA2)
fAvec2 = numeric(L2)+2                                                 ## weight for vanilla spectral method ##



##-----------------------------------------------------------------    ## compute matrix P ##
dval2 = 2*max(colSums(AA2))                                            ## value of d ##
P2 = matrix(0,n,n)
for(i in 1:n){
	for(j in 1:n){
		if(j != i){
			P2[i,j] = sum(AA2[,i]*AA2[,j]*WW2[,j]/fAvec2)/dval2        
		}
	}
	P2[i,i] = 1-sum(P2[i,])
}



##-----------------------------------------------------------------    ## solve theta and pi ##
tmp.P2 = t(t(P2)-diag(n))%*%(t(P2)-diag(n))
tmp.svd2 = svd(tmp.P2)
pihat2 = abs(tmp.svd2$v[,n])
thetahat2 = log(pihat2)-mean(log(pihat2))



##-----------------------------------------------------------------    ## output ##
RR2 = matrix(0,12,n)
colnames(RR2) = IdxJournal$x[-NullIdx]
RR2[1,] = thetahat2
RR2[2,] = n+1-rank(thetahat2)   



##-----------------------------------------------------------------       
Vmatrix2 = matrix(0,L2,n)
tauhatvec2 = numeric(n)
tmp.pimatrix2 = t(AA2)*pihat2
tmp.pivec2 = colSums(tmp.pimatrix2)
tmp.var2 = numeric(n)
for(oo in 1:n){
	tauhatvec2[oo] = sum(AA2[,oo]*(1-pihat2[oo]/tmp.pivec2)*pihat2[oo]/fAvec2)/dval2
	tmp.var2[oo] = sum(AA2[,oo]*(tmp.pivec2-pihat2[oo])/fAvec2/fAvec2)*pihat2[oo]/dval2/dval2/tauhatvec2[oo]/tauhatvec2[oo]
	Vmatrix2[,oo] = (AA2[,oo]*WW2[,oo]*tmp.pivec2-AA2[,oo]*pihat2[oo])/fAvec2  
}
sigmahatmatrix2 = matrix(tmp.var2,n,n)+t(matrix(tmp.var2,n,n))



##-----------------------------------------------------------------    ## Weighted bootstrap ##
Wmatrix2 = matrix(rnorm(L2*B),L2,B)
tmp.Vtau2 = (t(Vmatrix2)/tauhatvec2)%*%Wmatrix2







##-----------------------------------------------------------------    ##
R.left.m2 = numeric(n)
R.right.m2 = numeric(n)
R.left.one.m2 = numeric(n)
for(ooo in 1:n){
	print(ooo)
	tmpGMmatrix02 = matrix(rep(tmp.Vtau2[ooo,],n)-c(t(tmp.Vtau2)),B,n)
	tmpGMmatrix2 = abs(t(t(tmpGMmatrix02)/sqrt(sigmahatmatrix2[ooo,]))/dval2)
	tmpGMmatrixone2 = t(t(tmpGMmatrix02)/sqrt(sigmahatmatrix2[ooo,]))/dval2
	tmp.GMvecmax2 = apply(tmpGMmatrix2,1,max)
	tmp.GMvecmaxone2 = apply(tmpGMmatrixone2,1,max)
    cutval2 = quantile(tmp.GMvecmax2,0.95)
    cutvalone2 = quantile(tmp.GMvecmaxone2,0.95)
	tmp.theta.sd2 = sqrt(sigmahatmatrix2[ooo,])
	tmp.theta.sd2 = tmp.theta.sd2[-ooo]
	R.left.m2[ooo] = 1+sum(1*(((thetahat2[-ooo]-thetahat2[ooo])/tmp.theta.sd2)>cutval2))
	R.right.m2[ooo] = n-sum(1*(((thetahat2[-ooo]-thetahat2[ooo])/tmp.theta.sd2)<(-cutval2)))
	R.left.one.m2[ooo] = 1+sum(1*(((thetahat2[-ooo]-thetahat2[ooo])/tmp.theta.sd2)>cutvalone2))
}
## two-sided CI for rank ##
RR2[3,] = R.left.m2
RR2[4,] = R.right.m2
## left-sided CI for rank ##
RR2[5,] = R.left.one.m2







##-----------------------------------------------------------------    ## Weighted bootstrap ##
Wmatrix2 = matrix(rnorm(L2*B),L2,B)
tmp.Vtau2 = (t(Vmatrix2)/tauhatvec2)%*%Wmatrix2
Mval = n
GMvecmax2 = numeric(B)-1
GMvecmaxone2 = numeric(B)-Inf
tmpTMval2 = -1
tmpTMvalone2 = -Inf




##-----------------------------------------------------------------    ##
for(ooo in 1:n){
	tmpGMmatrix02 = matrix(rep(tmp.Vtau2[ooo,],n)-c(t(tmp.Vtau2)),B,n)
	tmpGMmatrix2 = abs(t(t(tmpGMmatrix02)/sqrt(sigmahatmatrix2[ooo,]))/dval2)
	tmpGMmatrixone2 = t(t(tmpGMmatrix02)/sqrt(sigmahatmatrix2[ooo,]))/dval2
	tmp.GMvecmax2 = apply(tmpGMmatrix2,1,max)
	tmp.GMvecmaxone2 = apply(tmpGMmatrixone2,1,max)
	GMvecmax2 = c(GMvecmax2,tmp.GMvecmax2)
	GMvecmaxone2 = c(GMvecmaxone2,tmp.GMvecmaxone2)
}




##-----------------------------------------------------------------    ##
GMmaxmatrixone2 = matrix(GMvecmaxone2,B)
GMmaxone2 = apply(GMmaxmatrixone2,1,max)
cutvalone2 = quantile(GMmaxone2,0.95)





##-----------------------------------------------------------------    ##
R.left.one2 = numeric(n)
for(oooo in 1:n){
	tmp.theta.sd2 = sqrt(sigmahatmatrix2[oooo,])
	tmp.theta.sd2 = tmp.theta.sd2[-oooo]
	R.left.one2[oooo] = 1+sum(1*(((thetahat2[-oooo]-thetahat2[oooo])/tmp.theta.sd2)>cutvalone2))
}
RR2[6,] = R.left.one2                                                  ## uniform left-sided CI for rank ##


















##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ## Two-stage spectal method ##
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ## 
fMLE2 = numeric(L2)
ftwo2 = rowSums(AA2*t(matrix(exp(thetahat2),n,nrow(AA2))))
fMLE2 = ftwo2                                                          ## weight for two-stage spectral method ##





##-----------------------------------------------------------------    ## compute matrix P ##
dval2 = 2*max(colSums(AA2))
PMLE2 = matrix(0,n,n)
for(i in 1:n){
	for(j in 1:n){
		if(j != i){
			PMLE2[i,j] = sum(AA2[,i]*AA2[,j]*WW2[,j]/fMLE2)/dval2
		}
	}
	PMLE2[i,i] = 1-sum(PMLE2[i,])
}




##-----------------------------------------------------------------    ## solve theta and pi ##
tmp.PMLE2 = t(t(PMLE2)-diag(n))%*%(t(PMLE2)-diag(n))
tmp.svd.MLE2 = svd(tmp.PMLE2)
pihatMLE2 = abs(tmp.svd.MLE2$v[,n])
thetahatMLE2 = log(pihatMLE2)-mean(log(pihatMLE2))
RR2[7,] = thetahatMLE2                                                 ## theta hat ##
RR2[8,] = n+1-rank(thetahatMLE2)                                       ## rank based on theta hat ##





##-----------------------------------------------------------------    ## compute sample xihat ##
pihat2 = pihatMLE2
thetahat2 = thetahatMLE2
P2 = PMLE2
fAvec2 = fMLE2
Vmatrix2 = matrix(0,L2,n)
tauhatvec2 = numeric(n)
tmp.pimatrix2 = t(AA2)*pihat2
tmp.pivec2 = colSums(tmp.pimatrix2)
tmp.var2 = numeric(n)
for(oo in 1:n){
	tauhatvec2[oo] = sum(AA2[,oo]*(1-pihat2[oo]/tmp.pivec2)*pihat2[oo]/fAvec2)/dval2
	tmp.var2[oo] = sum(AA2[,oo]*(tmp.pivec2-pihat2[oo])/fAvec2/fAvec2)*pihat2[oo]/dval2/dval2/tauhatvec2[oo]/tauhatvec2[oo]
	Vmatrix2[,oo] = (AA2[,oo]*WW2[,oo]*tmp.pivec2-AA2[,oo]*pihat2[oo])/fAvec2  
}
sigmahatmatrix2 = matrix(tmp.var2,n,n)+t(matrix(tmp.var2,n,n))




##-----------------------------------------------------------------    ## Weighted bootstrap ##
Wmatrix2 = matrix(rnorm(L2*B),L2,B)
tmp.Vtau2 = (t(Vmatrix2)/tauhatvec2)%*%Wmatrix2





##-----------------------------------------------------------------    ##
R.left.m2 = numeric(n)
R.right.m2 = numeric(n)
R.left.one.m2 = numeric(n)
for(ooo in 1:n){
	print(ooo)
	tmpGMmatrix02 = matrix(rep(tmp.Vtau2[ooo,],n)-c(t(tmp.Vtau2)),B,n)
	tmpGMmatrix2 = abs(t(t(tmpGMmatrix02)/sqrt(sigmahatmatrix2[ooo,]))/dval2)
	tmpGMmatrixone2 = t(t(tmpGMmatrix02)/sqrt(sigmahatmatrix2[ooo,]))/dval2
	tmp.GMvecmax2 = apply(tmpGMmatrix2,1,max)
	tmp.GMvecmaxone2 = apply(tmpGMmatrixone2,1,max)
    cutval2 = quantile(tmp.GMvecmax2,0.95)
    cutvalone2 = quantile(tmp.GMvecmaxone2,0.95)
	tmp.theta.sd2 = sqrt(sigmahatmatrix2[ooo,])
	tmp.theta.sd2 = tmp.theta.sd2[-ooo]
	R.left.m2[ooo] = 1+sum(1*(((thetahat2[-ooo]-thetahat2[ooo])/tmp.theta.sd2)>cutval2))
	R.right.m2[ooo] = n-sum(1*(((thetahat2[-ooo]-thetahat2[ooo])/tmp.theta.sd2)<(-cutval2)))
	R.left.one.m2[ooo] = 1+sum(1*(((thetahat2[-ooo]-thetahat2[ooo])/tmp.theta.sd2)>cutvalone2))
}
## two-sided CI for rank ##
RR2[9,] = R.left.m2
RR2[10,] = R.right.m2
## left-sided CI for rank ##
RR2[11,] = R.left.one.m2






##-----------------------------------------------------------------    ## Weighted bootstrap ##
Wmatrix2 = matrix(rnorm(L2*B),L2,B)
tmp.Vtau2 = (t(Vmatrix2)/tauhatvec2)%*%Wmatrix2
Mval = n
GMvecmax2 = numeric(B)-1
GMvecmaxone2 = numeric(B)-Inf
tmpTMval2 = -1
tmpTMvalone2 = -Inf




##-----------------------------------------------------------------    ##
for(ooo in 1:n){
	tmpGMmatrix02 = matrix(rep(tmp.Vtau2[ooo,],n)-c(t(tmp.Vtau2)),B,n)
	tmpGMmatrix2 = abs(t(t(tmpGMmatrix02)/sqrt(sigmahatmatrix2[ooo,]))/dval2)
	tmpGMmatrixone2 = t(t(tmpGMmatrix02)/sqrt(sigmahatmatrix2[ooo,]))/dval2
	tmp.GMvecmax2 = apply(tmpGMmatrix2,1,max)
	tmp.GMvecmaxone2 = apply(tmpGMmatrixone2,1,max)
	GMvecmax2 = c(GMvecmax2,tmp.GMvecmax2)
	GMvecmaxone2 = c(GMvecmaxone2,tmp.GMvecmaxone2)
}
GMmaxmatrixone2 = matrix(GMvecmaxone2,B)
GMmaxone2 = apply(GMmaxmatrixone2,1,max)
cutvalone2 = quantile(GMmaxone2,0.95)
R.left.one2 = numeric(n)
for(oooo in 1:n){
	tmp.theta.sd2 = sqrt(sigmahatmatrix2[oooo,])
	tmp.theta.sd2 = tmp.theta.sd2[-oooo]
	R.left.one2[oooo] = 1+sum(1*(((thetahat2[-oooo]-thetahat2[oooo])/tmp.theta.sd2)>cutvalone2))
}
RR2[12,] = R.left.one2                                                 ## uniform left-sided CI for rank ##


