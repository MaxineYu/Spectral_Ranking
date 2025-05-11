### Instructions 
## AA Matrix of the comparion graph 
## WW Matrix of winner indices: encodes observed winner information 
## Output reference guide: matrix RR
##-----------------------------------------------------------------    ## Table EC.9
## RR[1,] theta.hat
## RR[2,] rank based on theta hat
## RR{3,] & RR[4,] two-sided CI for rank                               ## RR[1,] to RR[6,] are based on the Vanilla spectral method
## RR[5,] left-sided CI for rank
## RR[6,] uniform left-sided CI for rank
##-----------------------------------------------------------------    ## Table 2
## RR[7,] theta.hat
## RR[8,] rank based on theta hat
## RR{9,] & RR[10,] two-sided CI for rank                              ## RR[7,] to RR[12,] are based on the two-stage theta estimators 
## RR[11,] left-sided CI for rank
## RR[12,] uniform left-sided CI for rank






library(Matrix)
XXX = read.csv('movie_ranking_table.csv')
XXX = as.matrix(XXX)
XXX = XXX[,-1]
Idx = read.csv('movie_name_by_id.csv')
Idx = as.matrix(Idx)
Num = read.csv('movie_ranking_count.csv')
Num = as.matrix(Num)





##-----------------------------------------------------------------    ## generate comparison graph AA & top choice WW ##
AA0 = matrix(0,196,0)
for(ss in 1:3000){
	print(ss)
	if(Num[ss,2]>0){
		AA0 = cbind(AA0,matrix(rep(XXX[ss,],Num[ss,2]),nrow=196))
	}
}
AA0 = as.matrix(t(AA0))
WW0 = as.matrix(1*(AA0 == 1))
AA0 = as.matrix(1*(AA0 > 0))

AA = as(AA0,'sparseMatrix')  ## comparison graph ##                                             
WW = as(WW0,'sparseMatrix')  ## top choice ##
 







##-----------------------------------------------------------------    ## Vanilla spectral method ##
B = 2000                                                               
n = ncol(AA)  
L = nrow(AA)   
fAvec = rowSums(AA>0)                                                  ## weight for vanilla spectral method ##



##-----------------------------------------------------------------    ## compute matrix P ##
dval = 2*max(colSums(AA))  ## value of d ##
P = matrix(0,n,n)
for(i in 1:n){
	print(i)
	for(j in 1:n){
		if(j != i){
			P[i,j] = sum(AA[,i]*AA[,j]*WW[,j]/fAvec)/dval              
		}
	}
	P[i,i] = 1-sum(P[i,])
}



##-----------------------------------------------------------------    ## solve theta and pi ##
tmp.P = t(t(P)-diag(n))%*%(t(P)-diag(n))
tmp.svd = svd(tmp.P)
pihat = abs(tmp.svd$v[,n])   ## pi hat##
thetahat = log(pihat)-mean(log(pihat))    ## theta hat##



##-----------------------------------------------------------------    ## output ##
RR = matrix(0,12,n)
colnames(RR) = Idx[,2]
RR[1,] = thetahat   ## theta hat ## 
RR[2,] = n+1-rank(thetahat)   ## rank based on theta hat ##     


                                                        
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
tmp.matrixV = t(Vmatrix)/tauhatvec
tmp.Vtau = matrix(0,n,B)
for(zzz in 1:B){
	print(zzz)
	tmp.Vtau[,zzz] = colSums(t(tmp.matrixV)*rnorm(L))
}



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
    cutval = quantile(tmp.GMvecmax,0.95)                                ## two-sided critical value ##
    cutvalone = quantile(tmp.GMvecmaxone,0.95)                          ## one-sided critical value ## 
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
tmp.matrixV = t(Vmatrix)/tauhatvec
tmp.Vtau = matrix(0,n,B)
for(zzz in 1:B){
	print(zzz)
	tmp.Vtau[,zzz] = colSums(t(tmp.matrixV)*rnorm(L))
}
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



##-----------------------------------------------------------------    ##
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
RR[6,] = R.left.one                                                   ## uniform left-sided CI for rank ##














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
dval = 2*max(colSums(AA))                                              ## value of d ##
PMLE = matrix(0,n,n)
for(i in 1:n){
	print(i)
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
RR[7,] = thetahatMLE                                                   ## theta hat ##
RR[8,] = n+1-rank(thetahatMLE)                                         ## rank based on theta hat ##





##-----------------------------------------------------------------    ## compute sample xihat ##
pihat = pihatMLE
thetahat = thetahatMLE
P = PMLE
fAvec = fMLE
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
tmp.matrixV = t(Vmatrix)/tauhatvec
tmp.Vtau = matrix(0,n,B)
for(zzz in 1:B){
	print(zzz)
	tmp.Vtau[,zzz] = colSums(t(tmp.matrixV)*rnorm(L))
}
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
RR[9,] = R.left.m
RR[10,] = R.right.m
## left-sided CI for rank ##
RR[11,] = R.left.one.m



##-----------------------------------------------------------------    ## Weighted bootstrap ##
Wmatrix = matrix(rnorm(L*B),L,B)
tmp.Vtau = (t(Vmatrix)/tauhatvec)%*%Wmatrix
Mval = n
GMvecmax = numeric(B)-1
GMvecmaxone = numeric(B)-Inf
tmpTMval = -1
tmpTMvalone = -Inf
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























