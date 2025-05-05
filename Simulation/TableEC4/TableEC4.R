##-----------------------------------------------------------------    ## #########Uncomment Line 146-149 to get the results for the oracle estimator#####
NN = 500; B = 300; n = 50; n1 = n*0.2; n2 = n*0.5
L = 24000
L1 = round(L/5); L2 = round(L*2/5)
LL = 27000
LL1 = round(LL/5); LL2 = round(LL*2/5)
##-----------------------------------------------------------------    ##
setM = c(10)         ## item of interest ##
##-----------------------------------------------------------------    ##
theta1 = seq(6,2,length=n); theta1 = theta1-mean(theta1) 
theta2 = seq(6,2,length=n); theta2 = theta2-mean(theta2)                                           ## null hypothesis  ##
#theta2 = theta1; theta2[1:9] = theta2[1:9]+0.2; theta2[11:n] = theta2[11:n]-0.2*9/40             ## null hypothesis  ## ##uncomment this when null holds###
theta2 = theta1; theta2[10] = theta1[13]; theta2[13] = theta1[10]                                ## alternative hypothesis  ## 
##-----------------------------------------------------------------    ##
statone = numeric(NN)
##statMLEone = numeric(NN)
##-----------------------------------------------------------------    ##
for(ww in 1:NN){
print(ww)
fAvec = numeric(L)
fMLE = numeric(L)
AA = matrix(0,L,n)
WW = matrix(0,L,n)
for(o in 1:L1){
	tmp.M = sample(c(2,3,4,5),1)                                     ## 
	fAvec[o] = tmp.M
	tmp.idx = sample(1:n1,tmp.M)
	AA[o,tmp.idx] = 1
	tmp.w = exp(theta1[tmp.idx])
	WW[o,tmp.idx] = rmultinom(1,1,tmp.w) 
	fMLE[o] = sum(tmp.w)
}
for(o in (L1+1):L2){
	tmp.M = sample(c(2,3,4,5),1)                                     ## 
	fAvec[o] = tmp.M
	tmp.idx = sample(1:n2,tmp.M)
	AA[o,tmp.idx] = 1
	tmp.w = exp(theta1[tmp.idx])
	WW[o,tmp.idx] = rmultinom(1,1,tmp.w) 
	fMLE[o] = sum(tmp.w)
}
for(o in (L2+1):L){
	tmp.M = sample(c(2,3,4,5),1)                                     ## 
	fAvec[o] = tmp.M
	tmp.idx = sample(1:n,tmp.M)
	AA[o,tmp.idx] = 1
	tmp.w = exp(theta1[tmp.idx])
	WW[o,tmp.idx] = rmultinom(1,1,tmp.w) 
	fMLE[o] = sum(tmp.w)
}
##-----------------------------------------------------------------    ##
fAvec2 = numeric(LL)
fMLE2 = numeric(LL)
AA2 = matrix(0,LL,n)
WW2 = matrix(0,LL,n)
for(o in 1:LL1){
	tmp.M = sample(c(2,3,4,5),1)                                     ## 
	fAvec2[o] = tmp.M
	tmp.idx = sample(1:n1,tmp.M)
	AA2[o,tmp.idx] = 1
	tmp.w = exp(theta2[tmp.idx])
	WW2[o,tmp.idx] = rmultinom(1,1,tmp.w) 
	fMLE2[o] = sum(tmp.w)
}
for(o in (LL1+1):LL2){
	tmp.M = sample(c(2,3,4,5),1)                                     ## 
	fAvec2[o] = tmp.M
	tmp.idx = sample(1:n2,tmp.M)
	AA2[o,tmp.idx] = 1
	tmp.w = exp(theta2[tmp.idx])
	WW2[o,tmp.idx] = rmultinom(1,1,tmp.w) 
	fMLE2[o] = sum(tmp.w)
}
for(o in (LL2+1):LL){
	tmp.M = sample(c(2,3,4,5),1)                                     ## 
	fAvec2[o] = tmp.M
	tmp.idx = sample(1:n,tmp.M)
	AA2[o,tmp.idx] = 1
	tmp.w = exp(theta2[tmp.idx])
	WW2[o,tmp.idx] = rmultinom(1,1,tmp.w) 
	fMLE2[o] = sum(tmp.w)
}
##-----------------------------------------------------------------    ## compute matrix P ##
dval = 2*max(colSums(AA))
P = matrix(0,n,n)
PMLE = matrix(0,n,n)
for(i in 1:n){
	for(j in 1:n){
		if(j != i){
			P[i,j] = sum(AA[,i]*AA[,j]*WW[,j]/fAvec)/dval        ## dval ##
			PMLE[i,j] = sum(AA[,i]*AA[,j]*WW[,j]/fMLE)/dval
		}
	}
	P[i,i] = 1-sum(P[i,])
	PMLE[i,i] = 1-sum(PMLE[i,])
}
##-----------------------------------------------------------------    ## solve theta and pi ##
tmp.P = t(t(P)-diag(n))%*%(t(P)-diag(n))
tmp.svd = svd(tmp.P)
pihat = abs(tmp.svd$v[,n])
thetahat = log(pihat)-mean(log(pihat))
##Theta.matrix[,ww] = thetahat
##sum((thetahat-thetavec)^2)
##-----------------------------------------------------------------    ## solve oracle MLE ##
tmp.PMLE = t(t(PMLE)-diag(n))%*%(t(PMLE)-diag(n))
tmp.svd.MLE = svd(tmp.PMLE)
pihatMLE = abs(tmp.svd.MLE$v[,n])
thetahatMLE = log(pihatMLE)-mean(log(pihatMLE))
##Theta.oracle[,ww] = thetahatMLE
##sum((thetahatMLE-thetavec)^2)
##-----------------------------------------------------------------    ## compute sample xihat ##
##-----------------------------------------------------------------    
##-----------------------------------------------------------------    ## compute matrix P ##
dval2 = 2*max(colSums(AA2))
P2 = matrix(0,n,n)
PMLE2 = matrix(0,n,n)
for(i in 1:n){
	for(j in 1:n){
		if(j != i){
			P2[i,j] = sum(AA2[,i]*AA2[,j]*WW2[,j]/fAvec2)/dval2        ## dval ##
			PMLE2[i,j] = sum(AA2[,i]*AA2[,j]*WW2[,j]/fMLE2)/dval2
		}
	}
	P2[i,i] = 1-sum(P2[i,])
	PMLE2[i,i] = 1-sum(PMLE2[i,])
}
##-----------------------------------------------------------------    ## solve theta and pi ##
tmp.P2 = t(t(P2)-diag(n))%*%(t(P2)-diag(n))
tmp.svd2 = svd(tmp.P2)
pihat2 = abs(tmp.svd2$v[,n])
thetahat2 = log(pihat2)-mean(log(pihat2))
##Theta.matrix[,ww] = thetahat
##sum((thetahat-thetavec)^2)
##-----------------------------------------------------------------    ## solve oracle MLE ##
tmp.PMLE2 = t(t(PMLE2)-diag(n))%*%(t(PMLE2)-diag(n))
tmp.svd.MLE2 = svd(tmp.PMLE2)
pihatMLE2 = abs(tmp.svd.MLE2$v[,n])
thetahatMLE2 = log(pihatMLE2)-mean(log(pihatMLE2))
##Theta.oracle[,ww] = thetahatMLE
##sum((thetahatMLE-thetavec)^2)
##-----------------------------------------------------------------    ## compute sample xihat ##
##-----------------------------------------------------------------    
##-----------------------------------------------------------------    
##-----------------------------------------------------------------    
##-----------------------------------------------------------------    ####Uncomment Line 146-149 to get the results for the oracle estimator
##pihat = pihatMLE
##thetahat = thetahatMLE
##P = PMLE
##fAvec = fMLE
##-----------------------------------------------------------------    
##-----------------------------------------------------------------    
##-----------------------------------------------------------------    
##-----------------------------------------------------------------    
##-----------------------------------------------------------------    
Vmatrix = matrix(0,L,n)
tauhatvec = numeric(n)
tmp.pimatrix = t(AA)*pihat
tmp.pivec = colSums(tmp.pimatrix)
tmp.var = numeric(n)
##tmp.V1 = colSums(t(AA)*pihat)
##-----------------------------------------------------------------    
Vmatrix2 = matrix(0,LL,n)
tauhatvec2 = numeric(n)
tmp.pimatrix2 = t(AA2)*pihat2
tmp.pivec2 = colSums(tmp.pimatrix2)
tmp.var2 = numeric(n)
##tmp.V12 = colSums(t(AA2)*pihat2)
##-----------------------------------------------------------------    
for(oo in 1:n){
	##xihatvec[oo] = sum(P[,oo]*pihat-P[oo,]*pihat[oo])
	tauhatvec[oo] = sum(AA[,oo]*(1-pihat[oo]/tmp.pivec)*pihat[oo]/fAvec)/dval
	tauhatvec2[oo] = sum(AA2[,oo]*(1-pihat2[oo]/tmp.pivec2)*pihat2[oo]/fAvec2)/dval2
	##tmp.var[oo] = sum(AA[,oo]*(tmp.pivec-pihat[oo])*pihat[oo]/fAvec/fAvec)*pihat[oo]/dval/dval/tauhatvec[oo]/tauhatvec[oo]
	tmp.var[oo] = sum(AA[,oo]*(tmp.pivec-pihat[oo])/fAvec/fAvec)*pihat[oo]/dval/dval/tauhatvec[oo]/tauhatvec[oo]
	tmp.var2[oo] = sum(AA2[,oo]*(tmp.pivec2-pihat2[oo])/fAvec2/fAvec2)*pihat2[oo]/dval2/dval2/tauhatvec2[oo]/tauhatvec2[oo]
	Vmatrix[,oo] = (AA[,oo]*WW[,oo]*tmp.pivec-AA[,oo]*pihat[oo])/fAvec  
	Vmatrix2[,oo] = (AA2[,oo]*WW2[,oo]*tmp.pivec2-AA2[,oo]*pihat2[oo])/fAvec2  
}
sigmahatmatrix = matrix(tmp.var,n,n)+t(matrix(tmp.var,n,n))
sigmahatmatrix2 = matrix(tmp.var2,n,n)+t(matrix(tmp.var2,n,n))
##-----------------------------------------------------------------    ## Bootstrap over setM ##
Wmatrix = matrix(rnorm(L*B),L,B)
tmp.Vtau = (t(Vmatrix)/tauhatvec)%*%Wmatrix
Mval = length(setM)
GMvecmax = numeric(B)-1
GMvecmaxone = numeric(B)-Inf
tmpTMval = -1
tmpTMvalone = -Inf
##-----------------------------------------------------------------    
Wmatrix2 = matrix(rnorm(LL*B),LL,B)
tmp.Vtau2 = (t(Vmatrix2)/tauhatvec2)%*%Wmatrix2
Mval2 = length(setM)
GMvecmax2 = numeric(B)-1
GMvecmaxone2 = numeric(B)-Inf
tmpTMval2 = -1
tmpTMvalone2 = -Inf
##-----------------------------------------------------------------    ##
##thetadiff = thetahat-thetavec
for(ooo in setM){
	tmpGMmatrix0 = matrix(rep(tmp.Vtau[ooo,],n)-c(t(tmp.Vtau)),B,n)
	tmpGMmatrix = abs(t(t(tmpGMmatrix0)/sqrt(sigmahatmatrix[ooo,]))/dval)
	tmpGMmatrixone = t(t(tmpGMmatrix0)/sqrt(sigmahatmatrix[ooo,]))/dval
	tmp.GMvecmax = apply(tmpGMmatrix,1,max)
	tmp.GMvecmaxone = apply(tmpGMmatrixone,1,max)
	GMvecmax = c(GMvecmax,tmp.GMvecmax)
	GMvecmaxone = c(GMvecmaxone,tmp.GMvecmaxone)
	##tmp.TMstat = max(abs((thetadiff[ooo]-thetadiff)/sqrt(sigmahatmatrix[ooo,])))
	##tmp.TMstatone = max((thetadiff[ooo]-thetadiff[-ooo])/sqrt(sigmahatmatrix[ooo,-ooo]))
	##tmpTMval = max(tmpTMval,tmp.TMstat) 
	##tmpTMvalone = max(tmpTMvalone,tmp.TMstatone) 
##-----------------------------------------------------------------    ##
	tmpGMmatrix02 = matrix(rep(tmp.Vtau2[ooo,],n)-c(t(tmp.Vtau2)),B,n)
	tmpGMmatrix2 = abs(t(t(tmpGMmatrix02)/sqrt(sigmahatmatrix2[ooo,]))/dval2)
	tmpGMmatrixone2 = t(t(tmpGMmatrix02)/sqrt(sigmahatmatrix2[ooo,]))/dval2
	tmp.GMvecmax2 = apply(tmpGMmatrix2,1,max)
	tmp.GMvecmaxone2 = apply(tmpGMmatrixone2,1,max)
	GMvecmax2 = c(GMvecmax2,tmp.GMvecmax2)
	GMvecmaxone2 = c(GMvecmaxone2,tmp.GMvecmaxone2)
	##tmp.TMstat2 = max(abs((thetadiff[ooo]-thetadiff)/sqrt(sigmahatmatrix[ooo,])))
	##tmp.TMstatone = max((thetadiff[ooo]-thetadiff[-ooo])/sqrt(sigmahatmatrix[ooo,-ooo]))
	##tmpTMval = max(tmpTMval,tmp.TMstat) 
	##tmpTMvalone = max(tmpTMvalone,tmp.TMstatone) 
}
GMmaxmatrix = matrix(GMvecmax,B)
GMmax = apply(GMmaxmatrix,1,max)
cutval = quantile(GMmax,0.975)
##GAthetatwo[ww] = 1*(tmpTMval>cutval)
##-----------------------------------------------------------------    ##
GMmaxmatrixone = matrix(GMvecmaxone,B)
GMmaxone = apply(GMmaxmatrixone,1,max)
cutvalone = quantile(GMmaxone,0.975)
##GAthetaone[ww] = 1*(tmpTMvalone>cutvalone)
##-----------------------------------------------------------------    ##
GMmaxmatrix2 = matrix(GMvecmax2,B)
GMmax2 = apply(GMmaxmatrix2,1,max)
cutval2 = quantile(GMmax2,0.975)
##GAthetatwo[ww] = 1*(tmpTMval>cutval)
##-----------------------------------------------------------------    ##
GMmaxmatrixone2 = matrix(GMvecmaxone2,B)
GMmaxone2 = apply(GMmaxmatrixone2,1,max)
cutvalone2 = quantile(GMmaxone2,0.975)
##GAthetaone[ww] = 1*(tmpTMvalone>cutvalone)
##-----------------------------------------------------------------    ##
##tmp.GAranktwo = 0
##tmp.GArankone = 0
##tmp.oo = 0
for(oooo in setM){
	##tmp.oo = tmp.oo+1
	##print(tmp.oo)
	tmp.theta.sd = sqrt(sigmahatmatrix[oooo,])
	tmp.theta.sd = tmp.theta.sd[-oooo]
	R.left = 1+sum(1*(((thetahat[-oooo]-thetahat[oooo])/tmp.theta.sd)>cutval))
	R.right = n-sum(1*(((thetahat[-oooo]-thetahat[oooo])/tmp.theta.sd)<(-cutval)))
	##CIlengthtwo[tmp.oo,ww] = R.right-R.left
	##tmp.GArankvaluetwo = sum(1*(oooo<R.left))+sum(1*(oooo>R.right))
	##tmp.GAranktwo = max(tmp.GAranktwo,tmp.GArankvaluetwo)
	##R.left.one = 1+sum(1*(((thetahat[-oooo]-thetahat[oooo])/tmp.theta.sd)>cutvalone))
	##Rankleftone[tmp.oo,ww] = R.left.one
	##tmp.GArankvalueone = sum(1*(oooo<R.left.one))
	##tmp.GArankone = max(tmp.GArankone,tmp.GArankvalueone)
	tmp.theta.sd2 = sqrt(sigmahatmatrix2[oooo,])
	tmp.theta.sd2 = tmp.theta.sd2[-oooo]
	R.left2 = 1+sum(1*(((thetahat2[-oooo]-thetahat2[oooo])/tmp.theta.sd2)>cutval2))
	R.right2 = n-sum(1*(((thetahat2[-oooo]-thetahat2[oooo])/tmp.theta.sd2)<(-cutval2)))
}
statone[ww] = 1*((max(R.right,R.right2)-min(R.left,R.left2))>(R.right+R.right2-R.left-R.left2))
write.csv(statone,"statone.csv")
##GAranktwo[ww] = tmp.GAranktwo
##GArankone[ww] = tmp.GArankone
##statone[ww] = (thetahat[10]-thetavec[10])/sqrt(tmp.var[10])
}
mean(statone) ###Power output in EC.4


