###Instructions
#AA Matrix of the Graph
#WW Matrix of Winner Indices:Encodes observed winner information.
#Note: Uncomment lines 123–126 to obtain the oracle versions of the estimators for all settings.
#Output Reference Guide
#•	Figure EC.1: Results are generated in lines 201–202.
#•	Figure EC.2: Results are generated in lines 219–221.
#•	Figure EC.3: Coverage levels are set in lines 169 and 174. You may modify the nominal level (e.g., change 0.95 to 0.85, 0.75, etc.). Output is saved in GAthetatwo.csv at line 206.
#•	Table EC.1: Outputs are generated at lines 206 and 208.
#•	Table EC.2: Set setM = c(3:10) in line 22. Outputs are generated at lines 210 and 212.
#•	Table EC.3: Set setM = 1:n in line 23. Outputs are again at lines 210 and 212.
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ##       ##  INPUT  ##
##-----------------------------------------------------------------    ##
NN = 500 #repeat 500 times
B = 300  #number of bootstrap samples
L = 36000 #|D| or L number of comparisons, can be changed under different settings
n = 50 #total number of items
n1 = n*0.2; n2 = n*0.5 
L1 = round(L/5); L2 = round(L*2/5)
thetavec = seq(6,2,length=n)  #\theta
thetavec = thetavec-mean(thetavec) #normalized theta
setM = c(8,20,30) #entries we are interested in
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ##       ##  OUTPUT  ##
##-----------------------------------------------------------------    ##
#Theta.matrix = matrix(0,n,NN)
#Theta.oracle = matrix(0,n,NN)
#GAthetatwo = numeric(NN)
#GAranktwo = numeric(NN)
#GAthetaone = numeric(NN)
#GArankone = numeric(NN)
#CIlengthtwo = matrix(0,length(setM),NN)
#Rankleftone = matrix(0,length(setM),NN)
#statone = numeric(NN)
#statone1=numeric(NN)
#statone2=numeric(NN)
Theta.matrix = matrix(0,n,NN)
#Theta_vector<-c()
#Theta_vector_o<-c()
Theta.oracle = matrix(0,n,NN)
GAthetatwo = numeric(NN)
GAranktwo = numeric(NN)
GAthetaone = numeric(NN)
GArankone = numeric(NN)
CIlengthtwo = matrix(0,length(setM),NN)
Rankleftone = matrix(0,length(setM),NN)
statone = numeric(NN)
statone1=numeric(NN)
statone2=numeric(NN)
##-----------------------------------------------------------------    ##
for(ww in 1:NN){
print(ww)
fAvec = numeric(L)
fMLE = numeric(L)
AA = matrix(0,L,n)
WW = matrix(0,L,n)
for(o in 1:L1){
	tmp.M = sample(c(2,3,4,5),1)                                     ## 
	##tmp.M = 4
	fAvec[o] = tmp.M
	tmp.idx = sample(1:n1,tmp.M)
	AA[o,tmp.idx] = 1
	tmp.w = exp(thetavec[tmp.idx])
	WW[o,tmp.idx] = rmultinom(1,1,tmp.w) 
	fMLE[o] = sum(tmp.w)
}
for(o in (L1+1):L2){
	tmp.M = sample(c(2,3,4,5),1)                                     ## 
	##tmp.M = 4
	fAvec[o] = tmp.M
	tmp.idx = sample(1:n2,tmp.M)
	AA[o,tmp.idx] = 1
	tmp.w = exp(thetavec[tmp.idx])
	WW[o,tmp.idx] = rmultinom(1,1,tmp.w) 
	fMLE[o] = sum(tmp.w)
}
for(o in (L2+1):L){
	tmp.M = sample(c(2,3,4,5),1)                                     ## 
	##tmp.M = 4
	fAvec[o] = tmp.M
	tmp.idx = sample(1:n,tmp.M)
	AA[o,tmp.idx] = 1
	tmp.w = exp(thetavec[tmp.idx])
	WW[o,tmp.idx] = rmultinom(1,1,tmp.w) 
	fMLE[o] = sum(tmp.w)
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
Theta.matrix[,ww] = thetahat
##sum((thetahat-thetavec)^2)
##-----------------------------------------------------------------    ## solve oracle MLE ##
tmp.PMLE = t(t(PMLE)-diag(n))%*%(t(PMLE)-diag(n))
tmp.svd.MLE = svd(tmp.PMLE)
pihatMLE = abs(tmp.svd.MLE$v[,n])
thetahatMLE = log(pihatMLE)-mean(log(pihatMLE))
Theta.oracle[,ww] = thetahatMLE
##sum((thetahatMLE-thetavec)^2)
##-----------------------------------------------------------------    ## compute sample xihat ##
##-----------------------------------------------------------------    
##-----------------------------------------------------------------    
##-----------------------------------------------------------------    
##-----------------------------------------------------------------    
#pihat = pihatMLE                                                      ###Uncomment this to get the oracle estimator.
#thetahat = thetahatMLE
#P = PMLE
#fAvec = fMLE
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
tmp.V1 = colSums(t(AA)*pihat)
for(oo in 1:n){
	##xihatvec[oo] = sum(P[,oo]*pihat-P[oo,]*pihat[oo])
	tauhatvec[oo] = sum(AA[,oo]*(1-pihat[oo]/tmp.pivec)*pihat[oo]/fAvec)/dval
	##tmp.var[oo] = sum(AA[,oo]*(tmp.pivec-pihat[oo])*pihat[oo]/fAvec/fAvec)*pihat[oo]/dval/dval/tauhatvec[oo]/tauhatvec[oo]
	tmp.var[oo] = sum(AA[,oo]*(tmp.pivec-pihat[oo])/fAvec/fAvec)*pihat[oo]/dval/dval/tauhatvec[oo]/tauhatvec[oo]
	Vmatrix[,oo] = (AA[,oo]*WW[,oo]*tmp.pivec-AA[,oo]*pihat[oo])/fAvec  
}
sigmahatmatrix = matrix(tmp.var,n,n)+t(matrix(tmp.var,n,n))
##-----------------------------------------------------------------    ## Bootstrap over setM ##
Wmatrix = matrix(rnorm(L*B),L,B)
tmp.Vtau = (t(Vmatrix)/tauhatvec)%*%Wmatrix
Mval = length(setM)
GMvecmax = numeric(B)-1
GMvecmaxone = numeric(B)-Inf
tmpTMval = -1
tmpTMvalone = -Inf
##-----------------------------------------------------------------    ##
thetadiff = thetahat-thetavec
for(ooo in setM){
	tmpGMmatrix0 = matrix(rep(tmp.Vtau[ooo,],n)-c(t(tmp.Vtau)),B,n)
	tmpGMmatrix = abs(t(t(tmpGMmatrix0)/sqrt(sigmahatmatrix[ooo,]))/dval)
	tmpGMmatrixone = t(t(tmpGMmatrix0)/sqrt(sigmahatmatrix[ooo,]))/dval
	tmp.GMvecmax = apply(tmpGMmatrix,1,max)
	tmp.GMvecmaxone = apply(tmpGMmatrixone,1,max)
	GMvecmax = c(GMvecmax,tmp.GMvecmax)
	GMvecmaxone = c(GMvecmaxone,tmp.GMvecmaxone)
	tmp.TMstat = max(abs((thetadiff[ooo]-thetadiff)/sqrt(sigmahatmatrix[ooo,])))
	tmp.TMstatone = max((thetadiff[ooo]-thetadiff[-ooo])/sqrt(sigmahatmatrix[ooo,-ooo]))
	tmpTMval = max(tmpTMval,tmp.TMstat) 
	tmpTMvalone = max(tmpTMvalone,tmp.TMstatone) 
}
GMmaxmatrix = matrix(GMvecmax,B)
GMmax = apply(GMmaxmatrix,1,max)
cutval = quantile(GMmax,0.95)
GAthetatwo[ww] = 1*(tmpTMval>cutval)
##-----------------------------------------------------------------    ##
GMmaxmatrixone = matrix(GMvecmaxone,B)
GMmaxone = apply(GMmaxmatrixone,1,max)
cutvalone = quantile(GMmaxone,0.95)
GAthetaone[ww] = 1*(tmpTMvalone>cutvalone)
##-----------------------------------------------------------------    ##
tmp.GAranktwo = 0
tmp.GArankone = 0
tmp.oo = 0
for(oooo in setM){
	tmp.oo = tmp.oo+1
	print(tmp.oo)
	tmp.theta.sd = sqrt(sigmahatmatrix[oooo,])
	tmp.theta.sd = tmp.theta.sd[-oooo]
	R.left = 1+sum(1*(((thetahat[-oooo]-thetahat[oooo])/tmp.theta.sd)>cutval))
	R.right = n-sum(1*(((thetahat[-oooo]-thetahat[oooo])/tmp.theta.sd)<(-cutval)))
	CIlengthtwo[tmp.oo,ww] = R.right-R.left
	tmp.GArankvaluetwo = sum(1*(oooo<R.left))+sum(1*(oooo>R.right))
	tmp.GAranktwo = max(tmp.GAranktwo,tmp.GArankvaluetwo)
	R.left.one = 1+sum(1*(((thetahat[-oooo]-thetahat[oooo])/tmp.theta.sd)>cutvalone))
	Rankleftone[tmp.oo,ww] = R.left.one
	tmp.GArankvalueone = sum(1*(oooo<R.left.one))
	tmp.GArankone = max(tmp.GArankone,tmp.GArankvalueone)
}
GAranktwo[ww] = tmp.GAranktwo
GArankone[ww] = tmp.GArankone
#statone[ww] = (thetahat[10]-thetavec[10])/sqrt(tmp.var[10])
statone[ww] = (thetahat[8]-thetavec[8])/sqrt(tmp.var[8])
statone1[ww] = (thetahat[20]-thetavec[20])/sqrt(tmp.var[20])
statone2[ww] = (thetahat[30]-thetavec[30])/sqrt(tmp.var[30])
write.csv(Theta.matrix,"Theta_matrix.csv")  #
write.csv(Theta.oracle,"Theta_oracle.csv")
#write.csv(Theta_vector,"Theta_vector.csv")
#write.csv(Theta_vector_o,"Theta_vector_o.csv")
#Theta.oracle = matrix(0,n,NN)
write.csv(GAthetatwo,"GAthetatwo.csv")   #coverage for theta in EC.1 and also pp-plot
#GAthetatwo = numeric(NN)
write.csv(GAranktwo,"GAranktwo.csv")  #EC(r) in Table EC.1
#GAranktwo = numeric(NN)
write.csv(GAthetaone,"GAthetaone.csv") #coverage for one-sided confidence intervals
#GAthetaone = numeric(NN)
write.csv(GArankone,"GArankone.csv") #coverage one-sided confidence intervals. used in EC2
#GArankone = numeric(NN)
write.csv(CIlengthtwo,"CIlengthtwo.csv") # #Length of two-sided confidence intervals in Table EC.1
#CIlengthtwo = matrix(0,length(setM),NN)
write.csv(Rankleftone,"Rankleftone.csv") #left boundary of one-sided confidence intervals. Used in Table EC.2
#Rankleftone = matrix(0,length(setM),NN)
#statone = numeric(NN)
write.csv(statone,"statone.csv")   #asymptotic normality entries 8,20,30 on Figure EC.2.
write.csv(statone1,"statone1.csv") 
write.csv(statone2,"statone2.csv")
}







