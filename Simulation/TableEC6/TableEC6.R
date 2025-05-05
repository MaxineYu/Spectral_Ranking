gradient.f = function(w,A,W1,W2,W3){
	wexp = exp(w)
	tmp.sum1 = colSums(t(A)*wexp)
	tmp.sum2 = colSums(t(A-W1)*wexp)
	f.theta = exp(colSums(t(W1+W2)*w))/(tmp.sum1)/(tmp.sum2)
	tmp.gradient1 = exp(colSums(t(W1+W2)*w))/(tmp.sum1)/(tmp.sum1)/f.theta
	tmp.gradient1 = W1*tmp.gradient1
	tmp.gradient2 = exp(colSums(t(W1+W2)*w))*(tmp.sum1*tmp.sum2-exp(colSums(t(W2)*w))*tmp.sum2-exp(colSums(t(W2)*w))*tmp.sum1)/(tmp.sum1)/(tmp.sum1)/(tmp.sum2)/(tmp.sum2)/f.theta
	tmp.gradient2 = W2*tmp.gradient2
	tmp.gradient3 = -exp(colSums(t(A)*w))*(tmp.sum1+tmp.sum2)/(tmp.sum1)/(tmp.sum1)/(tmp.sum2)/(tmp.sum2)/f.theta
	tmp.gradient3 = W3*tmp.gradient3
	return(colSums(tmp.gradient1+tmp.gradient2+tmp.gradient3))
}
MLE = function(AA,WW1,WW2,WW3){
	wvec0 = rnorm(n); wvec0 = wvec0-mean(wvec0)	
	tmp.Gradient = gradient.f(wvec0,AA,WW1,WW2,WW3)
	G.val = max(abs(tmp.Gradient))
	tmp.wvec = wvec0
	Kval = 1
	while(G.val>0.0000001 && Kval <= 1000){
		print(G.val)
		tmp.wvec = tmp.wvec+0.001*tmp.Gradient
		tmp.Gradient = gradient.f(tmp.wvec,AA,WW1,WW2,WW3)
		G.val = max(abs(tmp.Gradient))
		Kval = Kval+1
	}
	return(tmp.wvec)
}
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ##              ##  INPUT  ##
##-----------------------------------------------------------------    ##
p = 0.05  #choose between $0.02, 0.05,0.08, 0.11,0.14$.
n = 50
L = 10
NN = 100          ## number of independent replications
thetavec = seq(6,2,length=n)                                    
thetavec = thetavec-mean(thetavec)      ## centered true theta vector 
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ##              ##  OUTPUT  ##
##-----------------------------------------------------------------    ##
L2matrix = matrix(0,4,NN)             ## L2 error for four methods 
Linftymatrix = matrix(0,4,NN)         ## Linfty error for four methods 
theta.spectral = matrix(0,n,NN)       ## collection of theta hat for spectral method 
theta.twostep = matrix(0,n,NN)        ## collection of theta hat for two-step spectral method 
theta.oracle = matrix(0,n,NN)         ## collection of theta hat for spectral method with oracle MLE weight
theta.MLE = matrix(0,n,NN)            ## collection of theta hat for MLE
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ##
##-----------------------------------------------------------------    ##
for(ww in 1:NN){
print(ww)
AA = matrix(0,0,n)
WW1 = matrix(0,0,n)
WW2 = matrix(0,0,n)
for(i in 1:(n-2)){
	for(j in (i+1):(n-1)){
		for(k in (j+1):n){
			tmp.graph = rbinom(1,1,p)
			if(tmp.graph==1){
			for(l in 1:L){	
				tmp.AA = numeric(n)
				tmp.AA[c(i,j,k)] = 1
				AA = rbind(AA,tmp.AA)
			tmp.comp1 = rmultinom(1,1,exp(thetavec[c(i,j,k)]))
			if(tmp.comp1[1,1]==1){
				tmp.WW1 = numeric(n)
				tmp.WW1[i] = 1	
				WW1 = rbind(WW1,tmp.WW1)
				tmp.comp2 = rmultinom(1,1,exp(thetavec[c(j,k)]))
				tmp.WW2 = numeric(n)
				tmp.WW2[c(j,k)] = tmp.comp2
				WW2 = rbind(WW2,tmp.WW2)
			}
			if(tmp.comp1[2,1]==1){
				tmp.WW1 = numeric(n)
				tmp.WW1[j] = 1	
				WW1 = rbind(WW1,tmp.WW1)
				tmp.comp2 = rmultinom(1,1,exp(thetavec[c(i,k)]))
				tmp.WW2 = numeric(n)
				tmp.WW2[c(i,k)] = tmp.comp2
				WW2 = rbind(WW2,tmp.WW2)
			}
			if(tmp.comp1[3,1]==1){
				tmp.WW1 = numeric(n)
				tmp.WW1[k] = 1	
				WW1 = rbind(WW1,tmp.WW1)
				tmp.comp2 = rmultinom(1,1,exp(thetavec[c(i,j)]))
				tmp.WW2 = numeric(n)
				tmp.WW2[c(i,j)] = tmp.comp2
				WW2 = rbind(WW2,tmp.WW2)
			}
			}
			}
		}
	}
}
##-----------------------------------------------------------------    ##
dval = 2*max(colSums(AA))
fMLE = rowSums(AA*t(matrix(exp(thetavec),n,nrow(AA))))
P = matrix(0,n,n)
PMLE = matrix(0,n,n)
for(i in 1:n){
	for(j in 1:n){
		if(j != i){
			tmp.AA = AA[,-c(i,j)]
			tmp.Pval = sum(tmp.AA*AA[,i]*AA[,j]*WW1[,j])/3
			tmp.Pval2 = sum(tmp.AA*AA[,i]*AA[,j]*WW1[,-c(i,j)]*WW2[,j])/2
			P[i,j] = (tmp.Pval+tmp.Pval2)/dval                   
			wvec = exp(thetavec)
			tmp.wvec = (wvec[i]+wvec[j])+wvec[-c(i,j)]
			tmp.wMatrix = t(matrix(tmp.wvec,(n-2),nrow(AA)))
			tmp.PvalMLE = sum((tmp.AA/tmp.wMatrix)*AA[,i]*AA[,j]*WW1[,j])
			tmp.Pval2MLE = sum(tmp.AA*AA[,i]*AA[,j]*WW1[,-c(i,j)]*WW2[,j])/(exp(thetavec[i])+exp(thetavec[j]))
			PMLE[i,j] = (tmp.PvalMLE+tmp.Pval2MLE)/dval 
		}
	}
	P[i,i] = 1-sum(P[i,])
	PMLE[i,i] = 1-sum(PMLE[i,])
}
##-----------------------------------------------------------------    ##
tmp.P = t(t(P)-diag(n))%*%(t(P)-diag(n))
tmp.svd = svd(tmp.P)
pihat = abs(tmp.svd$v[,n])
thetahat = log(pihat)-mean(log(pihat))
L2matrix[1,ww] = sqrt(sum((thetahat-thetavec)^2))
Linftymatrix[1,ww] = max(abs(thetahat-thetavec))
theta.spectral[,ww] = thetahat-thetavec
##-----------------------------------------------------------------    ##
tmp.PMLE = t(t(PMLE)-diag(n))%*%(t(PMLE)-diag(n))
tmp.svd.MLE = svd(tmp.PMLE)
pihatMLE = abs(tmp.svd.MLE$v[,n])
thetahatMLE = log(pihatMLE)-mean(log(pihatMLE))
L2matrix[2,ww] = sqrt(sum((thetahatMLE-thetavec)^2))
Linftymatrix[2,ww] = max(abs(thetahatMLE-thetavec))
theta.oracle[,ww] = thetahatMLE-thetavec
##-----------------------------------------------------------------    ##
dval = 2*max(colSums(AA))
ftwo = rowSums(AA*t(matrix(exp(thetahat),n,nrow(AA))))
Ptwo = matrix(0,n,n)
for(i in 1:n){
	for(j in 1:n){
		if(j != i){
			tmp.AA = AA[,-c(i,j)]
			wvec = exp(thetahat)
			tmp.wvec = (wvec[i]+wvec[j])+wvec[-c(i,j)]
			tmp.wMatrix = t(matrix(tmp.wvec,(n-2),nrow(AA)))
			tmp.PvalMLE = sum((tmp.AA/tmp.wMatrix)*AA[,i]*AA[,j]*WW1[,j])
			tmp.Pval2MLE = sum(tmp.AA*AA[,i]*AA[,j]*WW1[,-c(i,j)]*WW2[,j])/(exp(thetahat[i])+exp(thetahat[j]))
			Ptwo[i,j] = (tmp.PvalMLE+tmp.Pval2MLE)/dval 
		}
	}
	Ptwo[i,i] = 1-sum(Ptwo[i,])
}
##-----------------------------------------------------------------    ##
tmp.Ptwo = t(t(Ptwo)-diag(n))%*%(t(Ptwo)-diag(n))
tmp.svd.two = svd(tmp.Ptwo)
pihattwo = abs(tmp.svd.two$v[,n])
thetahattwo = log(pihattwo)-mean(log(pihattwo))
L2matrix[3,ww] = sqrt(sum((thetahattwo-thetavec)^2))
Linftymatrix[3,ww] = max(abs(thetahattwo-thetavec))
theta.twostep[,ww] = thetahattwo-thetavec
##-----------------------------------------------------------------    ##
WW3 = AA-WW1-WW2
thetaMLE = MLE(AA,WW1,WW2,WW3)
L2matrix[4,ww] = sqrt(sum((thetaMLE-thetavec)^2))
Linftymatrix[4,ww] = max(abs(thetaMLE-thetavec))	
theta.MLE[,ww] = thetaMLE-thetavec

write.csv(L2matrix,"L2matrix.csv") #MLE \ell_2
write.csv(Linftymatrix,"Linftymatrix.csv")  #MLE infinity
write.csv(theta.spectral,"theta_spectral.csv") #\theta_vanilla
write.csv(theta.twostep,"theta_twostep.csv") #theta_two-step
write.csv(theta.oracle,"theta_oracle.csv") #\theta_oracle
write.csv(theta.MLE,"theta_MLE.csv") #\theta_{MLE}
}














































































































