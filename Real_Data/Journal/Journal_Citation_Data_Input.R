## Download and load "AuthorPaperInfo.RData" ##
## https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/V7VUFO ## 

xx = matrix(0,0,36)
Idx = unique(AuPapMat$journal)
year = numeric(0)
ww = matrix(0,0,36)
for(o in 1:522062){
	tmp.xx = numeric(36)
	tmp.ww = numeric(36)
	tmp.v1 = AuPapMat$journal[(AuPapMat$idxPap==(PapPapMat$FromPap[o]))]
	tmp.v2 = AuPapMat$journal[(AuPapMat$idxPap==(PapPapMat$ToPap[o]))]
	if(tmp.v1[1] != tmp.v2[1]){
	print(o)
	year = c(year,PapPapMat$FromYear[o])
	tmp.xx[(Idx==tmp.v1[1])] = 1
	tmp.xx[(Idx==tmp.v2[1])] = 1
	tmp.ww[(Idx==tmp.v2[1])] = 1
	xx = rbind(xx,tmp.xx)
	ww = rbind(ww,tmp.ww)
	}
}



tmp.yy = c(xx)
yy = matrix(tmp.yy,ncol=36)
tmp.zz = c(ww)
zz = matrix(tmp.zz,ncol=36)



write.table(yy,file='Ranking_AA')
write.table(Idx,file='Ranking_Journal')
write.table(year,file='Ranking_Year')
write.table(zz,file='Ranking_WW')










































