
#eigs = eigenvalues;   n.ind = number of individuals; k = number of axes to test
TWcalc<-function(eigs, n.inds, k){
	#make sure eigenvalues are sorted
	eigs = sort(eigs, decr=T)
	#load a table for test statistics
	test<-read.table("twtable.txt", header=FALSE)
	m0=n.inds
	k=k
	# calculate the statistic for dimensions m0
	# equations are from Patterson et al. 2006
	TWres_<-c()
	for (i in 1:k){
		m=(m0)-i+1
		m=m
		#use the eigenvalues from k to mprime (step 3 of "a test for population structure")
		eiv<-eigs[i:(m-1)]
		#nprime is equation 10
		nprime=((m+1)*sum(eiv)^2)/((m-1)*sum(eiv^2)-sum(eiv)^2)
		#S is equation 6
		S=((sqrt(nprime-1)+sqrt(m)))/nprime*(1/sqrt(nprime-1)+1/sqrt(m))^(1/3)
		#u_ is equation 5
		u_=(sqrt(nprime-1)+sqrt(m))^2/nprime
		#L1 is step 5 in "a test for population structure"
		L1<-(m-1)*eiv[1]/sum(eiv)
		#the test statistic: equation 7
		TW_<-(L1-u_)/S
		# collate test statistics
		TWres_<-c(TWres_,TW_)
	}
	
	#calculate p values
	if (require(RMTstat)){
		pvals = c()
		for (i in 1:length(TWres_)){
			p = ptw(TWres[i], lower.tail=F)
			pvals = c(pvals, p)
		}
	}
	else{
		pvals = NA
	}
	pres<-c()
	for (i in 1:length(TWres_)){
		TW<-TWres_[i]
		dif<-abs(test[,1]-TW)
		p<-test[dif==min(dif),2]
		pres<-c(pres,p)
	}
	
	res<-list(TWstats=TWres_, Sig=pres, pval=pvals)
	return(res)
}


