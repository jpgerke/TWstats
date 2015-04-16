Mmake<-function(data){
	pcdat<-data
	mv<-apply(pcdat,2,mean,na.rm=TRUE)
	#make vector of sqrt(p*(1-p))
	vv<-mv/2
	vv<-sqrt(vv*(1-vv))
	#redo data and vectors to take out monomophic sites
	dat<-pcdat[,vv>0&is.na(vv)==FALSE]
	mv<-mv[vv>0&is.na(vv)==FALSE]
	vv<-vv[vv>0&is.na(vv)==FALSE]
	M<-t((t(dat)-mv)/vv)
	M[is.na(M)]<-0
	return(M)	
}
	
TWcalc<-function(dat,k){
	test<-read.table("twtable.txt", header=FALSE)
	x<-dat
	y<-x%*%t(x)
	eig<-eigen(y)
	eig<-eig$values
	eig<-eig[1:(length(eig))]
	n<-ncol(x)
	m0=nrow(x)
	k=k
	TWres_<-c()
	for (i in 1:k){
		m=(m0)-i+1
		m=m
		eiv<-eig[i:(m-1)]
		#n_=n
		n_=((m+1)*sum(eiv)^2)/((m-1)*sum(eiv^2)-sum(eiv)^2)
		S=((sqrt(n_-1)+sqrt(m)))/n_*(1/sqrt(n_-1)+1/sqrt(m))^(1/3)
		u_=(sqrt(n_-1)+sqrt(m))^2/n_
		L1<-(m-1)*eiv[1]/sum(eiv)
		TW_<-(L1-u_)/S
		TWres_<-c(TWres_,TW_)
	}

	#TWres_
	pres<-c()
	for (i in 1:length(TWres_)){
		TW<-TWres_[i]
		dif<-abs(test[,1]-TW)
		p<-test[dif==min(dif),2]
		pres<-c(pres,p)
	}
	#pres
	res<-list(TWres_,pres)
	return(res)
}

#usage, data as 0,1 and 2s
M<-Mmake(dat)
TWcalc(M, 10)
TWcalc(M, 25)

