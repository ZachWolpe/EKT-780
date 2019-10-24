options ls=132 ps=1000 nodate pageno=1 nocenter;

libname ekt "c:\departement\ekt720\knearest" ;

proc iml ;

use ekt.yx ;
read all into xy ;
*xy=xy[1:15,] ;

start knn ;
	y=xy[,1] ;
	x=xy[,2:3] ;

	* print x ;

	dist= (distance(x))[1:ntr,] ;

*     print dist ;

    res=J(n,1,.) ;
	do i = 1 to n ;
		 dist_c = dist[,i] || y[1:ntr,] ; 
*		 print i dist_c ;

		 call sort(dist_c,{1}) ;

		 dist_c = dist_c[1:k,] ;

*		 print i dist_c ;

		 cnt= (dist_c[,2]=1)[+] || (dist_c[,2]=2)[+] || (dist_c[,2]=3)[+] ;
*		 print cnt ;
		 class = cnt[<:>] ;
		 res[i] = class ;
	 end ;

*	print res ;

	mod = xy || res ;
	nm={"y" "x1" "X2" "gr" "yh"} ;

*	print mod[colname=nm] ;

	pactr = ((mod[,1]=mod[,5])[1:ntr,])[+] / ntr ;
	pacte = ((mod[,1]=mod[,5])[ntr+1:n])[+]/ nte ;

*	print k n ntr nte pactr pacte ;
finish knn ;

k=50 ;
kfold=10;
prop= (kfold-1)/kfold ;

n=nrow(xy) ; 
ntr=round(n*prop) ;
nte=n-ntr ;

print ntr nte ;
xys = ranuni(J(n,1,0)) || xy ;
call sort(xys,{1}) ;
xy =  xys[,2:ncol(xys)] ;
xystart=xy ;
*print xy ;

sobs = (1:n) ;

do ii = 1 to kfold ;

*print "**************************************" ;
*print "Start kfold:" ii ; 

nst = (ii-1)*nte+1 ;
nen = ii*nte ;
nsel = (nst:nen) ;
nselid = J(1,n,0) ;
nselid[nsel]=1 ;
seltest = loc(nselid=1) ;
seltrain = loc(nselid=0) ;

xy = (xystart[seltrain,] || J(ntr,1,1) ) //  
     (xystart[seltest,]  || J(nte,1,2) );

*print nst nen nsel nselid seltest seltrain ; *xystart xy ;
call knn ;

res_kfold = res_kfold //
     (ii || k ||pactr || pacte) ;
end ;

nm={"fold" "k" "pactr" "pacte"} ;
print res_kfold[colname=nm] ;

ares_kfold = res_kfold[:,2:ncol(res_kfold)] ;
nm1={"k" "pactr" "pacte"} ;
print ares_kfold[colname=nm1] ;

quit ;


