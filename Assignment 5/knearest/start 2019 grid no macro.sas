
quit ;
libname ekt "c:\departement\ekt720\knearest" ;
quit ;
symbol1 value=dot 
        height=2 i=none;
symbol2 value=dot 
        height=2 i=none;
symbol3 value=dot 
        height=2 i=none;

proc gplot data=ekt.yx ;
plot x2*x1=group ;
title;
run ;

proc freq data=ekt.yx ;
tables group;
run ;

proc means data=ekt.yx ; ;
run ;


proc iml ;

*own distance function - memory problems ;
start ddist ;
dist = J(n,nrow(xyg),.) ;
do i = 1 to  nrow(xyg) ;
  do j = 1 to n  ;
   dist[j,i] = sqrt( (xyg[i,2] - xyg[j,2])**2 + (xyg[i,3] - xyg[j,3])**2 ) ;
  end;
end;
finish ddist ;

use ekt.yx  ;
read all  into xy ;
xy=xy[1:50,] ;

n=nrow(xy) ; 
print xy ;

y=xy[,1] ;
x=xy[,2:3] ;

print x ;

 k=5 ;

 r=10 ;
 step=10;
 x1min=round(min(xy[,2])) - r ;
 x1max=round(max(xy[,2])) + r ;
 x2min=round(min(xy[,3])) - r ;
 x2max=round(max(xy[,3])) + r ;

 print x1min x1max x2min x2max ;

 do i = x1min to x1max by step ;
    do j = x2min to x2max by step ;
	 xyg = xyg // (555  ||   i || j) ;
    end ; 
 end ;

xyg = (xy || J(n,1,1)) //  
      (xyg || J(nrow(xyg),1,2)) ;
print xyg ,(nrow(xyg)) ;

*dist= distance(xyg[,2:3])[1:n,] ;

call ddist ;

print dist ;

 do i = 1 to ncol(dist) ;
	 dist_c = dist[,i] || y ; 
if i <=5 then	 print i dist_c ;

	 call sort(dist_c,{1}) ;

	 dist_c = dist_c[1:k,] ;

if i <=5 then 		 print "Sorted " dist_c ;

	 cnt= (dist_c[,2]=1)[+] || (dist_c[,2]=2)[+] 
       || (dist_c[,2]=3)[+] ;
if i <=5 then 		 print cnt ;
	 class = cnt[<:>] ;
	 res = res // class ;
if i <=5 then 		 print class ;
 end ;

print res ;

mod = xyg || res ;
nm={"y" "x1" "X2" "gr" "yh"} ;

print mod[colname=nm] ;

pac = (mod[1:n,1]=mod[1:n,5])[+]/n ;
print pac ;

create kdat from mod[colname=nm] ;
append from mod ;
close kdat ;


quit ;

proc gplot data=kdat ;
plot x1*x2=yh ;
run ;
quit ;

data kdat ;
 set kdat ;
 yh1=yh ;
 if gr=1 then yh1=y+10;
run ;


symbol1 interpol=none width=3
 		color=blue
        value=dot
        height=1;

symbol2 interpol=none width=3
 		color=red
        value=dot
        height=1;

symbol3 interpol=none width=3
 		color=green
        value=dot
        height=1;

symbol4 interpol=none width=3
 		color=black
        value=dot
        height=3;

symbol5 interpol=none width=3
 		color=orange
        value=dot
        height=3;

symbol6 interpol=none width=3
 		color=purple
        value=dot
        height=3;

proc gplot data=kdat ;
plot x1*x2=yh1 ;
title "KNN " ;
run ;

proc gplot data=kdat ;
plot x1*x2=yh1 ;
where yh1 <=3 ;
title "KNN GRID " ;
run ;
quit ;

title ;

quit ;

