libname mid "/folders/myfolders/sasuser.v94/EKT 720/midterm" ;
run ;

data a;
set mid.q2;
run;


symbol1 i=none color=blue ;
symbol2 i=none color=green;

proc means data=a ;
run ;


proc gplot data=a ;
 plot x1*x2;
run ;
quit ;


/* proc gplot data=a ;
 plot x1*x2=y ;
run ; */
proc sgplot data=a;
	scatter x=x1 y=x2 / group=y;
run;

quit ;















proc iml ;
use a ;
read all into xy ;
*xy=xy[1:10,] ;
n=nrow(xy) ;
k=50;

d=distance(xy[,1:2]) ;
res=J(n,2,.) ;

do i = 1 to n ;
 d_k = d[,i] || xy[,3] ;
 call sort(d_k,{1});
 d_k = d_k[2:k,] ;
 w =  1 / d_k[,1];
 class = sum(w#d_k[,2])/sum(w) ; 
* print i d_k w class;
 res[i,] = xy[i,3] ||   round(class) ;
* print res ;
end ;

*print res ;

predacc = (res[,1]=res[,2]) ;
predacc = sum(predacc)/n ;
print  k predacc;

res1 = xy || res || (res[,1]=res[,2]) ;
nm = {"x1" "x2"  "y" "y2" "yh" "corr"} ;

create res_p from res1[colname=nm] ;
append from res1 ;

quit ; 

symbol1 interpol=none width=3
 		color=purple
        value=dot
        height=3;
symbol2 interpol=none width=3
 		color=red
        value=dot
        height=3;
        
proc means data=res_p;
        
proc sgplot data=res_p;
	scatter 

/* proc gplot data=res_p ;
 plot x2*x1=corr ;
run ;
quit ;


