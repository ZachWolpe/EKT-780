options ls=72 nodate pageno=1 ;

data a (keep=x y) ;
 s=.5 ;
 do i = 1 to 100 ;
   x= i ;
   y = sqrt((2.35*x+5.89))+rannor(5675)*s ;
   output;
 end ;
run ;


title ;
title2 ;

proc gplot data=a ;
 plot y*x ;
run ;

/*proc print data=a ;*/
/*run;*/

proc iml ;
use a;
read all into xy ;
print xy ;

n = nrow(xy) ;
y = xy[,2] ;

results = J(4,1,9999) ;
yh = J(n,1,9999) ;

msebm = 10000000000 ;
*mse =10000000000 ;

step=.01 ;

gg=0 ;

do th1 = 1 to 3 by step;
  do th2 = 2.01 to 7 by step;
    
 do i = 1 to n ;
  yh[i,1] = sqrt(th1*xy[i,1]+th2) ;
 end;

  res = (y-yh) ;
  mse= res`*res/(n-2) ;

 * print  th1 th2  mse ;

 if mse <= msebm then do ;
  results =  (th1 // th2 // mse) ;
  yhf = yh ;
  msebm = mse ;
 end ;

*print mse ;

 gg=gg+ 1 ;
end ;
end;

nm = {"th1" "th2" "MSE"} ;

print gg results[rowname=nm] ;

pdat = xy || yhf ;

nm1= { x y yhf} ;
print pdat[colname=nm1] ;

create pdat from pdat[colname=nm1] ;
append from pdat ;

quit ;

proc gplot data=pdat ;
 plot y*x yhf*x / overlay;
run ;


