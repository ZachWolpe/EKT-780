quit ;
dm 'odsresults;clear';


options ls=72 nodate pageno=1 nocenter ; 
*ods pdf file='C:\departement\ekt720\nonlin\n1.pdf';


title1 'Non linear estimation, Newton Raphson'  ;

* data generation for illustrative purposes ;
data a ;
  do i = 1 to 20 ;
   x=i ;
   y = 30 - 80*x/(exp(x)) +rannor(11567)*1 ;
  output ;
  end ;
run;



*print of the data ;
 proc print data =a ;
  run ;

*plot of data ;

proc gplot data=a ;
 plot y*x ;
run ; 


*estimation of paramters using proc nlin - not part of the assignment - test results ;

proc nlin method=newton;
parms th1 = 1 th2 =1 ;
model y = th1 + th2*x/(exp(x)) ;
run ;

proc iml ;

use a ;
 read all into yx ;

 print yx ;

 y=yx[,3] ;
 x=yx[,2] ;

 n=nrow(yx) ;

 print y x n ;

 theta={ 10 , 1000} ;

mdiff=10000000 ;

 do i = 1 to 100  untill (diff<0.000005);

 d1 = theta[2,1]*x ;
 d2 = exp(-x) ;
 
*print x d1 d2 ;

 yh = J(n,1,theta[1,1])+(d1#d2) ;

* print i yh ;

 print i theta ;

*derivatives ;

*first ;
 q1 = -2*sum(y-yh) ;
 q2 = -2*sum((y-yh)#(d2#x)) ;

 q = q1 || q2 ;

*hessian ;

 h11 = 2*n;
 h12 = 2*sum(x#d2) ;
 h21 = h12 ;
 h22 = 2*sum((x#d2)#(x#d2)) ;

 h = (h11 || h12) // (h21 || h22) ;

 print i q h;

 thetanew = theta - inv(h)*q` ;

 print i thetanew ;


 diff = max(abs(theta-thetanew)) ;
 
 if diff < mdiff then mdiff = diff ;

 theta=thetanew ;
    
 print diff ;

end ;

 d1 = theta[2,1]*x ;
 *d1=x ;
 d2 = exp(-x) ;

 yh = J(n,1,theta[1,1])+(d1#d2) ;

print yh ; 

rr = y || yh || x ;

nm1={"y" "yh" "x"};
print rr[colname=nm1];


create yd from rr[colname=nm1] ;
append from rr ;


 quit ;

symbol1 colour=red value=dot height = 1;
symbol2 colour=green  i=join width=2;

proc gplot data=yd;
plot (y yh)*x /overlay;
run;

quit; 

proc sgplot data=yd;
	scatter y=y x=x;
	series y=yh x=x;
run;


*ods pdf close;
