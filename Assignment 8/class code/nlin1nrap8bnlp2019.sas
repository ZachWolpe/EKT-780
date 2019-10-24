quit ;
dm 'odsresults;clear';


options ls=72 nodate pageno=1 ;

data a (keep= x y);
s=0.5;
do i=1 to 100;
x=i;
y=sqrt((2.35*x+5.89)) + rannor(5675)*s ;
*y2=y*y;
output;
end;
run;

title ;
title2 ;

proc print data=a ;
run;

proc nlin data=a  ; *list  method=newton;
 parms a=100, b=100 ;
 model y = sqrt(a*X+b) ;
run ;




proc iml ;
use a;
read all into xy ;

n = nrow(xy) ;
y = xy[,2] ;
x = xy[,1] ;

start sse(parms) global (x,y);
   bh=parms ;
   yh= sqrt(bh[1]*x+bh[2]) ;
   sse = (y-yh)`*(y-yh) ; 
	/*   print sse ;*/
   return (sse);
finish;

parms = {2,6} ;
sse= sse(parms) ;

print "Initial values " parms sse ;

opt = {0,    /* find maximum of function   */
       4};   /* print a LOT of output      */
      
call nlpnra(rc, parm_estimates, "sse", parms, opt);

print parm_estimates ;
quit ;


* question 6








