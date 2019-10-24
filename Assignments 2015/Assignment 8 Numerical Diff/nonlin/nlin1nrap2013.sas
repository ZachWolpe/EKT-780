options ls=72 nodate pageno=1 ;

data a (keep= x y y2);
s=0.5;
do i=1 to 100;
x=i;
y=sqrt((2.35*x+5.89)) + rannor(5675)*s ;
y2=y*y;
output;
end;
run;

title ;
title2 ;

proc print data=a ;
run;

proc nlin data=a ;
 parms a=100, b=100 ;
 model y = sqrt(a*X+b) ;
run ;


proc iml ;
use a;
read all into xy ;


n = nrow(xy) ;
y = xy[,2] ;
x = xy[,1] ;


bho = {1,1} ;
diff=1010101010101010101;

do i = 1 to 10 while (diff>0.0001) ;

t1 = sqrt(bho[1,1]*x+bho[2,1]) ;
t2 = bho[1,1]*x+bho[2,1] ;


sa = -sum(x#y/t1)+sum(x) ;
sb = -sum(y/t1)+n ;

q= sa // sb ;


saa = sum(0.5*(x#x)#y#(t2##(-3/2))) ;
sab = sum(0.5*x#y/(t2##(3/2))) ;
sbb = sum(0.5*y/(t2##(3/2))) ;

h = (saa || sab) // (sab || sbb ) ;


bhn = bho - inv(h)*q ;

diff = max(abs(bhn-bho)) ;

sse = (y-sqrt(bho[1,1]*x+bho[2,1]))`*(y-sqrt(bho[1,1]*x+bho[2,1])) ;

print i  q h bho bhn diff sse;


res = res // (i || bhn` || sse );

bho=bhn ;

end ;


nm = {"i" "th1" "th2" "SSE"} ;
print res[colname=nm] ;

create ps from res[colname=nm] ;
 append from res ;




pdat  = x || y || t1 ;
nm1 = {"x" "y" "yh"} ;

create pdat from pdat[colname=nm1] ;
append from pdat ;



quit ;


proc gplot data=pdat;
plot y*x yh*x/overlay;
symbol1 v=o i=none c=red;
symbol2 v=dot i=join c=blue ;
run;


proc gplot data=ps;
plot sse*i ;
symbol1 v=dot i=join c=blue ;
run;
