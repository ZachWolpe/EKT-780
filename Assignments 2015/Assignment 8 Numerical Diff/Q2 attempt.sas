/*MI Masetla 12333710*/
data a ;
do i = 1 to 20 ;
x=i ;
y = 30 - 80*x/(exp(x)) +rannor(11567)*1 ;
output ;
end ;
run;
*proc nlin data=a; *check;
*parms t1=20 t2=-90;
*model y = t1 + t2*x/(exp(x));
*run;

*proc print data =a noobs;
*run ;

*proc plot data=a ;
*plot y*x ;
*run ;

proc iml;
use a;
read all into xy;
print xy;

x=xy[,2];
y=xy[,3];
*print x y;
n=max(xy[,1]);
*print n;
step=1;
error=10##17;

***************************Question 2***********************;
/*Grid search*/
do theta1=20 to 40 by step;
do theta2=-90 to -70 by step;
yhat=theta1+((theta2#x)/exp(x));
residual=(y-yhat);
mse=(residual`*residual)/(n-2);
if mse<=error then do;
error=mse;
est=theta1||theta2;
y_h=yhat;
end;
end;
end;

xyyh=xy||y_h;

print mse est xyyh;

/*Question 3;
Numerical Differentiation

start num_diff;
yhat=theta1+((theta2#x)/exp(x));
yhdiff1=theta2#exp(-x)-theta2#x#exp(x);
yhdiff2=-theta2#exp(-x)-theta2#exp(-x)+theta2#x#exp(-x);
finish num_diff;*/

create q2 from xyyh[colname={'x' 'y' 'yh'}];
append from xyyh;
quit;
symbol1 colour=red width=1 i=j; 
symbol2 colour=blue width=2 i=j;
proc gplot data=q2;
plot y*x yh*x/ overlay legend;
run;
