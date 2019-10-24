/*MI Masetla 12333710*/
data a (keep=x y) ;
s=0 ;
do i = 1 to 100 ;
x= i ;
y = sqrt((2.35*x+5.89))+rannor(5675)*s ;
output;
end ;
run ;
*proc print noobs;
*run;
proc iml;
use a;
read all into xy;

x=xy[,1];
y=xy[,2];
*print x y;
n=nrow(xy);
*print n;
step=0.01;
error=10##12;
*result=j(4,1,.);

/*Grid search*/
do theta1=1 to 3 by step;
do theta2=2 to 7 by step;
*grid=grid//(theta1||theta2);
yhat=sqrt((theta1#x)+theta2);
residual=(y-yhat);
mse=(residual`*residual)/(n-2);
if mse<=error then do;
error=mse;
*result=(theta1//theta2//mse);
est=theta1||theta2;
y_h=yhat;
end;
end;
end;

xyyh=xy||y_h;

*print result; 
print mse est xyyh /*grid*/;

create q1 from xyyh[colname={'x' 'y' 'yh'}];
append from xyyh;
quit;
symbol1 value=star colour=red width=5 i=j; 
symbol2 colour=black width=2 i=j;
proc gplot data=q1;
plot y*x yh*x/ overlay legend;
run;
