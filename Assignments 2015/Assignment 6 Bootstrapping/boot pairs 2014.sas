options ls=72 nodate pageno=1 ;

data a (keep=x y) ;

do i = 1 to 100 ;
 x=ranuni(567)*200+300 ;
 u=rannor(4567)*800 ;
 y = -15 + 25*x + u ; /*benchmark model created by SM*/
 output ;
end ;

run ; 


proc print data=a ;
run ;

proc reg data=a ;
model y = x / clb;
run ;


proc iml ;

print "Bootstrap X,Y pairs" ;

use a ;
read all into xy ;

*print xy ;

n=nrow(xy) ;
x=J(n,1,1) || xy[,1] ;
y= xy[,2] ;

*regression module - sub_routine;
start reg ;
k=ncol(x);
bh = (inv(x`*x)*x`*y) ;
sse =ssq(y-x*bh) ;
cssy=ssq(y-(sum(y)/n)); /*total sum of squares*/         
r2 =  (cssy-sse)/cssy;
F = (r2/(k-1))/((1-r2)/(n-k)) ;
if pp=1 then do ; /*prints if pp=1*/
print bh r2 f ; 
end;
finish reg ;

*bootstrap routine ;
start bs1 ;
u=j(n,1,0);
u=uniform(u) ;
u=int(u*n)+1 ;/*for the case when u*n=0*/
bs = in[u,] ;/*selection matrix -> resample a percentage of n*/
if pp=1 then do; 
print u bs; 
end;
finish bs1 ;


*position routine ; /*determines upper and lower, percentiles etc*/
start pos;
clpos = (alpha/2)/100*l ; /*lower position of value*/
cupos = l - (alpha/2)/100*l ; /*upper position of value*/
  
lp1 = int(clpos) ;
lp2 = int(clpos+1) ;
if lp1 < 1 then lp1=1 ; /*if position is less than 1*/

up1 = int(cupos) ;
up2 = int(cupos+1) ;
if up2 > l then up2=l ;

call sort(bres,st) ; /*st=number of columns in bres*/

cll = (bres[lp2,st]-bres[lp1,st])*(clpos-lp1)+bres[lp1,st] ; /*lower limit*/
clu = (bres[up2,st]-bres[up1,st])*(cupos-up1)+bres[up1,st] ; /*upper limit*/
finish pos ;


pp=1 ;
print "Model" ;
call reg ;


pp=0 ;

it=550 ; /*iterations*/

do i = 1 to it ;
 in=xy ;
 call bs1 ;
 x= J(n,1,1) || bs[,1] ;
 y= bs[,2] ;
* print x y ;
* print "BS" i;
 call reg ;
 bh_sum = bh_sum // (bh` || r2) ;
end;

nm={"bh1" "bh2" "r2"} ;

print bh_sum[colname=nm]  ;

nm={"bh1" "bh2" "r2"} ;
create b from bh_sum[colname=nm] ;
append from bh_sum ;

l=it ;
alpha=5;
print l alpha ;
bres=bh_sum ;
do st = 1 to ncol(bres) ;
call pos ;
cc = nm[,st] ;
print st cc cll clu ;
end ;


quit ;

proc univariate data=b normal ;
 var bh1 bh2 r2;
 histogram bh1 bh2 r2;
run ;


proc sort data =b ;
 by r2 ;
 run ;

