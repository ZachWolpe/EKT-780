quit ;
dm 'odsresults;clear';

options ls=85 nodate pageno=1 ;

libname ekt "c:\departement\ekt720\boot" ;

symbol1 interpol=none width=4
 		color=blue
        value=dot
        height=1;

symbol3 interpol=join width=3
 		color=green
        value=dot
        height=.3;

symbol4 interpol=join width=3
 		color=green
        value=dot
        height=.3;


proc gplot data=ekt.cdata;
plot y*x / href=175 255;
title "Structural break data" ;
run;

%macro sb(x1s,x2s,bs) ;

proc iml;
use ekt.cdata; 
read all into xy;
xy_orig=xy ;

start reg ;
k=ncol(x) ;
bh=inv(x`*x)*x`*y ;
res=y-x*bh ;
ess=res`*res ;
tss=(y-y[:,])[##,] ;
rss=tss-ess ;
df_m = ncol(x)-1 ;
df_e = n-ncol(x) ;
df_t = n-1 ;
ms_m = rss/df_m ;
ms_e = ess/df_e ;
f=ms_m/ms_e ;
p_f = 1-probf(f,df_m,df_e) ;

rmse = sqrt(ms_e) ;
meany = y[:,] ;
cv=rmse/meany*100 ;

r2=rss/tss ;
ar2 = 1 - ((1-r2)*((n-1)/(n-k))) ;


seb=sqrt(vecdiag(ms_e*inv(x`*x))) ;
t=bh/seb ;
p_t = 2*(1-probt(abs(t),df_e)) ;

cll = bh - tinv(0.975,df_e)*seb ;
clu = bh + tinv(0.975,df_e)*seb ;
finish reg ;

n=nrow(xy);
k=ncol(xy);
x=xy[,1];
y=xy[,2];

x1s=&x1s ;
x2s=&x2s ;

xm=J(n,1,1)||x||(x-x1s)#((x-x1s)>0)|| 
       (x-x2s)#((x-x2s)>0);
x=xm ;
call reg ;

print y xm bh r2 ;

yh_orig = x*bh ;
err = y - yh_orig ;

res = xy || (x*bh) ;
nm={"x" "y" "yh"} ;
	
create res from res[colname=nm] ;
append from res ;
close res ;

dt = x || y ;
nm1={"int" "x1" "xs1" "xs2" "y"} ;
create dt from dt[colname=nm1] ;
append from dt ;
close dt ;


it=&bs ;
do i = 1 to it ;

bserr= (sample(err,n))` ;
y = yh_orig +bserr ;
*print x y ;
call reg ;
*print bh ;

resbs = resbs // (bh` || r2) ;

resbs_yh =  resbs_yh || x*bh ;

end ;

nm={"bh0" "bh1" "bh2" "bh3" "r2"} ;
*print resbs[colname=nm] ;
mresbs = resbs[:,] ;
print mresbs[colname=nm] ;

create bres from resbs[colname=nm] ;
append from resbs ;
close bres;

alpha = 0.05 ;
perc = alpha/2 ||(1-alpha/2) ;

call qntl(ci,resbs,perc) ;

print perc ci ;

call qntl(clres_yh,resbs_yh`,perc) ;
clres_yh=clres_yh` ;



*print clres_yh[colname=nm1] ;

clresp =  xy_orig || yh_orig || clres_yh ;
nm2 = {"x" "y" "yh" "lcl" "ucl"} ;


create clresp from clresp[colname=nm2] ;
append from clresp ;
close clresp ;


quit;


proc univariate data=bres ;
histogram ;
run ;


proc sort data=clresp ;
 by x ;
run ;


proc gplot data=clresp;
plot (y yh lcl ucl)*x / overlay;
title "Structural break data" ;
title2 "x1s=&x1s and x2s=&x2s : Bootstrap iterations &bs" ;
run;


%mend ;



%sb(175,255,1000) ;
/*%sb(150,300,100) ;*/
/*%sb(200,220,100) ;*/
/*%sb(150,170,100) ;*/
/*%sb(280,290,100) ;*/


*the below is used for checking purposes ;

proc reg data =dt noprint ;
 model y = x1 xs1 xs2  / p r clm;
 output out=bbb  lclm=lcl uclm=ucl p=pred; ;
run ;

proc sort data=bbb;
by x1;
run ;

proc gplot data=bbb ;
plot (y pred lcl ucl)*x1 / overlay;
title "Structural break data" ;
title "Using normality assumption" ;
run;



quit ;



