options ls=85 nodate pageno=1 ;

quit ;
dm 'odsresults;clear';
title ;


libname ekt "c:\departement\ekt720\boot" ;

symbol1 interpol=none width=4
 		color=blue
        value=dot
        height=1;

symbol2 interpol=join width=3
 		color=red
        value=dot
        height=.3;


proc gplot data=ekt.cdata;
plot y*x;
title "Structural break data" ;
run; 

%macro sb(x1s,x2s) ;

proc iml;
use ekt.cdata; 
read all into xy;

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

xm=J(n,1,1)||x||
       (x-x1s)#((x-x1s)>0)|| 
       (x-x2s)#((x-x2s)>0);
x=xm ;
call reg ;

print y xm bh r2 , cll clu ;

res = xy || (x*bh) ;
nm={"x" "y" "yh"} ;
	
create res from res[colname=nm] ;
append from res ;

quit;


proc sort data=res ;
 by x ;
run ;

proc gplot data=res;
plot (y yh)*x / overlay;
title "Structural break data" ;
title2 "x1s=&x1s and x2s=&x2s" ;
run;
QUIT ;
%mend ;


%sb(175		,	255.02) ;
%sb(170.041	,	280.666) ;
%sb(201.13	,	223.253) ;
%sb(153.694	,	170.04) ;
%sb(280.666	,	290.579) ;




