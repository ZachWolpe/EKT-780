
* __________________________ SB memo Q3 BS 2019 __________________________ ;

dm 'odsresults; clear';
options ls=85 nodate pageno=1;
libname ekt '/folders/myfolders/sasuser.v94/EKT 720/Assignment 6';

symbol1 interpol=none
	width=4
	color=blue
	value=dot
	height=1;
		
symbol2 interpol=join width=3
	color=green
	value=dot
	height=0.3;
	
symbol3 interpol=join width=3
	color=green
	value=dot
	height=0.3;
	


proc sgplot data=ekt.cdata;
	scatter y=y x=x;
title 'Structural Break Data';
run;


%macro sb(x1s, x2s, bs);

	proc iml;
	use ekt.cdata;
	read all into xy;
	xy_orig = xy;
	
	start reg;
		k = ncol(x);
		bh = inv(x`*x)*x`*y;
		res = y - x*bh;
		ess = res`*res;
	 	tss = (y-y[:,])[##,]; * [##,] sum of squares;
	 	rss = tss-ess;
	 	df_m = k-1;
	 	df_e = n-k;
	 	df_t = n-1;
	 	ms_m = rss/df_m;
	 	ms_e = ess/df_e;
	 	f = ms_m/ms_e;
	 	p_f = 1 - probf(f,df_m,df_e);
	 	rmse = sqrt(ms_e);
	 	meany = y[:,];
	 	cv = rmse/meany*100;
	 	r2 = rss/tss;
	 	ar2 = 1 - ((1-r2)*((n-1)/(n-k)));

		seb = sqrt(vecdiag(ms_e*inv(x`*x)));
		t = bh/seb;
		p_t = 2*(1-probt(abs(t),df_e)) ;
		cll = bh - tinv(0.975, df_e)*seb;
		clu = bh + tinv(0.975, df_e)*seb;
	finish reg;

	
	n = nrow(xy);
	k = ncol(xy);
	x = xy[,1];
	y = xy[,2];
		
	x1s = &x1s;
	x2s = &x2s;
	
	xm = J(n,1,1)||x||(x-x1s)#((x-x1s)>0)||
		(x-x2s)#((x-x2s)>0);
	x = xm;
	call reg;
	
	print y xm bh r2;
	
	yh_orig = x*bh;
	err = y - yh_orig;
	
	res = xx || (x*bh);
	nm = {'x' 'y' 'yh'};
	
	create res from res[colname=nm];
	append from res;
	close res;
	
	dt = x||y;
	nm1 = {'int' 'x1' 'xs1' 'xs2' 'y'};
	
	create dt from dt[colname=nm];
	append from dt;
	close dt;
	
	it = &bs;
	
	do i=1 to it;
		bserr = (sample(err, n))`;
		y = yh_orig + bserr;
		call reg;
		resbs = resbs // (bh`||r2);
		resbs_yh = resbs_yh || x*bh;
	end;
	
	nm={'bh0' 'bh1' 'bh2' 'bh3' 'r2'};
	mresbs = resbs[:,];
	print mresbs[colname=nm];
	
	create bres from resbs[colname=nm];
	append from resbs;
	close bres;

	alpha = 0.05;
	perc = alpha/2 || (1-alpha/2);
	
	call qntl(ci,resbs,perc);
	print perc ci;

	call qntl(clres_yh, resbs_yh`, perc);
	clres_yh = clres_yh`;
	
	*print clres_yh[colname=nm1] ;
	
	clresp = xy_orig || yh_orig || clres_yh;
	nm2 = {'x' 'y' 'yh' 'lcl' 'ucl'};

	create clresp from clresp[colname=nm2];
	append from clresp;
	close clresp;
	quit;
	
	proc univariate data=bres;
	histogram;
	run;
	
	proc sort data=clresp;
		by x;
	run;

	proc sgplot data=clresp;
		scatter y=y x=x;
		scatter y=yh x=x;
		scatter y=lcl x=x;
		scatter y=ucl x=x;
		band x=x upper=ucl lower=lcl;
		title "Structural break data" ;
	title2 "x1s=&x1s and x2s=&x2s : Bootstrap iterations &bs" ;
		

	proc gplot data=clresp;
	plot (y yh lcl ucl)*x / overlay;
	title "Structural break data" ;
	title2 "x1s=&x1s and x2s=&x2s : Bootstrap iterations &bs" ;
	run;


%mend;


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

proc sgplot data=clresp;
	scatter y=y x=x;
	scatter y=yh x=x;
	scatter y=lcl x=x;
	scatter y=ucl x=x;
	band x=x upper=ucl lower=lcl;
	title "Structural break data";
	title2 "Using normality assumption";
run;	

quit ;











* __________________________ SB memo Q3 BS 2019 __________________________ ;

options ls=85 nodate pageno=1 ;

quit ;
dm 'odsresults;clear';
title ;


libname ekt '/folders/myfolders/sasuser.v94/EKT 720/Assignment 6';

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











* __________________________ SB memo Q4 NO BS 2019 __________________________ ;
quit;

dm 'odsresults; clear';
options ls=85 nodate pageno=1;
libname ekt '/folders/myfolders/sasuser.v94/EKT 720/Assignment 6';


symbol1 interpol=none width=4
 		color=blue
        value=dot
        height=1;

symbol2 interpol=join width=3
 		color=red
        value=dot
        height=.3;
        
proc sgplot data=ekt.cdata;
	scatter y=y x=x ;
title "Structural break data" ;
run;

proc gplot data=ekt.cdata;
plot y*x;
title "Structural break data" ;
run;



proc iml;
use ekt.cdata;
read all into xy;

start reg;
	k=ncol(x);
	bh=inv(x`*x)*x`*y;
	res=y-x*bh;
	ess=res`*res;
	tss=(y-y[:,])[##,];
	rss=tss-ess;
	df_m = ncol(x)-1;
	df_e = n-ncol(x);
	df_t = n-1 ;
	ms_m = rss/df_m ;
	ms_e = ess/df_e ;
	f = ms_m/ms_e;
	p_f = 1-probf(f,df_m,df_e);
	rmse=sqrt(ms_e);
	meany=y[:,];
	cv=rmse/meany*100;
	
	r2 = rss/tss;
	ar2 = 1 - ((1-r2)*((n-1)/(n-k)));
	seb = sqrt(vecdiag(ms_e*inv(x`*x)));
	p_t = 2*(1-probt(abs(t),df_e));
	
	cll = bh - tinv(0.975, df_e)*seb;
	clu = bh + tinv(0.975, df_e)*seb;
finish reg;


start bs1;
	u = sample(1:n,n);
	bs = in[u,];
	if pp=1 then do; print u bs; end;
finish bs1;


call sort(xy,{1});

	n=nrow(xy);
	k=ncol(xy);
	x=xy[,1];
	y=xy[,2];
	yinit=y ;
	xinit=x ;
	n=nrow(xy) ;

	c=10 ;
	print c n;
	essb=999999999999999999999999999999999999999999999;

do i=(c+1) to (n-2*c);
	do j=i+c+1 to n-c;
		grid = grid // (i||j);
	end;
end;	

print grid ;

*x1s=175 ;
*x2s=255 ;

do i=1 to nrow(grid);
	x1s_p = (grid[i,1]);
	x2s_p = (grid[i,2]);
	x1s = xinit[x1s_p,1];
	x2s = xinit[x2s_p,1];
	
	xm = J(n,1,1)||xinit||
			(xinit-x1s)#((xinit-x1s)>0)||
			(xinit-x2s)#((xinit=x2s)>0);
		
	x = xm;
	call reg;
	
	if ess < essb then do;
		x1s_pm = x1s_p;
		x2s_pm = x2s_p;
		x1s_m = x1s;
		x2s_m = x2s;
		x_m = x;
		r2_m = r2;
		bh_m = bh;
		essb=ess;
	end;
	
	anim_res = anim_res // (xinit||y||xm*bh||J(n,1,i));
end;
 
print 'BBBBB' x1s_pm x1s_m x2s_pm x2s_m r2_m bh_m x_m;
print essb;

res = xy || (x_m*bh_m);
nm = {'x' 'y' 'yh'};

create res from res[colname=nm];
append from res;

nm2 = {'x' 'y' 'yh' 'it'};
create anim_res from anim_res[colname=nm2];
append from anim_res;
quit;


proc sort data=res;
	by x;
run;

proc sgplot data=res;
	scatter y=y x=x;
	scatter y=yh x=x;
title "Structural break data";
run;

proc gplot data=res;
plot (y yh)*x / overlay;
title "Structural break data" ;
run;





