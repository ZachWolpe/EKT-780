PROC IMPORT OUT= WORK.cdata 
            DATAFILE= "C:\Users\IrvinK\Desktop\EKT(final)\DONE\Ass 6\cdata.xls" 
            DBMS=EXCEL REPLACE;
     RANGE="lreg"; 
     GETNAMES=YES;
     MIXED=NO;
     SCANTEXT=YES;
     USEDATE=YES;
     SCANTIME=YES;
RUN;
data a6;
set cdata;
run;
*proc print data=a6;
*run;
goptions reset=all;
axis1 label=(angle=90 'Y');
axis2 label=('X');
symbol1 color=navy line=3 value=dot width=1;
title1 'Plot of Cdata';
proc gplot data=a6;
plot y*x / vaxis=axis1 haxis=axis2;
run;
proc iml;
use a6;
read all into xy;

n=nrow(xy);
x=xy[,1];
y=xy[,2];
step=10;

/*using the graph printed in gplot statement*/
x_star1min=150;
x_star1max=255;
x_star2min=170;
x_star2max=310;

**Standard regression Module;
start std_reg;
bh=inv(xm`*xm)*xm`*y;
yhat=xm*bh;
sse=(y-yhat)`*(y-yhat);
finish std_reg;

**Bootstrap Module;
start boot;
	u=J(n,1,0);
	u=uniform(u);
	u=int(u#n)+1; /*correction if u=0*/
	bs=in[u,]; 
	if pp=1 then do;
	print u bs; 
	end;
finish boot;

/*Position routine from Buwang;
instead, can use call qntl('mtrx name', 'variable',{'interval'}, 'SAS option')*/
start pos;
	clpos = (alpha/2)/100*l ;
	cupos = l - (alpha/2)/100*l ;
	 
	lp1 = int(clpos) ;
	lp2 = int(clpos+1) ;
	if lp1 < 1 then lp1=1 ;

	up1 = int(cupos) ;
	up2 = int(cupos+1) ;
	if up2 > l then up2=l ;

	call sort(bres,st) ;

	cll = (bres[lp2,st]-bres[lp1,st])*(clpos-lp1)+bres[lp1,st] ;
	clu = (bres[up2,st]-bres[up1,st])*(cupos-up1)+bres[up1,st] ;
finish pos ;


start val; /*value search*/
ssm=9999999999999999999;
**value search of x pairs;
do x1_star=x_star1min to x_star1max by step;
do x2_star=x_star2min to x_star2max by step;
  xm=J(n,1,1)||x||(x-x1_star)#((x-x1_star)>0)||(x-x2_star)#((x-x2_star)>0); /*X matrix with conditions for dummies*/
**call std_reg;
bh=inv(xm`*xm)*xm`*y;
sse=(y-xm*bh)`*(y-xm*bh);

if sse < ssm then do;
	bh_min=bh;
	x_star1min=x1_star;
	x_star2min=x2_star;
	ssm=sse;
	end;
end;
end;
finish val;

start regression;
	call val ;
	if pp=1 then print "Regression - phase 1" ,
	step x_star1min x_star1max x_star2min x_star2max , bh_min x_star1min x_star2min ssm;

	sv=5;
	x_star1min=x_star1min-sv;
	x_star1max=x_star1min+sv;
	x_star2min=x_star2min-sv;
	x_star2max=x_star2min+sv;
	step=0.2;

	call val;
	if pp=1 then print "Regression - phase 2" ,
	step x_star1min x_star1max x_star2min x_star2max , bh_min x_star1min x_star2min ssm;
	sv=0.3;
	x_star1min=x_star1min-sv;
	x_star1max=x_star1min+sv;
	x_star2min=x_star2min-sv;
	x_star2max=x_star2min+sv;
	step=0.01;

	call val ;
	if pp=1 then print "Regression - phase 3" ,
	step x_star1min x_star1max x_star2min x_star2max , bh_min x_star1min x_star2min ssm;
finish regression;


pp=1 ;
print "Initial regression" ;
call regression ;

xyyh = x || y || yhat;
create xyyh from xyyh[colname={x y yhat}] ;
append from xyyh ;


pp=0; 
it=250 ;
do ii=1 to it ;
	in=xy ;
	step=1 ;
x_star1min=150;
x_star1max=255;
x_star2min=170;
x_star2max=310;
	call boot ;
	*print boot ;
	x=bs[,1] ;
	y=bs[,2] ;
	*print "Bootstrap" ;
	call regression ;
	bres1 = bres1 // (ii ||  bh_min` || x_star1min || x_star2min || ssm) ;
end ;


l=it ;
alpha=5;
nn={x1_star x2_star} ;
print l alpha ;
bres=bres1[,6:7] ;
do st = 1 to ncol(bres) ;
call pos ;
**cc = nn[,st] ;
*print st cc cll clu ;
end ;


print bres1[colname={nr bh1 bh2 bh3 bh4 x1_star x2_star sse}] ;

create assgn6_data from bres1[colname=nm] ;
append from bres1 ;
quit ;
PROC SGPLOT DATA =assgn6_data;
HISTOGRAM x1s ;
HISTOGRAM x2s;

TITLE "Histogram of structural break points";
RUN;
title ;

/*proc univariate data=cidat ;*/
/*var x1s x2s ;*/
/*histogram x1s x2s ;*/
/*run ;*/
goptions reset=all;
axis1 label=(angle=90 'Y');
axis2 label=('X');
symbol1 color=navy line=3 i=join value=dot width=1;
title1 'Plot of Cdata Regression';

proc gplot data=assgn6_data;
plot (y yh)*x / overlay ;
run ;

quit ;
