options ls=120 nodate ps=1000 pageno=1 nocenter ;
libname ekt "c:\departement\ekt720\non parametric" ;
title ;

data xy_data; 
input x y;
datalines;
10 44
11 34
12 19
13 14
14 24
15 29
16 26
17 22
18 6
19 18
20 19
21 23
22 31
23 31
24 38
25 39
26 41
27 42
28 36
29 28
30 24
31 18
32 8
33 25
34 29
35 30
36 30
37 34
38 32
39 36
;
run;

%macro lpr(mpar) ;

proc iml ;
*use ekt.npdat2017 ;
use xy_data;
read all into xy ;

ix=xy[,1] ;
iy=xy[,2] ;
n=nrow(ix) ;

xm = repeat(ix,1,n) ;
* print xm;
d=(xm-xm`) ;
* print d;


m=&mpar ;
s=m/n ;

do i = 1 to n ;
xw=abs(d[,i]) || d[,i] || iy ;
call sort(xw,{1}) ;
xww=xw[1:m,] ;

print xw xww ;

call sort(xww,{2}) ;
hxww = xww[m,2]-xww[1,2] ;
h=hxww/2 ;
*h=8 ;

*print hxww h ;
az=abs(d[,i]/h) ;

print az ;


tc = ((1-az##3)##3)#(az<1) ;
print tc ;

w=diag(tc) ;
x=J(n,1,1) ||d[,i] ;
k=ncol(x) ;
print k ;
bh = inv(x`*w*x)*x`*w*iy ;
n1= (tc>0)[+,] ; print n1 ;
s2h = (iy-x*bh)`*w*(iy-x*bh)/(n1-k) ;
sebh= sqrt(vecdiag(s2h*inv(x`*w*x))) ;
lcl = bh[1,1 ]-sebh[1,1]*tinv(0.975,n1-k);
ucl = bh[1,1 ]+sebh[1,1]*tinv(0.975,n1-k);


*print i bh;

fv = fv // (bh[1,1 ] || lcl || ucl) ;
print (bh[1,1 ] || lcl || ucl) ;
end ;

print fv ;



moddat = xy || fv ;
nm={"X" "Y" "YH" "lcl" "ucl"};

create ekt.moddat from moddat[colname=nm] ;
append from moddat ;



quit ;


symbol1 value=dot
        height=.5 i=join;

proc sgplot data=ekt.moddat ;
band x=x lower=lcl upper=ucl / transparency=.5 legendlabel="Confidence Band" fillattrs=(color=purple);
series x=x y=y / lineattrs=(color=blue);
series x=x y=yh / lineattrs=(color=red);
Title "LPR: m= &mpar" ;
run ;

%mend ;

%lpr(10) ;
