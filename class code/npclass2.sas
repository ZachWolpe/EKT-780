options ls=72 nocenter ;

libname ekt "c:\departement\ekt720\non parametric" ;
title ;

proc gplot data=ekt.npdat2017 ;
plot y*x ;
run ;

proc iml ;
 use ekt.npdat2017 ; 
 read all into xy ;
 x1=xy[,1] ;
 y1=xy[,2] ;
 n =nrow(xy) ;
 m=7 ;

 print n m  ;


 do i = 1 to n ;
	  *	iteration selecting each observation as a focal
	 	point ;

	 f_point = x1[i,] ;
	 *print "AAAA" i   f_point  ;

	 
	 x0=x1-f_point ;
	 print x0 ;
	 axy=abs(x0) || x0 || xy ;

	 call sort(axy,{1}) ;
	 print axy ;

	 x=J(m,1,1) || axy[1:m,2] ;
	 y= axy[1:m,4] ;

	 bh=inv(x`*x)*x`*y ;
	 yh_fp = bh[1,1] ;
	 
	 print x y ;
	 print i bh  yh_fp ;
	 
	 yh = yh // yh_fp ;

 end ;

 res = xy || yh ;
 nm={"x" "y" "yh"} ;
 print res[colname=nm];

 create npres from res[colname=nm] ;
 append from res ;

quit ;

symbol1 interpol=none width=4
 		color=blue
        value=dot
        height=2;

symbol2 interpol=line width=3
 		color=red
        value=dot
        height=.3;

proc gplot data=npres ;
plot (y yh)*x / overlay;
run ;

quit ;
