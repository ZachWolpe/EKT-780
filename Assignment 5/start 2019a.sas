libname ekt "c:\departement\ekt720\knearest" ;

symbol1 value=dot 
        height=2 i=none;
symbol2 value=dot 
        height=2 i=none;
symbol3 value=dot 
        height=2 i=none;

proc gplot data=ekt.yx ;
plot x2*x1=group ;
title;
run ;

proc iml ;
 use ekt.yx ;
 read all  into xy ;
 xy=xy[1:20,];
 n=nrow(xy) ; 
 print xy ;

 y=xy[,1] ;
 x=xy[,2:3] ;
 *x=(x-x[:,])/std(x) ;
 n=nrow(xy) ;

 print x ;

 dist= distance(x) ;
 print dist ;

quit ;
