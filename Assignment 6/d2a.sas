options nosource ;
title ;


%macro mcint(a,b,ss) ;

proc iml ;

start eval(x) ;
y=  10 +6.32972*x	-1.72728*(x##2)
     +0.2017*(x##3)	-0.00996*(x##4) +	0.00017 *(x##5) ;
return(y) ;
finish eval ;

x=((1:2000)/100)`;

y=eval(x) ;


xy = x || y ;
nm = {"x" "y"} ;

create a from xy[colname=nm] ;
append from xy ;


maxy = round((xy[,2])[<>]+1) ;
*print maxy ;

a=&a ;
b=&b ;


ss=&ss ;
ds = J(ss,1,0);

xc = ranuni(ds)*(b-a)+a ;
yc = ranuni(ds)*maxy ;

yy = eval(xc) ;

cprop = (yc<yy)[+]/ss ;
area = cprop*((b-a)*maxy) ;

*print  xc yc yy cprop area;

datall = (xy || J(nrow(x),1,5) )  //
          (xc || yc || (yc<yy) ) ;

print area ;

nm1 = {"x" "y" "gr"} ;
create dall from datall[colname=nm1] ;
append from datall ;

create area from area[colname="area"] ;
append from area;

quit ; 
%mend ;


%mcint(2,10,10000) ;
* lower boundary , upper boundary , number of points to generate ;


proc gplot data=dall ;
plot y*x=gr ;
run ;

quit ;
