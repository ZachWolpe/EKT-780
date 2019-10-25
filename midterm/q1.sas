data a ;
do t = 0 to 200 ;
a= 5160;
b= 0.5 ;
c= 0.2 ;
 p1 = a*(1-(b*exp(-c*t)));
 output ;
end ;
run ;

symbol1 i=join color=blue  height=.5;
symbol2 i=join color=green;
symbol3 i=join color=purple;


/* proc gplot data=a ;
plot (p1)*t /overlay;
run ; */
proc sgplot data=a;
	scatter y=p1 x=t;
run;
quit ;

proc iml ;
n=120 ;
t=(0:n)`/4 ;
a= 5160 ;
b= 0.8 ;
c= 0.1 ;

pt = t || a*(1-b*exp(-c*t)) ;

*print pt ;

l=2 ;
u=20 ;

h =a ;
w = u-l ;

ss = 10000 ;

x = (ranuni(J(ss,1,0))*w + J(ss,1,l)) // t ; 
y = ranuni(J(ss,1,0))*h  // pt[,2] ; 

pti =  y || a*(1-b*exp(-c*x)) ;
cnt = (pti[,1] <=  pti[,2]) ;
t_cnt = cnt[+] ;

print  t_cnt ;*x y pti  cnt t_cnt ;
aa=(h*w) ;

int = cnt[:]*(h*w) ;

print t_cnt aa  int ;

nm={"t" "p" "pti" "cnt"} ;
pdat = x || pti || cnt ;


pmax = a*(1-b*exp(-c*u)) ;
print pmax ;

create pdat from pdat[colname=nm] ;
append from pdat ;

quit ;

symbol1 i=none color=blue  height=.5;
symbol2 i=join color=green  height=.5;

proc sort data=pdat ;
 by t ;
run ;


/* proc gplot data=pdat ;
plot (p pti ) * t / overlay;
run ; */
proc sgplot data=a;
	scatter y=p x=t;
	scatter y=pti x=t;
	



quit ;
