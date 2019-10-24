quit ;
dm 'odsresults;clear';


options ls=85 nodate pageno=1 ;

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

start bs1 ;
u=sample(1:n,n) ;
bs = in[u,] ;
if pp=1 then do; print u bs; end;
finish bs1 ;

call sort(xy,{1}) ;


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

do i = (c+1) to (n-2*c) ;
  do j = i+c+1 to n-c ;
  grid = grid // (i ||j) ;
  end ;
end ;

print grid ;

*x1s=175 ;
*x2s=255 ;


do i = 1 to nrow(grid) ;


x1s_p = (grid[i,1]) ;
x2s_p = (grid[i,2]) ;
x1s = xinit[x1s_p,1] ;
x2s = xinit[x2s_p,1] ;


*print x1s_p x1s x2s_p x2s ;


xm= J(n,1,1) || xinit ||
    (xinit-x1s)#((xinit-x1s)>0)|| 
    (xinit-x2s)#((xinit-x2s)>0);

* Just to illustrate ;
*xm=J(n,1,1)||xinit||xinit##2 || xinit##3 ||
       ((xinit-x1s)#((xinit-x1s)>0))##3 ||
       ((xinit-x2s)#((xinit-x2s)>0))##3 ;

x=xm ;
call reg ;

*print i  x1s x2s bh r2 ;

if ess < essb then do ;
 x1s_pm = x1s_p ;
 x2s_pm = x2s_p ;
 x1s_m = x1s ;
 x2s_m = x2s ;
 x_m = x ;
 r2_m = r2 ;
 bh_m = bh ;
 essb=ess ;
 *print "BBBBB" x1s_m x2s_m r2_m bh_m ;
end ;

anim_res = anim_res // (xinit || y || xm*bh  || J(n,1,i))  ;

end ;

print "BBBBB" x1s_pm x1s_m  x2s_pm x2s_m r2_m bh_m  x_m;

print essb ;

res = xy || (x_m*bh_m) ;
nm={"x" "y" "yh"} ;

create res from res[colname=nm] ;
append from res ;

nm2 = {"x" "y" "yh" "it"} ;
create anim_res from anim_res[colname=nm2] ;
append from anim_res ;

quit;


proc sort data=res ;
 by x ;
run ;

proc gplot data=res;
plot (y yh)*x / overlay;
title "Structural break data" ;
run;





quit ;


