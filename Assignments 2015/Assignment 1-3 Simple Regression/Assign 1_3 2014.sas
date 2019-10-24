options nodate;
data a;
input income1 income2 house stand double prop;
cards;
  521.502    118.348    735.779     920.53      0     1215.86
   14.116    457.801    413.522     690.15      0      917.67
  308.237    205.341    567.238     903.11      1     1201.16
  449.589    157.470    496.226     659.05      1     1099.73
   47.286    555.871    414.292     769.92      0     1077.25
   12.702    400.744    283.223     539.62      1      973.28
  303.539    360.630    671.121     934.58      0     1206.01
  325.548    369.610    284.473     707.85      0     1071.03
  328.079     17.192    492.552     699.23      0      834.60
  479.735     34.212    767.408    1097.32      0     1102.11
   70.381    319.148    373.140     760.19      0      774.12
  232.232    255.517    238.515     577.39      0      828.64
   56.125    326.705    589.865     930.42      0     1020.10
  510.569     36.773    461.059     920.65      0     1044.39
   15.890    353.851    345.385     655.05      1      932.65
  298.906    126.398    531.592    1093.24      0     1036.68
  280.401    105.089    497.296     727.87      0      867.00
  188.411    419.229    383.097     903.32      0     1114.49
   11.004    462.602    351.969     575.44      0      833.46
  408.952    119.757    650.882     950.26      0     1044.39
  114.999    253.868    439.853     849.32      0      831.60
  200.932    141.234    400.907     571.64      0      773.70
  276.907    350.366    554.191     948.33      1     1273.24
  271.076    109.235    734.862     970.72      1     1124.66
  357.141    324.151    507.147     686.02      0     1146.86
   74.029    403.535    372.881     520.79      0      836.79
  112.752    195.755    550.987    1048.71      0     1023.95
  189.496    273.100    400.458     550.31      0      834.02
  283.516    395.697    445.404     600.35      0     1064.84
  255.701    154.743    535.123    1078.51      0     1075.30
  ;
run ;
**Assignment 2 coding**;
proc iml;
use a;
read all into b;
n=nrow(b);
y=b[,6];
income=b[,1]+b[,2];
house=b[,3];
stand=b[,4];
double=b[,5];
ratio=house/stand;
x=j(n,1,1)||income||stand||ratio||double;
**print x;

print 'Assignment 2a';
xpxinv=inv(x`*x);
bhat=xpxinv*x`*y;
print bhat;

print 'Assignment 2b';
k=5 /*number of parameters(betas)*/;
j=j(n,1,1)*j(n,1,1)`;
ssto=y`*y-(1/n)*y`*j*y;
sse=y`*y-bhat`*x`*y;
ssr=ssto-sse;
rsq=ssr/ssto;
radsq=1-((n-1)/(n-k))*(sse/ssto);
print rsq, radsq;

print 'Assignment 2c';
msr=ssr/(k-1);
mse=sse/(n-k);
fsig=msr/mse;
print fsig;

print 'Assignment 2d';
covm=mse*inv(x`*x);
t_b3=bhat[3,1]/covm[3,3]**(0.5);
print t_b3;

print 'Assignment 2e(i)';
uplim=bhat[4]+tinv(0.95,n-k)*covm[4,4]**(0.5);
lolim=bhat[4]-tinv(0.95,n-k)*covm[4,4]**(0.5);
print lolim uplim;

print 'Assignment 2e(ii)';
ybar=sum(y)/n;
p=3 /*number of parameters in restricted model*/;
x_res=x[,1]||x[,3]||x[,5];
bhat_res=inv(x_res`*x_res)*x_res`*y;
yhat_res=x_res*bhat_res;
sse_res=y`*y-bhat_res`*x_res`*y;
**ssr_res1=ssq(yhat_res-ybar);
rsq_res=1-sse_res/ssto;
**rsq_res1=ssr_res1/ssto;
f_res=((rsq - rsq_res)/2)/((1-rsq)/(n-k));
print f_res;

**Assignment 3 coding**;
proc iml;
use a;
read all into c;
n=nrow(c);
y=c[,6];
income=c[,1]+c[,2];
house=c[,3];
stand=c[,4];
double=c[,5];
ratio=house/stand;
x=j(n,1,1)||income||stand||ratio||double;
yx=y||x[,2:5];

print '3a';
bhat=inv(x`*x)*x`*y;
print bhat;

result=j(n,5,.);

print '3b';
m=y||x[,2:5];
do i=1 to n;
if i=1 then lout=m[2:n,];  /*if removing the first obs*/
else if i=n then lout=m[1:n-1,];  /*if removing the last obs*/
else lout=m[1:i-1,]//m[i+1:n,]; /*if removing obs in between, concatinate vertically*/

y_rem=lout[,1];
x_rem=j(n-1,1,1)||lout[,2:5];
bhat_rem=inv(x_rem`*x_rem)*x_rem`*y_rem;
**print bhat_rem;


result[i,]=bhat_rem`;
end;
print result;

print '3c';
aresult=result[:,]; /*average*/
print aresult;

print '3d';
result_d=j(n,5,.);
**u=sample(1:30,30);
**s=x[u];
do i=1 to n;
b=j(n,1,0);
draw=int(uniform(b)*n)+1;
sample=m[draw,];
**print draw sample;
x_sam=j(n,1,1)||sample[,2:5];
y_sam=sample[,1];
bhat_d=inv(x_sam`*x_sam)*x_sam`*y_sam;
result_d[i,]=bhat_d`;
end;
print result_d;

print '3d average';
aresult_d=result_d[:,];
print aresult_d;

print '3e';
result_e=j(100,5,.);

do i=1 to 100;
b1=j(n,1,0);
draw1=int(uniform(b1)*n)+1;
sample1=m[draw1,];

**print draw sample;
xe_sam=j(n,1,1)||sample1[,2:5];
ye_sam=sample1[,1];
bhat_e=inv(xe_sam`*xe_sam)*xe_sam`*ye_sam;
result_e[i,]=bhat_e`;
end;
print result_e;

print '3e average';
aresult_e=result_e[:,];
print aresult_e;


