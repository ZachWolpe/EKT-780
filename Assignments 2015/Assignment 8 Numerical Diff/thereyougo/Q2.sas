options ls=200 ps=1000;
data a (keep=x y) ;
do i = 1 to 20 ;
x=i ;
y = 30 - 80*x/(exp(x)) +rannor(11567)*1 ;
output ;
end ;
run;

proc IML;
use A;
read all into M;
print M;

/******************regression****************/
start OLS;
do t1=t1min to t1max by step;
	do t2= t2min to t2max by step;

yhat=t1 + t2*x/(exp(x));

SSE=ssq(InputY-Yhat);
if SSE<BestSSE then; 
	do;
	BestSSE=SSE;
 	break1=t1;
 	break2=t2;
	end;
	end;
end;
 	
finish;
/***********************End of regression********************/
/******************Three Stage Regression****************/
start ThreeStageReg;
BestSSE=99999999;
regoutput=0;
step=2;
run OLS;

h=5;
t1min=break1-h;
t1max=break1+h;
t2min=break2-h;
t2max=break2+h;
step=1;
run OLS;


h=3;
t1min=break1-h;
t1max=break1+h;
t2min=break2-h;
t2max=break2+h;
step=0.5;
run OLS;

finish ThreeStageReg;
/******************End of Three Stage Regression****************/

t1min=20;
t1max=40;
t2min=-90;
t2max=-70;
InputY=M[,2];
X=M[,1];
run threestagereg;
print break1 break2;


quit;


proc IML;
use A;
read all into M;
X=M[,1];
Y=M[,2];

thetainit={1,1};
thOld=thetainit;

start eval(x,y,theta1,t2,funcval) ;
yhat=theta1+t2#x/(exp(x));
funcval=(y-yhat)`*(y-yhat);
finish eval ;


diff=999999;
do i=1 to 20 while (diff>0.000001);
h=0.00001;
theta1=thOld[1];
t2=thOld[2];
run eval(x,y,theta1,t2,SSE);
print i SSE;

run eval(x,y,theta1,t2,Baseval);
run eval(x,y,theta1+h,t2,difft1);
dt1=(difft1-baseval)/h;
run eval(x,y,theta1,t2+h,difft2);
dt2=(difft2-baseval)/h;
m=dt1//dt2;
run eval(x,y,theta1+2*h,t2,ddt1component);
ddt1=(ddt1component-2*difft1+baseval)/(h**2);
run eval(x,y,theta1,t2+2*h,ddt2component);
ddt2=(ddt2component-2*difft2+baseval)/(h**2);
run eval(x,y,theta1+h,t2+h,dt1dt2component);
dt1dt2=(dt1dt2component-difft1-difft2+baseval)/(h**2);

h=(ddt1||dt1dt2)//(dt1dt2||ddt2);

thNew=thOld-inv(h)*m;

diff=sum(abs(thnew-thold));
thOld=thNew;

end;

print thNew;

quit;

proc nlin data=a method=newton list;
parms theta1=1 theta2=1;
model y=theta1+theta2*x/(exp(x));
run;


proc IML;
use a;
read all into M;
X=M[,1];
Y=M[,2];

bInit={1,1};

start NRprep;
XX=j(nrow(m),1,1)||(x/exp(x));
XXinv=inv(XX`*XX);
finish NRprep;

bOld=bInit;

SSEdiff=99999999;
CrapSSE=99999999;
Tolerance=10##-8;
do runs=1 to 20 while (SSEdiff>Tolerance);
theta1=bOld[1];
theta2=bOld[2];
run NRprep;
yhat=theta1+theta2*x/(exp(x));
SSE=ssq(y-yhat);
SSEdiff=abs(CrapSSE-SSE);
bNew=bOld+XXinv*XX`*(y-yhat);
print runs SSE SSEdiff bOld bNew ;
bOld=bNew;
CrapSSE=SSE;
end;
runs=runs-1;
print 'converges after' runs SSE bOld bNew;

start ols(beta) global(y,x);
residu=y-(beta[1]+beta[2]*x/(exp(x))); 
scr=residu[##];
return(scr);
finish; 
beta0={1,1}; 
optn={0 4}; 
call nlpnra(rc,beta,'ols',beta0,optn); 
print beta; 


quit;
