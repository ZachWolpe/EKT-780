data a (keep=x y) ;
s=0 ;
do i = 1 to 100 ;
x= i ;
y = sqrt((2.35*x+5.89))+rannor(5675)*s ;
output;
end ;
run ;
title ;
title2 ;
proc print data=a ;
run;
/*
proc gplot data=a;
plot y*x;
run;
*/
proc IML;
use A;
read all into M;

/******************regression****************/
start OLS;
do t1=t1min to t1max by step;
	do t2= t2min to t2max by step;

yhat=sqrt(t1#x+t2);

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
step=0.05;
run OLS;

h=0.1;
t1min=break1-h;
t1max=break1+h;
t2min=break2-h;
t2max=break2+h;
step=0.05;
run OLS;


h=0.01;
t1min=break1-h;
t1max=break1+h;
t2min=break2-h;
t2max=break2+h;
step=0.005;
run OLS;

finish ThreeStageReg;
/******************End of Three Stage Regression****************/

t1min=1;
t1max=3;
t2min=2;
t2max=7;
InputY=M[,2];
X=M[,1];
run threestagereg;
print break1 break2;


quit;


proc nlin data=a method=newton list;
parms theta1=1 theta2=1;
model y=sqrt(theta1*x+theta2);
run;

proc IML;********************* using XPX*****************************;
use a;
read all into M;
X=M[,1];
Y=M[,2];
obs=nrow(m);
bInit={1,1};

start NRprep;
XX=(0.5#((theta1#X+theta2)##-0.5)#X)||(0.5#((theta1#X+theta2)##-0.5));
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
yhat=sqrt((theta1*x+theta2));
SSE=ssq(y-yhat);
SSEdiff=abs(CrapSSE-SSE);
bNew=bOld+XXinv*XX`*(y-yhat);
print runs SSE SSEdiff bOld bNew ;
bOld=bNew;
CrapSSE=SSE;
end;
runs=runs-1;
print 'converges after' runs SSE bOld bNew;

quit;

proc IML;******************with bootstrap*********************;
use a;
read all into M;

print M;
start BS;
rannums=j(obs,1,0);
selection=int(uniform(rannums)#obs)+1;
Finish BS;
obs=nrow(m);

start NRprep;
XX=(0.5#((theta1#X+theta2)##-0.5)#X)||(0.5#((theta1#X+theta2)##-0.5));
XXinv=inv(XX`*XX);
finish NRprep;

rep=50;
free collect;
do j=1 to rep;
run bs;
mnew=M[selection,];
print Mnew;

X=Mnew[,1];
Y=Mnew[,2];

bInit={1,1};
bOld=bInit;

SSEdiff=99999999;
CrapSSE=99999999;
Tolerance=10##-8;
do runs=1 to 20 while (SSEdiff>Tolerance);
theta1=bOld[1];
theta2=bOld[2];
run NRprep;
yhat=sqrt((theta1*x+theta2));
SSE=ssq(y-yhat);
SSEdiff=abs(CrapSSE-SSE);
bNew=bOld+XXinv*XX`*(y-yhat);
print runs SSE SSEdiff bOld bNew ;
bOld=bNew;
CrapSSE=SSE;
end;
runs=runs-1;
print 'converges after' runs SSE bOld bNew;

collect=collect//(rep||runs||SSE||bNew`);
end;
print collect;
quit;
