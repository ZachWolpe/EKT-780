options ls=150 ps=300;
PROC IMPORT OUT= WORK.a4xl 
            DATAFILE= "C:\SAS\User Data\Hons\EKT 720\xy_2015.xls" 
            DBMS=EXCEL REPLACE;
     SHEET="Use$"; 
     GETNAMES=YES;
     MIXED=NO;
     SCANTEXT=YES;
     USEDATE=YES;
     SCANTIME=YES;
RUN;

data ass4;
set a4xl;
run;

proc iml;
use ass4;
read all into XY;
n=nrow(XY);
k=ncol(XY);
Trainperc=0.8;
span=0.1;

m= span*n;
*****************Kernal Estimation**************;
rangen=j(n,1,0);
sampling=ranuni(rangen)||xy;
call sort(sampling,{1});
TrainingSet=sampling[1:round(Trainperc*n),2:3];
TestSet=sampling[round(Trainperc*n):n,2:3];
call sort(TrainingSet,{1});
print XY trainingset;
mini= min(trainingset[,1]);
maxi= max(trainingset[,1]);
nsample=nrow(trainingset);
print mini maxi nsample;
do method= 1 to 2;
do i= 1 to nsample;
xdev=abs(trainingset[,1]-trainingset[i,1]);
Step1=xdev||TrainingSet;
call sort(Step1,{1});
*print xdev Step1;
Window=Step1[1:m,];
*print Window;
h=0.5*(max(window[,2])-min(window[,2]));
*print h;
Step2=window[,1]/h;
Wt=((1-Step2##3)##3)#(Step2<1);
if method=1 then 
do;
*print wt;
yhat=yhat//(Wt`*Window[,3])/sum(Wt);
end;
if method=2 then 
do;
Wmat=diag(Wt);
X=J(nrow(window),1,1)||Window[,1];
Y=Window[,3];
*print X Wmat Y;
bhat=inv(X`*Wmat*X)*X`*Wmat*y;
print bhat;
yhatLPR=X*bhat;
yhatLPRs=yhatLPRs//yHatLPR[1];
print yhatLPRs;
end;
end;
end;
Lowessplotdata=Trainingset||yhat;
LPRplotdata=Trainingset||yhatLPRs;
print LPRplotdata Lowessplotdata;

create lowess from LowessPlotData[colname={'X','Y','Ynew'}];
append from LowessPlotData;

create LPR from LPRPlotData[colname={'X','Y','Ynew'}];
append from LPRPlotData;

quit;
/*
proc print data=lowess;
run;
*/

goptions reset = all;
symbol1 v=dot c=black h=1;
symbol2 i=join c=red h=1;
Title2 'Lowess plot';
proc gplot data=lowess;
plot y*x Ynew*x /overlay;
run;

goptions reset = all;
symbol1 v=dot c=black h=1;
symbol2 i=join c=red h=1;
Title2 'LPR plot';
proc gplot data=lPR;
plot y*x Ynew*x /overlay;
run;
