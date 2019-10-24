proc iml;
use ass4;
read all into EKT4;
n =nrow(ekt4);
in=J(n,1,0);
ektin = ranuni(in) || ekt4;
call sort(ektin,{1}); /*sorts values in column 1 of xys in ascending order*/  
ektin1 = ektin[ 1: round(0.8*n),2:3]; /* round is for giving an integer value and thus this takes
										80% of the random selection of coloumns 2 and 3 */
x1 =ektin1[ ,1]; 
y1 = ektin1[ ,2];
ektin2 = ektin[round(0.8*n)+1:n,2:3];  /* for the testing matrix will give us remaining 20% of the matrix */
call sort(ektin1,{1}); 
call sort(ektin2,{1}); 
n2 =nrow(ektin1);
x =ektin1[ ,1];
y = ektin1[ ,2];
s =0.4;   /* we choose a small span as a high span would defeat the purpose of the NPR. Note a
				span of 0.4 means that we taking 40% of the observations close to our x0 */ 
xmax = max(x); 
xmin = min(x); 
m=round(s*n2); 
wght=J(m,1,.); 

/*part a */
x11 =ekt4[ ,1];
y11 = ekt4[ ,2];
title " Scatter plot of the initial x and y values ";
run Scatter ( x11 ,y11) label = { " X " " Y "} ;

                                                                                                                                    
 do x0=xmin to xmax by 1 ;                                                                                                           
 nwx = x - x0;
 nwxd = abs(nwx);
 md = ( ektin1 || nwx || nwxd); /* neighbourhood */
 call sort(md,{4});
 md = md[1:m,];
 h = 0.5* (max(md[,1])-min(md[ ,1]));
 	do j =1 to m;
	con = abs(md[j,1] - x0)/h;
	wgt = (1- con**3)**3;
	wght[j,] = (con<1)#wgt;
	end;

	ker = ker // ( x0 || ( wght` * md[, 2] )/ (J(1,m,1)*wght));
  /*LPR being implemented here */  
	xlp = J(m,1,1) || md[ ,3];
	ylp = md[ ,2];
	wlp = diag(wght);
	bhlp = inv(xlp` *wlp*xlp)*xlp`*wlp*ylp;
	poly = poly //(x0 || bhlp[1,1]); /* where bhlp[1,1] is the y estimate*/                                                                   
    lyp = lyp //(x0 || bhlp[1,1]); 
	xplx = xlp` * xlp;
	wxplwx = xlp`*wlp*xlp;
	ls = xplx[1,1];
	sig = (ylp-(xlp*bhlp))`*wlp* (ylp-(xlp*bhlp))/(m-2);
	cov = sig*inv(wxplwx);
	se = sqrt(cov[1,1]);
	cond =abs(tinv(0.975,m-2));
	upp = bhlp[1,1]+(cond*se);
	low = bhlp[1,1]-(cond*se);
	mdu = mdu //(x0 || low || 4)//(x0 || upp|| 5);
end;                                                                                                                        
r1 =nrow(ekt4);
r2 =nrow(lyp);
r3 = nrow(mdu);

you = { x yhat z};
F = (ker || J(nrow(ker),1,1))// ( lyp || J(nrow(lyp),1,2))// (x || y || J(nrow(x),1,3))// mdu;


 create dat1 from F[colname = you];
 append from F;
 quit;

 data dat1;
 set dat1;
 if z=1 then d="Kernel Estimaion";
 if z=2 then d="Local Polynomial Regression" ;
 if z=3 then d="Observed data" ;
 if z=4 then d="95% Lower Prediction Limit" ;
 if z=5 then d="95% Upper Prediction Limit" ;
run;
goptions reset all;
symbol1 value=square interpol=join height=1.1;
symbol2 value=star interpol=join height=1.1;
symbol3 value=dot  interpol=join height=1.1;

symbol1 value=square interpol=join height=1.1;
symbol2 value=star interpol=join height=1.1;
symbol3 value=dot  interpol=join height=1.1;

proc gplot data=dat1;
 plot yhat*x=d;
 where (z=1 or 4<=z<=5) ;
 title " Kernel estimation and Prediction Intervals" ;
 run;
quit;

proc gplot data=dat1;
 plot yhat*x=d;
 where (z=2 or 4<=z<=5);
 title " Local Polynomial Regression and Prediction Intervals";
 run;
quit;


proc gplot data=data1;
 plot yhat*x=d;
 where 3<=z<=5;
 title "Observed Data and Prediction Intervals";
 run;
quit;


 
