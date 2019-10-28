
* ____________ ____________________ Question 1: Monte Carlo Integration ____________________ ____________;

proc iml;
title 'Question 1';

start function(x);
	fx = 10 + 6.32972#x - 1.72728*(x##2) + 0.2017#(x##3) 
			- 0.00996*(x##4) + 0.00017#(x##5);
	return fx;
finish function;


start monte_carlo_integration;

	res = J(iters,4,999);	
	do i=1 to iters;
	
		* locate y max;
		t = do(x_low, x_up,0.05)`;
		y_max = max(function(t))*(1+0.15);
	
		* sample;
		y = rand('uniform')*(y_max - 0);
		x = rand('uniform')*(x_up - x_low);
		
		fx = function(x);
		
		if fx > y then; res[i,] = (0||x||y||fx);
		if fx<=y then; res[i,] = (1||x||y||fx);
		
	end;
	
	total_area = (y_max - 0)*(x_up - x_low);
	mcmc_int = mean(res[,1])*total_area;
	
finish monte_carlo_integration;
	


iters = 10000;
title ' Question 1.1';
x_low = 3; 
x_up = 8;
call monte_carlo_integration;
print 'Area Under Curve: ' (mcmc_int);

title 'Question 1.2'; 
x_low = 1; 
x_up = 10;
q12 = monte_carlo_integration;
print 'Area Under Curve: ' (mcmc_int);

title 'Question 1.3';
x_low = 0; 
x_up = 20;
call monte_carlo_integration;
print 'Area Under Curve: ' (mcmc_int);


create res from res[colname={'score' 'x' 'y' 'fx'}];
append from res;
quit;


proc sgplot data=res;
	scatter x=x y=y / group=score;
	scatter y=fx x=x;
title 'Monte Carlo Integration';
run;






* ____________ ____________________ Question 2: Bootstrap ____________________ ____________;

data cdata;
set '/folders/myfolders/sasuser.v94/EKT 720/Assignment 6/cdata.sas7bdat';
run;
proc print data=cdata (obs=10);
run;	




proc iml;
use cdata; read all into xy_orig;
n = nrow(xy_orig);


start regression;
	b = inv(x`*x)*x`*y;
	yh = x*b;
	err = y - yh;
	mse = ((y-yh)`*(y-yh))/(nrow(x) - ncol(x));
	ess = (yh - mean(y))`*(yh - mean(y));
	tss = (y - mean(y))`*(y - mean(y));
	r2 = ess/tss;
finish regression;



* ____ ____ ____ Bootstrap Residuals ____ ____ ____;
x = J(n,1,1) || xy_orig[,1]; y=xy_orig[,2];
call regression;
e = err;
yh_true = yh;

do i=1 to 3000;
	y = yh_true + sample(e, n, 'replace')`;
	call regression;
	res_r2 = res_r2 // r2;
end;




* ____ ____ ____ Bootstrap Residuals ____ ____ ____;
do i=1 to 3000;
	in = sample((1:n), n, 'replace')`;
	x = J(n,1,1) || xy_orig[in,1];
	y =  xy_orig[in,2];
	call regression;
	res_beta = res_beta // b`;
end;





* ____ ____ ____ Confidence Intervals ____ ____ ____;
n= nrow(res_r2);

call sort(res_r2);
lcl = res_r2[0.0275*n,];
ucl = res_r2[0.975*n,];

print (mean(res_r2)) lcl ucl;

res_r2 = res_r2 || J(nrow(res_r2), 1, lcl) || J(nrow(res_r2), 1, ucl);
create res_r2 from res_r2[colname={'r2' 'lcl' 'ucl'}];
append from res_r2;



call sort(res_beta, {1});
lcl_0 = res_beta[0.0275*n,1];
ucl_0 = res_beta[0.975*n,1];

call sort(res_beta, {2});
lcl_1 = res_beta[0.0275*n,2];
ucl_1 = res_beta[0.975*n,2];


res_beta = res_beta || 
	J(nrow(res_beta), 1, lcl_0) || J(nrow(res_beta), 1, ucl_0) ||
	J(nrow(res_beta), 1, lcl_1) || J(nrow(res_beta), 1, ucl_1);
create res_beta from res_beta[colname=
	{'b0' 'b1' 'lcl_0' 'ucl_0' 'lcl_1' 'ucl_1'}];
append from res_beta;



* ____ ____ ____ Results ____ ____ ____;
mean_r2 = mean(res_r2);
print mean_r2 lcl ucl;

mean_b0 = mean(res_beta[,1]);
print mean_b0 lcl_0 ucl_0;

mean_b1 = mean(res_beta[,2]);
print mean_b1 lcl_1 ucl_1;
quit;







* ____ ____ ____ Visualization ____ ____ ____;
proc sgplot data=res_r2;
	histogram r2; 
	refline lcl / axis=x lineattrs=(color=red);
	refline ucl / axis=x lineattrs=(color=red);
	title 'Bootstrap Errors: R2 Distribution';
run;

proc sgplot data=res_beta;
	histogram b0;
	refline lcl_0 / axis=x lineattrs=(color=red);
	refline ucl_0 / axis=x lineattrs=(color=red);
	title 'Bootstrap Pairs: Beta 0';
run;

proc sgplot data=res_beta;
	histogram b1;
	refline lcl_1 / axis=x lineattrs=(color=red);
	refline ucl_1 / axis=x lineattrs=(color=red);
	title 'Bootstrap Pairs: Beta 1';
run;











* ____________ ____________________ Question 3: Structural Breaks ____________________ ____________;

data cdata;
set '/folders/myfolders/sasuser.v94/EKT 720/Assignment 6/cdata.sas7bdat';
run;
proc print data=cdata (obs=10);
run;	


proc sgplot data=cdata;
	scatter x=x y=y;
	title 'Structural Break Data';
run;






proc iml;
use cdata; read all into xy;
x_orig = xy[,1]; y_orig = xy[,2];



start regression;
	b = inv(x`*x)*x`*y;
	yh = x*b; yb = mean(y);
	err = y-yh;
	ssr = (yh-yb)`*(yh-yb);
	sse = (y-yh)`*(y-yh);
	sst	= (y-yb)`*(y-yb);
	r2 = ssr/sst;
finish regression;


start structural_break;
	y = y;
	n = nrow(y);
	x = J(n,1,1) || x ||
		((x-x1)#((x-x1)>0)) || ((x-x2)#((x-x2)>0));
	call regression;
finish structural_break;


x1=177; x2=255;

* _______ fit original model _______;
x=x_orig; y=y_orig;
call structural_break;
e = err; 
yh_true = yh;

vis = xy || yh;
create vis from vis[colname={'x' 'y' 'yh'}]; append from vis;


* _______ get results _______;
do i=1 to 1000;
	y = yh + sample(e, nrow(e), 'replace')`;
	x = x_orig;
	call structural_break;
	res = res // (r2 || b`);
end;



* calculate confidence intervals - will not work in loop?;
call sort(res, {1});
res = res || J(nrow(res),1,res[nrow(res)*0.025,1]) || J(nrow(res),1,res[nrow(res)*0.975,1]);

call sort(res, {2});
res = res || J(nrow(res),1,res[nrow(res)*0.025,2]) || J(nrow(res),1,res[nrow(res)*0.975,2]);

call sort(res, {3});
res = res || J(nrow(res),1,res[nrow(res)*0.025,3]) || J(nrow(res),1,res[nrow(res)*0.975,3]);

call sort(res, {4});
res = res || J(nrow(res),1,res[nrow(res)*0.025,4]) || J(nrow(res),1,res[nrow(res)*0.975,4]);

call sort(res, {5});
res = res || J(nrow(res),1,res[nrow(res)*0.025,5]) || J(nrow(res),1,res[nrow(res)*0.975,5]);




create res from res[colname=
	{'r2' 'b0' 'b1' 'b2' 'b3' 'lcl_r' 'ucl_r' 'lcl_0' 'ucl_0' 'lcl_1' 'ucl_1' 'lcl_2' 'ucl_2' 'lcl_3' 'ucl_3'}];
append from res;
quit; 

proc sgplot data=vis;
	scatter y=y x=x;
	scatter y=yh x=x;
	title 'Fitted Values: Peicewise Regression';
run;



proc sgplot data=res;
	histogram r2;
	refline lcl_r / axis=x lineattrs=(color=blue);
	refline ucl_r / axis=x lineattrs=(color=blue);
	title 'Bootstrap Residuals: R2';
run;

proc sgplot data=res;
	histogram b0;
	refline lcl_0 / axis=x lineattrs=(color=orange);
	refline ucl_0 / axis=x lineattrs=(color=orange);
	title 'Bootstrap Residuals: Beta 0';
run;

proc sgplot data=res;
	histogram b1;
	refline lcl_1 / axis=x lineattrs=(color=red);
	refline ucl_1 / axis=x lineattrs=(color=red);
	title 'Bootstrap Residuals: Beta 1';
run;

proc sgplot data=res;
	histogram b2;
	refline lcl_2 / axis=x lineattrs=(color=red);
	refline ucl_2 / axis=x lineattrs=(color=red);
	title 'Bootstrap Residuals: Beta 2';
run;

proc sgplot data=res;
	histogram b3;
	refline lcl_3 / axis=x lineattrs=(color=green);
	refline ucl_3 / axis=x lineattrs=(color=green);
	title 'Bootstrap Residuals: Beta 3';
run;


	
proc sgplot data=res;
	scatter x=x y=y;
	series y=yh x=x / lineattrs=(color=blue thickness=2);
	series y=lcl x=x / lineattrs=(color=red thickness=0.25);
	series y=ucl x=x / lineattrs=(color=red thickness=0.25);
title 'Local Polynomial Regression';
run;




















* ____________ ____________________ Question 4: Structural Breaks Observational Search ____________________ ____________;

data cdata;
set '/folders/myfolders/sasuser.v94/EKT 720/Assignment 6/cdata.sas7bdat';
run;





proc iml;
use cdata; read all into xy;
x_orig = xy[,1]; y_orig = xy[,2];
x=x_orig; y=y_orig; n=nrow(xy);


* _______ metric to evaluate fit: MSE _______;
start mse;
	* check for singularity;
	if inv(x`*x) ^= 0 then do;
		b = inv(x`*x)*x`*y;
		yh = x*b;
		mse = ((y-yh)`*(y-yh))/(nrow(x)-ncol(x));
	end;
finish mse;



* _______ Model: define x _______;
start structural_design;
	x = J(n,1,1) || x_orig || 
		(x_orig-x1)#((x_orig-x1)>0) || ((x_orig-x2)#((x_orig-x2)>0));
finish structural_design;



* _______ padding: mini obs above & below _______;
p = 3;






* _______ observational search _______;
do i=p to n-2*p by 1;
	do j=i+p to n-p by 1;
		x1=x_orig[i]; x2=x_orig[j];
		call structural_design;
		call mse;
		res = res // (mse || x1 || x2);
	end;
end;

print res;



















