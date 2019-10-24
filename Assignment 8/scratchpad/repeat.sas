

title2 'Non linear estimation: Newton Raphson';
data a (keep=x y);
	do i=1 to 20;
		x=i;
		y=30 - 80*x/(exp(x))+rannor(11567)*1;
		output;
	end;
run;
proc print data=a (obs=5);


proc sgplot data=a;
	scatter y=y x=x;
run;




proc iml;
use a; read all into xy; x=xy[,1]; y=xy[,2]; n = nrow(xy);

theta = {10,1000};
diff = 100;

do i=1 to 15 while (diff>0.0001);
	
	
	q1 = (2*(X*theta[2]/(exp(x)) - y + theta[1]))[+];
	q2 = (-2*X/(exp(x))#(-X*theta[2]/(exp(x)) + y - theta[1]))[+];
	
	h11 = 2*n;
	h12 = (2*X/(exp(x)))[+];
	h21 = h12;
	h22 = (2#(X##2)/(exp(2#x)))[+];
	
	theta_n = theta - inv((h11||h12)//(h21||h22))*(q1//q2);
	diff = max(abs(theta_n - theta));
	print i theta theta_n;
	theta = theta_n;
end;

* predict;
yh = theta[1] + theta[2]#X#exp(-x);

res =  x || y || yh;
create res from res[colname={'x' 'y' 'yh'}];
append from res;

proc sgplot data=res;
	scatter y=y x=x;
	series y=yh x=x;
run;



* _____________________ Question 1: NLPNRA _____________________;


proc iml;
use a; read all into xy; x=xy[,1]; y=xy[,2]; n = nrow(xy);

theta = {10,1000};
parms = theta;

start sse(parms);
	yh = parms[1] + parms[2]*X*exp(-x);
	sse = (y-yh)`*(y-yh);
   return sse;
finish;

parms = {2,6} ;
sse= sse(parms);

opt = {0, 		/* find maximum of function   */
	   4}; 		/* print a LOT of output      */
	   
call nlpnra(rc, parm_estimates, 'sse', parms, opt);

print parm_estimates ;






* _____________________ Question 2 _____________________;


data a (keep = x y);
s=0.5;
do i=1 to 100;
	x=i;
	y=sqrt(2.35*x + 5.89) + rannor(5675)*s;
	output;
end;

proc sgplot data=a;
	scatter y=y x=x;





proc iml;
use a; read all into xy; x=xy[,1]; y=xy[,2]; n=nrow(xy);

theta = {1,1};
diff = 100;

do i=1 to 10 while (diff>0.0001);
	t = sqrt(x#theta[1] + theta[2]);
	tt = x#theta[1] + theta[2];	
	
	q1 = -((x#(y-t))/t)[+];
	q2 = -((y-t)/t)[+];
	
	h11 = sum(((X##2)#(y-t)/(2*(tt##(3/2)))) + X##2/(2#tt));
	h12 = sum((X#(y-t)/(2*(tt##(3/2)))) + X/(2#tt));
	h21 = h12;
	h22 = ((y-t)/2#(tt##(3/2)) + 1/(2#tt))[+];
	
	Q = q1//q2; H = (h11||h12)//(h21||h22);
	theta_n = theta - inv(H)*Q;
	print i theta theta_n;
	diff = max(abs(theta_n - theta));
	theta = theta_n;

end;

* predict;
yh = sqrt(X*theta[1] + theta[2]);

res =  x || y || yh;
create res from res[colname={'x' 'y' 'yh'}];
append from res;

proc sgplot data=res;
	scatter y=y x=x;
	series y=yh x=x;
run;





* _____________________ Bootstrap _____________________;

proc iml;
use a; read all into xy; x_orig=xy[,1]; y_orig=xy[,2]; n=nrow(xy);


do j=1 to 1000;

	* ______________ bootstrap Sample ______________;
	x = sample(x_orig, n, 'replace')`;
	y = sample(y_orig, n, 'replace')`;
	* ______________ bootstrap Sample ______________;
	
	theta = {1,1};
	diff = 100;
	
	do i=1 to 10 while (diff>0.0001);
		t = sqrt(x#theta[1] + theta[2]);
		tt = x#theta[1] + theta[2];	
		
		q1 = -((x#(y-t))/t)[+];
		q2 = -((y-t)/t)[+];
		
		h11 = sum(((X##2)#(y-t)/(2*(tt##(3/2)))) + X##2/(2#tt));
		h12 = sum((X#(y-t)/(2*(tt##(3/2)))) + X/(2#tt));
		h21 = h12;
		h22 = ((y-t)/2#(tt##(3/2)) + 1/(2#tt))[+];
		
		Q = q1//q2; H = (h11||h12)//(h21||h22);
		theta_n = theta - inv(H)*Q;
		*print i theta theta_n;
		diff = max(abs(theta_n - theta));
		theta = theta_n;
	
	end;
	* predict;
	yh = yh // theta`;
end;


create yh from yh[colname={'theta_1' 'theta_2'}];
append from yh;




proc sgplot data=yh;
	histogram theta_1;
run;

proc sgplot data=yh;
	histogram theta_2;
run;


proc iml;
use yh; read all into yh;
theta_1 = yh[,1]; theta_2 = yh[,2]; n=nrow(yh);
call sort(theta_1);
call sort(theta_2);

q025 = round(0.025*n);
q975 = round(0.975*n);

r_1 = theta_1 || J(n,1,theta_1[q025]) || J(n,1,theta_1[q975]);
r_2 = theta_2 || J(n,1,theta_2[q025]) || J(n,1,theta_2[q975]);

create r_1 from r_1[colname={'theta_1' 'lower' 'upper'}];
append from r_1;

create r_2 from r_2[colname={'theta_2' 'lower' 'upper'}];
append from r_2;


proc sgplot data=r_1;
	histogram theta_1;
	refline lower / axis=x;
	refline upper / axis=x;
title "Bootstrap Distribution";
run;









