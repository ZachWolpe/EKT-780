

* ________________________________________ Assignment 8 ________________________________________;

dm 'odsresults;clear'; options ls=72 nodate pageno=1 center ; 



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


proc nlin method=newton;
	parms th1=1 th2=1;
	model y = th1 + th2*x/(exp(x));
run ;



* ______________________________________ Question 1 ______________________________________;

data a (keep=x y);
	do i=1 to 20;
		x=i;
		y=30 - 80*x/(exp(x))+rannor(11567)*1;
		output;
	end;
run;

proc iml;
title "Algorithm 1";
title2 'Non linear estimation: Newton Raphson';
use a; read all into xy;
x = xy[,1]; y = xy[,2];
n = nrow(xy);


theta={10,1000};
diff = 100;
pp=1;

do i=1 to 20 while (diff>0.0001);
	
	yh = theta[1] + (theta[2]#x)/exp(x);
	q1 = -2*(y-yh)[+];
	q2 = -2*(((y-yh)#X)/exp(x))[+];
	
	* print q1 q2;
	
	h11 = 2*n;
	h22 = 2*(X##2#exp(-2*x))[+];
	h12 = 2*(x#exp(-X))[+];
	h21 = 2*(x#exp(-X))[+];
	
	* print h11 h12 h21 h22;
	
	t = q1//q2;
	b = (h11 || h12) // (h21 || h22);
	theta_n = theta - inv(b)*t;
	
	if pp=1 then; print i theta_n theta;
	diff = max(abs(theta_n - theta));
	theta = theta_n;

end;

* fit model;
yh = theta[1] + (theta[2]#X)/exp(X);
results = xy || yh;
print results;

create results from results[colname={'x' 'y' 'yh'}];
append from results;

proc sgplot data=results;
	series y=yh x=x;
	scatter y=y x=x;
run;




* ______________________________________ Question 2 ______________________________________;

data a (keep=x y);
	s=0.5;
	do i=1 to 100;
		x=i;
		y=sqrt((2.35*X+5.89)) + rannor(5675)*s;
		output;
	end;
title; 
title2;


proc iml;
title "Algorithm 2";
title2 'Non linear estimation: Newton Raphson';
use a; read all into xy;
x = xy[,1]; y=xy[,2]; n=nrow(xy);

theta = {1,1};
diff = 100;
pp=1;



do i=1 to 20 while (diff>0.0001);

	t = X#theta[1]+theta[2];
	
	q1 = (-(X#(y-sqrt(t)))/sqrt(t))[+];
	q2 = (-(y-sqrt(t))/sqrt(t))[+];
	
	h11 = ((X##2#(y-sqrt(t)) / (2#(t)##(3/2))) + (X##2/(2#t)))[+];
	h12 = ((X#(y-sqrt(t)) / (2#(t)##(3/2))) + (X/(2#t)))[+];
	h21 = ((X#(y-sqrt(t)) / (2#(t)##(3/2))) + (X/(2#t)))[+];
	h22 = ((y-sqrt(t) / (2#(t)##(3/2))) + (1/(2#t)))[+];
	
	q = q1//q2; h=(h11||h12)//(h21||h22);
	theta_n = theta -inv(h)*q;
	if pp=1 then; print i theta_n theta;
	diff = max(abs(theta_n-theta));
	theta = theta_n;
	
end;

* fit model;
yh = sqrt(theta_n[2] + theta_n[1]#X);
results = xy || yh;

create results from results[colname={'x' 'y' 'yh'}];
append from results;

proc sgplot data=results;
	series y=yh x=x;
	scatter y=y x=x;
run;




















