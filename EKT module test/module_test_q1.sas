
* ____________________________________ MODULE TEST ____________________________________;

* ____________________ Question 1: Monte Carlo Integration _________________________;


proc iml;
title "Question 1: Monte Carlo Integration";
title2 "question 1.1";
start function(t);
	fx = 5160*(1 - (0.8*exp(-0.1*t)));
	return fx;
finish function;


start MC_integration(function,x1,x2,iters);

	t = do(x1,x2,0.05)`;
	y = function(t);
	
	y2 = max(y) + 20;
	y1 = 0;
	
	* question 1.1;
	
	do i=1 to iters;
		sample_y =  (y2 - y1)*rand('uniform');
		sample_x = x1 + (x2 - x1)*rand('uniform'); 
		
		* evaluate point;
		fx = function(sample_x);
		if sample_y < fx then; results = results // (sample_x || sample_y ||1);
		if sample_y > fx then; results = results // (sample_x || sample_y ||0); 
	end;
	
	total_area = (x2-x1)*(y2-y1);
	
	q1_1 = mean(results[,3])*total_area;
	print q1_1;
finish MC_integration;

call MC_integration(function,2,20,10000);


title "Question 1: Monte Carlo Integration";
title2 "question 1.2";

t = do(200,400,0.05)`;
y = function(t);
max = max(y);

tt= t||y||J(nrow(t), 1, max);


print max;

create xy from tt[colname={'x' 'y' 'max'}];
append from tt;

print "Based on the functional form it is clear that the growth rate is 
	maximized when the second terms is zero;
	The max value for y is Y=5160";

proc sgplot data=xy;
	series x=x y=y;
	refline max / axis=x;








