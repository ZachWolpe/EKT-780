
************** ******************** Question 1: Logic - Monte Carlo Integration ******************** **************;
title 'Question 1';

proc iml;

start function_x(x);
	fx = 10 + 6.32972*x - 1.72728*(x**2) + 0.2017*(x**3) 
			- 0.00996*(x**4) + 0.00017*(x**5);
	return fx;
finish function_x;



start MC_integration(function_x, iterations, x_low, x_up);
	
	* approximate Y range (domain space);
	range = do(x_low, x_up, 0.5)`;
	
	do i=1 to nrow(range);
		yt = yt // function_x(range[i]);
	end;
	
	y_up = max(yt) + 10;
	y_low = 0;
	
	 
	do i=1 to iterations;
		sample_y = y_low + (y_up - y_low)*rand('uniform');
		sample_x = x_low + (x_up - x_low)*rand('uniform');

		* evaluate point;
		fx = function_x(sample_x);
		if sample_y < fx then; results = results // 1;
		if sample_y > fx then; results = results // 0;
	end;
	
	total_area = (y_up - y_low)*(x_up - x_low);
	results = mean(results);
	area_under_curve = results*total_area;
	
	return area_under_curve;	
finish MC_integration;





n = 10000;
title ' Question 1.1';
q11 = MC_integration(function_x, n, 3, 8);
print 'Area Under Curve: ' (q11);

title 'Question 1.2'; 
q12 = MC_integration(function_x, n, 1, 10);
print 'Area Under Curve: ' (q12);

title 'Question 1.3';
q13 = MC_integration(function_x, n, 0, 20);
print 'Area Under Curve: ' (q13);








* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ASSIGNMENT 6 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ;    

************** ******************** Question 1: Logic - Monte Carlo Integration ******************** **************;

title 'Question 1';
proc iml;

start function_x(x);
	fx = 10 + 6.32972*x - 1.72728*(x**2) + 0.2017*(x**3) 
		- 0.00996*(x**4) + 0.00017*(x**5);
	return fx;	
finish function_x;



start MC_integration(function_x, n, x_lower_limit, x_upper_limit);


	t = do(x_lower_limit,x_upper_limit,.1)`;
	do i=1 to nrow(t);
		yt = yt // function_x(t[i]);
	end;

	approx_height_u = round(max(yt))+3;
	approx_height_l = round(min(yt))-3;
	approx_height_l = 0;
	
	results = J(n,3,999);
	do i=1 to n;
		x = rand('integer', round(x_lower_limit), round(x_upper_limit));
		y = rand('integer', approx_height_l, approx_height_u);
		y_bound = function_x(x);
		results[i,1] = x; 
		results[i,2] = y;
		if y > y_bound then do; 
			results[i,3] = 0; 
			end;
		else do; 
			results[i,3] = 1; 
			end;
	end;
	
	total_area = (approx_height_u - approx_height_l) * 
		(round(x_upper_limit) - round(x_lower_limit));
 	area_under_curve = mean(results[,3])*total_area;
	return area_under_curve;
finish MC_integration;


n = 10000;
* Question 1.1;
q11 = MC_integration(function_x, n, 3, 8);
print 'Area Under Curve: ' (q11);

* Question 1.2; 
q12 = MC_integration(function_x, n, 1, 10);
print 'Area Under Curve: ' (q12);

* Question 1.3;
q13 = MC_integration(function_x, n, 0, 20);
print 'Area Under Curve: ' (q13);
















