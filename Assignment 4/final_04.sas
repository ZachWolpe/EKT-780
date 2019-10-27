data xy_data; 
input x y;
datalines;
10 44
11 34
12 19
13 14
14 24
15 29
16 26
17 22
18 6
19 18
20 19
21 23
22 31
23 31
24 38
25 39
26 41
27 42
28 36
29 28
30 24
31 18
32 8
33 25
34 29
35 30
36 30
37 34
38 32
39 36
;
run;

proc sgplot data=xy_data;
	scatter y=y x=x;
run;


options center;
* _____________________ Local Polynomial Regression _____________________;
proc iml;
use xy_data; read all into xy;
n = nrow(xy);


d = distance(xy[,1]);
k = 5;



do i=1 to n;
* for each datapoint;

	focal = xy[i,];
	
	col = d[,i] || xy;
	call sort(col, {1});
	col = col[2:k,];

	* local polynomial regression;
	x = J(nrow(col),1,1) || col[,2]; y = col[,3];	
	b = inv(x`*x)*x`*y;
	yh = b[1] + focal[1]*b[2]; yh_all = x*b;

	kk = ncol(x); nn=nrow(x);
	
	mse = sum((y-yh_all)##2)/(nn-kk);
	
	*sebh = sqrt(mse#inv(x`*x));
	lcl = yh - tinv(0.975, nn-kk)*sqrt(mse);
	ucl = yh + tinv(0.975, nn-kk)*sqrt(mse);
	
	yh = yh || lcl || ucl;
	reg = reg // yh;

end;

print reg;
res = xy || reg;
create res from res[colname={'x' 'y' 'yh' 'lcl' 'ucl'}];
append from res;
quit;



proc sgplot data=res;
	scatter x=x y=y;
	series y=yh x=x / lineattrs=(color=blue thickness=2);
	series y=lcl x=x / lineattrs=(color=red thickness=0.25);
	series y=ucl x=x / lineattrs=(color=red thickness=0.25);
title 'Local Polynomial Regression';
run;





