
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

options ls=72 nocenter ;

/* proc plot data=xy;
	plot y*x;
run;
*/

proc sgplot data=xy_data;
	scatter y=y x=x;
run;


proc iml;
use xy_data;
read all into xy;
x1 = xy[,1];
y1 = xy[,2];
n = nrow(xy);
m = 7;
print n m;



do i=1 to n;
	* select each focal point;
	f_point = x1[i,];
	
	* compute distances;
	x0 = x1 - f_point;
	* print x0;

	axy = abs(x0) || x0 || xy;
	* print axy;

	* Kernel Regression;
	x = j(m,1,1) || axy[1:m, 2];
	y = axy[1:m, 4];
	
	bh = inv(x`*x)*x`*y;
	yh_fp = bh[1,1]; 		* centered thus intercept is sufficient;
	
	* print x y ;
	* print i bh  yh_fp ;
	
	yh = yh // yh_fp;

 end ;
 * print yh;
 
 
 res = xy || yh;
 nm={'x' 'y' 'yh'};
 print res[colname=nm];
 
create npres from res[colname=nm];
	append from res;
quit;



symbol1 interpol=none 
		width=4
		color=aquamarine
		value=dot
		height=2;
	
symbol2 interpol=none 
		width=3
		color=firebrick
		value=dot
		height=.3;
		


/* proc gplot data=npres;
	plot (y yh)*x / overlay;
run;
*/
proc sgplot data=npres;
	scatter y=y x=x;
	series y=yh x=x;
run;

quit;






