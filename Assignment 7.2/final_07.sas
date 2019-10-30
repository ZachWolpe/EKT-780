* ________________________________________ LOGISTIC REGRESSION ________________________________________;


data a (keep = age chd )  ;
input nr age chd @@;
cards ;
1 20 0 2  23 0 3  24 0 4  25 0
5 25 1 6  26 0 7  26 0 8  28 0
9 28 0 10 29 0 11 30 0 12 30 0
13 30 0 14 30 0 15 30 0 16 30 1
17 32 0 18 32 0 19 33 0 20 33 0
21 34 0 22 34 0 23 34 1 24 34 0
25 34 0 26 35 0 27 35 0 28 36 0
29 36 1 30 36 0 31 37 0 32 37 1
33 37 0 34 38 0 35 38 0 36 39 0
37 39 1 38 40 0 39 40 1 40 41 0
41 41 0 42 42 0 43 42 0 44 42 0
45 42 1 46 43 0 47 43 0 48 43 1
49 44 0 50 44 0 51 44 1 52 44 1
53 45 0 54 45 1 55 46 0 56 46 1
57 47 0 58 47 0 59 47 1 60 48 0
61 48 1 62 48 1 63 49 0 64 49 0
65 49 1 66 50 0 67 50 1 68 51 0
69 52 0 70 52 1 71 53 1 72 53 1
73 54 1 74 55 0 75 55 1 76 55 1
77 56 1 78 56 1 79 56 1 80 57 0
81 57 0 82 57 1 83 57 1 84 57 1
85 57 1 86 58 0 87 58 1 88 58 1
89 59 1 90 59 1 91 60 0 92 60 1
93 61 1 94 62 1 95 62 1 96 63 1
97 64 0 98 64 1 99 65 1 100 69 1
;
run ;

proc means data=a;

 
 
proc logistic data=a;
	model chd(event='1') = age / outroc=rocl lackfit;
	output out=n p=pred;
run;



* _________________ create categorical variable for PROC LOGISTIC _________________;
data aa;
set a;
group = 0;
if age <= 30 then group=1;
if age <= 45 and age > 30 then group=2;
run;


proc logistic data=aa;
	class group;
	model chd(event='1') = group / lackfit;
	title 'Categorical Independent Variable';
run;


* _________________ create categorical variable for PROC IML _________________;
data aa;
set a;
b1=-1; b2=-1;
if age <= 30 then b1=1; 
if age <= 30 then b2=0;
if age <= 45 and age > 30 then b1=0; 
if age <= 45 and age > 30 then b2=1;
run;

proc print data=aa;



*** ---------------------  Dummy Variable Logistic Regression --------------------- ***;
title "Logistic Regression";

proc iml;
use aa;
read all into xy;
n=nrow(xy);
y = xy[,2];
b1 = xy[,3]; b2 = xy[,4];
x = J(n,1,1)||b1||b2;




* initialize Beta=0 ==> p=1;
bo={0,0,0};
diff=999;




start logistic_regression;
	do i=1 to 20 while (diff>0.00001);
		lo = x*bo;
		o = exp(lo);
		p = o/(1+o);
		w = diag(p#(1-p));
		
		bn = bo + inv(x`*w*x)*x`*(y-p);
		print i bo bn;
		diff = max(abs(bn-bo));
		bo = bn;
	end;
finish logistic_regression;


call logistic_regression;
print bn;
quit;











* ______________________________  ______________________________;
























