options ls=72  nodate pageno=1 ;
quit ;
dm 'odsresults;clear';
title ;

*Agec C: 51+
      B: 36-50
      A: Up to 35 ; 

data a  ;
input ID chd agec$ @@;
d1=-1 ; d2=-1 ;
if agec ="A" then do ; d1=1 ; d2=0 ; end ;
if agec ="B" then do ; d2=1 ; d1=0 ; end ;
cards;
1	0	A	26	0	A	51	1	B	76	1	C
2	0	A	27	0	A	52	1	B	77	1	C
3	0	A	28	0	B	53	0	B	78	1	C
4	0	A	29	1	B	54	1	B	79	1	C
5	1	A	30	0	B	55	0	B	80	0	C
6	0	A	31	0	B	56	1	B	81	0	C
7	0	A	32	1	B	57	0	B	82	1	C
8	0	A	33	0	B	58	0	B	83	1	C
9	0	A	34	0	B	59	1	B	84	1	C
10	0	A	35	0	B	60	0	B	85	1	C
11	0	A	36	0	B	61	1	B	86	0	C
12	0	A	37	1	B	62	1	B	87	1	C
13	0	A	38	0	B	63	0	B	88	1	C
14	0	A	39	1	B	64	0	B	89	1	C
15	0	A	40	0	B	65	1	B	90	1	C
16	1	A	41	0	B	66	0	B	91	0	C
17	0	A	42	0	B	67	1	B	92	1	C
18	0	A	43	0	B	68	0	C	93	1	C
19	0	A	44	0	B	69	0	C	94	1	C
20	0	A	45	1	B	70	1	C	95	1	C
21	0	A	46	0	B	71	1	C	96	1	C
22	0	A	47	0	B	72	1	C	97	0	C
23	1	A	48	1	B	73	1	C	98	1	C
24	0	A	49	0	B	74	0	C	99	1	C
25	0	A	50	0	B	75	1	C	100	1	C
;
run ;





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



data a;
set a;
b1=-1; b2=-1;
if age <= 30 then b1=1; 
if age <= 30 then b2=0;
if age <= 45 and age > 30 then b1=0; 
if age <= 45 and age > 30 then b2=1;
run;









proc print data=a;
run ;

proc freq data = a;
tables agec*chd / chisq expected;
run;

proc means data=a ;
var chd   ;
run ;

proc logistic data = a outdesign=des;
class agec ;
model chd(event="1") = agec / outroc=roc1 lackfit; 
output out=b p=pred ;
run ;



 * ---------------------  Dummy Variable Logistic Regression --------------------- *;
title "Logistic Regression";

* ______ create categorical variable for PROC IML ______;


proc iml;
use a;
read all into xy;

n=nrow(xy);
y = xy[,1];
b1 = xy[,3]; b2 = xy[,4];
x = J(n,1,1)||b1||b2;
* print x y;


* initialize Beta=0 ==> p=1;
bho={0,0,0};
diff=999;


start logistic_regression;
	do i=1 to 20 while (diff>0.00001);
	lo = x*bho ;
	 o = exp(lo);
	 p = o/(1+o);
	 w = diag(p#(1-p));

	 bhn = bho + inv(x`*w*x)*x`*(y-p) ;
	 print i bho bhn ;

	 diff= max(abs(bhn-bho)) ;
	 bho=bhn ;
	end;
finish logistic_regression;


call logistic_regression;
print bho;




 
