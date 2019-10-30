options ls=72 ps=1000 nodate pageno=1 ;
quit ;
dm 'odsresults;clear';


title "Logistic regression" ;
title2 "A SAS overview, including goodness of fit" ;


data a (keep = age chd chd2 agec)  ;
input ID	AGE	CHD @@ ;
v=rannor(3466) ;
chd2=0 ;
if v > 0 then chd2=1 ;
                  agec="c:51-- ";
if age <= 50 then agec="b:36-50" ;
if age <= 35 then agec="a:-35  " ;
cards ;
1	20	0 2	23	0 3	24	0 4	25	0
5	25	1 6	26	0 7	26	0 8	28	0
9	28	0 10	29	0 11	30	0 12	30	0
13	30	0 14	30	0 15	30	0 16	30	1
17	32	0 18	32	0 19	33	0 20	33	0
21	34	0 22	34	0 23	34	1 24	34	0
25	34	0 26	35	0 27	35	0 28	36	0
29	36	1 30	36	0 31	37	0 32	37	1
33	37	0 34	38	0 35	38	0 36	39	0
37	39	1 38	40	0 39	40	1 40	41	0
41	41	0 42	42	0 43	42	0 44	42	0
45	42	1 46	43	0 47	43	0 48	43	1
49	44	0 50	44	0 51	44	1 52	44	1
53	45	0 54	45	1 55	46	0 56	46	1
57	47	0 58	47	0 59	47	1 60	48	0
61	48	1 62	48	1 63	49	0 64	49	0
65	49	1 66	50	0 67	50	1 68	51	0
69	52	0 70	52	1 71	53	1 72	53	1
73	54	1 74	55	0 75	55	1 76	55	1
77	56	1 78	56	1 79	56	1 80	57	0
81	57	0 82	57	1 83	57	1 84	57	1
85	57	1 86	58	0 87	58	1 88	58	1
89	59	1 90	59	1 91	60	0 92	60	1
93	61	1 94	62	1 95	62	1 96	63	1
97	64	0 98	64	1 99	65	1 100	69	1
;
run ;

proc means data=a ; 
var chd age  ;
run ;
proc freq data=a ;
tables  chd chd2 ;
run ;

/*data a ;*/
/* set a ;*/
/* chd=chd2 ;*/
/*run ;*/

proc iml ;
 use a ;
  read all into xy ;

 n= nrow(xy) ;
 x = J(n,1,1) || xy[,1];
 y = xy[,2] ;

 bhi = {0 , 0} ;
 bho=bhi ;
 diff=11111 ;

do i = 1 to 15 while (diff>0.00001) ;
	 lo = x*bho ;
	 o = exp(lo);
	 p = o/(1+o) ;
	 w = diag(p#(1-p)) ;

	 bhn = bho + inv(x`*w*x)*x`*(y-p) ;
	 print i bho bhn ;

	 diff= max(abs(bhn-bho)) ;
	 bho=bhn ;
  end ;

 quit ;

 * PREDICTED PROBABILITIES FOR CHD;
proc logistic data = a ;
model chd(event="1") = age / outroc=roc1 lackfit; 
*model chd2(event="1") = age / outroc=roc1 lackfit; 
output out=b p=pred ;
run ;



*--------------------------------------------------------------------------------------------------------------------------------------------------;
proc means data=b min p10 p20  p30 p40 p50 p60 p70 p80 max ;
 var pred ;
 output out=cc min= p10= p20=  p30= p40= p50= p60= p70= p80= max= /autoname;
run ;

proc print data = b;
run;

* TAKE PREDICTED PROB AND GROUP THEM CATEGORICALLY;
data dd (drop = _TYPE_ _FREQ_);
 set b ;
 if _N_ =1 then set cc ; 
 if pred <= pred_max then pred_c = "p80 - max" ;
 if pred <= pred_p80 then pred_c = "p60 - p80" ;
 if pred <= pred_p60 then pred_c = "p40 - p60" ;
 if pred <= pred_p40 then pred_c = "p20 - p40" ;
 if pred <= pred_p20 then pred_c = "min - p20" ;
run ;

proc freq data=dd;
 tables pred_c*chd ;
run ;

* CUM PERCENTAGE FROM PROC FREQ FOR CHD=0;
proc freq data=dd;
 tables pred_c*chd / out=ee1;
 where chd=0 ;
run ;

proc print data = ee1;
run;
****************** now make cumper0**********************;
data ee1 (keep = pred_c freq0 perc0 cumper0);
 set ee1 ;
 freq0 = count ;
 perc0 = percent ;
 retain cumper0;
 cumper0+perc0;
run ;

proc print data = ee1;
run;


* CUM PERCENTAGE FROM PROC FREQ FOR CHD=1;
proc freq data=dd;
 tables pred_c*chd / out=ee2;
 where chd=1 ;
run ;
****************** now make cumper1**********************;
data ee2 (keep = pred_c freq1 perc1 cumper1);
 set ee2 ;
 freq1 = count ;
 perc1 = percent ;
 retain cumper1;
 cumper1+perc1;
run ;





* CUM MUST START FROM 0 SO MAKE AN EXTRA ROW WITH ZEROS TO ADD TO THE REST;
data ff ;
 pred_c = "min      " ;
 cumper0 = 0 ;
 cumper1 = 0 ;
run ;

proc print data=ff;
run;

* FINALLY CUMPER1 AND CUMPER0 THAT START FROM 0;
data gg ;
 merge ff ee1 ee2 ;
 by pred_c ;
run ;

proc print data = gg;
run;

* GET THE DIFF BETWEEN CHD=1 AND CHD=0 FOR EACH PERCENTILE GROUP ---> AND MAKE IT A PROBABILITY;
data gg ; 
 set gg ;
 ks = abs(cumper1-cumper0);
 ksr = ks/100 ;
run ;

proc print data=gg ;
run ;

proc means data=gg max;
 var ks ksr;
run ;






*------------------------------------------------------------------------------------------------------------------------------------;

* MAKING A THRESHOLD OF 0.5 to divide predicted prob for each obs into chd=0/1;
********* b is proc logistic's output*********;
data b ; 
set b ;
chd_p = 0 ;
if pred > 0.5 then chd_p= 1 ;
run ;

* CONFUSION MATRIX WITH ACTUAL VS PREDICTED CHD;
proc freq data=b ;
 tables chd*chd_p ;
run ;







* PLOTS * ------------------------------------------------------------------------------------------;
symbol1 interpol=join width=4 color=blue value=dot  height=1;
symbol2 interpol=join width=3 color=red  value=dot  height=1;


proc gplot data=gg;
 plot ( cumper0 cumper1)*cumper1 / overlay href=(4.6511627907 16.279069767 37.209302326
       72.093023256) chref=purple ;
run ;


proc gplot data=gg;
 plot (cumper1 cumper0 )*cumper1 / overlay href=(4.6511627907 16.279069767 37.209302326
       72.093023256) chref=purple areas=3;
run ;

proc gplot data=gg;
 plot (cumper1 cumper0)*pred_c / overlay ;
run ;


proc print data=gg ;
run ;

data hh (keep = cumper0 cumper1); 
set gg ;
run ;

proc iml ;
 use hh ;
 read all into gini ;
 gini = gini/100 ;
  

 print gini ;
 smp = 10000;
 *s1 = J(smp,1,1241240) ;
 
 do i = 2 to nrow(gini) ;
  slope = (gini[i,1] - gini[i-1,1])  / (gini[i,2] - gini[i-1,2]  ) ;
  int = gini[i,1] - slope*gini[i,2] ;
  print i slope int;
     
    n1=round(smp*(gini[i,2] - gini[i-1,2] )) ;
    s1 = J(n1,1,2552350) ;
    y=ranuni(s1) ;
    x=gini[i-1,2]+ ranuni(s1)*(gini[i,2] - gini[i-1,2]  ) ;
	fv = int + slope*x ;  
    check = (y < fv) ;
	check2 = (y > x)*2;    
	check3= check+check2 ;
*	print i y x fv check;
  res = J(n1,1,i) || x || y || fv || check || check3;
  fsim = fsim // res ; ;
 end ;

 *print fsim ;

 nm={"g" "x" "y" "fv" "check" "check3"} ;
 create gini1 from fsim[colname=nm] ;
 append from fsim ;
 close gini1 ;

 ginicoef = ((fsim[,5]/nrow(fsim))[+]-0.5)/0.5 ;
 print ginicoef ;

quit ;

symbol1 interpol=none width=4 color=blue value=dot  height=.7;
symbol2 interpol=line width=3 color=red  value=dot  height=.7;

proc sort data=gini1 ;
 by g x ;
run ;

proc gplot data=gini1 ;
 plot (y fv)*x  fv*fv  / overlay ;
 *where g<=6 ;
run ;


symbol1 interpol=none width=4 color=blue value=dot  height=.7;
symbol2 interpol=none width=3 color=red  value=dot  height=.7;
symbol3 interpol=none width=3 color=green  value=dot  height=.7;

proc gplot data=gini1 ;
 plot  y*x=check3 / areas=3;
run ;
quit ;

proc means data=gini1 ;
var check ;
run ;



quit ;
