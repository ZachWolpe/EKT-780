libname sev "C:\Users\Mbavhalelo\Documents\Academics\BCom Hons Stats\2nd Semester\EKT 720\Buwang\ass 5";

proc logistic data=sev.model outdesign=des1;
*class d3sc d8s gmesc q55as q50b_7s nd5s ;
*model rq34_1s=d3sc d8s gmesc q55as q50b_7s nd5s ;

class  d3sc d8s gmesc nd5s q50b_7s q55as ;
  model rq34_1s = d3sc d8s gmesc nd5s q50b_7s q55as/lackfit outroc=sev.graph ;
run;

data new (keep=rq34); 
set sev.model;
rq34 = rq34_1s + 0 ; 
run;

proc iml;
use new; 
read all into y;
use des1; 
read all into x;

c=ncol(x); 
B_0=j(c,1,0);*starting point; 
tol=10##-12; 
error=1000; 
yn=1-y;*important;
*newton-raphson estimation; 
do while (error>tol); 

p=exp(x*B_0)/(1+exp(x*B_0)); 
vec_w=p#(1-p); 
w=diag(vec_w); 
B_new=B_0 + inv(x`*w*x)*x`*(yn-p);
error=max(abs(B_new-B_0));
B_0=B_new; 

end; 

nm={ 'intercept' 'd3sc 10-11' 'd3sc 12-13' 'd3s 14-15' 'd8s 1' 'd8s 2' 'd8s 3' 'gmesc High' 'gmesc Low' 'q55as 0' 'q50b_7s 0' 'nd5s 1' 'nd5s 2'};

print B_new[rowname=nm] error;

quit;
