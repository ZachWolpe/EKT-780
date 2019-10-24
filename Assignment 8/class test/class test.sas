libname lib '/folders/myfolders/sasuser.v94/EKT 720/Assignment 8/class test';

proc print data=lib.classnlin (obs=10);

proc sgplot data=lib.classnlin;
	scatter y=y x=y;
run;

data a;
set lib.classnlin;
run;

proc iml;
title "Algorithm 1";
title2 'Non linear estimation: Newton Raphson';
use a; read all into xy;
x = xy[,1]; y = xy[,2];
n = nrow(xy);
n = nrow(xy);
y = xy[,2];
x = xy[,1];


start sse(parms) global (x,y);
	bh=parms;
	yh = 1 - exp(bh[1]+bh[2]#X);
	sse=(y-yh)`*(y-yh);
	return (sse);
finish;



parms = {-0.1,-5.5} ;
sse= sse(parms);

print "Initial values " parms sse;
opt = {0,    /* find maximum of function  */
      0};
      /* 4};   /* print a LOT of output     */
      
call nlpnra(rc, parm_estimates, "sse", parms, opt);
print parm_estimates ;



do i=1 to n;
	yh = yh // 1 - exp(parm_estimates[1] + parm_estimates[2]*x[i]);
end;

res = xy || yh;

create res from res[colname={"x" 'y' 'yh'}];
append from res;

* _____________ Question 6 ______________;


print 'PREDICTION QUESTION 6';
q6_a =  1 - exp(parm_estimates[1] + parm_estimates[2]*3.5);
q6_b =  1 - exp(parm_estimates[1] + parm_estimates[2]*5/100);
print q6_a q6_b;


* _____________ Question 7 ______________;


title 'Bootstrap';


do i=1 to 1000;

	xy_set = sample(xy, n, 'replace');
	x = xy_set[,1];
	y = xy_set[,2];
	
	parms = {-0.1,-5.5} ;
	sse= sse(parms);
	
	*print "Initial values " parms sse;
	opt = {0,    /* find maximum of function  */
	       4};   /* print a LOT of output     */
	
	call nlpnra(rc, parm_estimates, "sse", parms, opt);
	*print parm_estimates ;

	

	dd = dd // parm_estimates;

end;

create dd from dd[colname={'b0' 'b1'}];
create from dd;












proc sgplot data=res;
	scatter y=y x=x;
	scatter y=yh x=x;
run;




* _________________ Question 4 _________________;


* CANNOT BE COMPUTED ;


















