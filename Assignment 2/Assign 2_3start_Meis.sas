options nocenter ;
data a;
input income1 income2 house stand double prop;
income = income1 + income2;
ratio = house/stand;
cards;
  521.502    118.348    735.779     920.53      0     1215.86
   14.116    457.801    413.522     690.15      0      917.67
  308.237    205.341    567.238     903.11      1     1201.16
  449.589    157.470    496.226     659.05      1     1099.73
   47.286    555.871    414.292     769.92      0     1077.25
   12.702    400.744    283.223     539.62      1      973.28
  303.539    360.630    671.121     934.58      0     1206.01
  325.548    369.610    284.473     707.85      0     1071.03
  328.079     17.192    492.552     699.23      0      834.60
  479.735     34.212    767.408    1097.32      0     1102.11
   70.381    319.148    373.140     760.19      0      774.12
  232.232    255.517    238.515     577.39      0      828.64
   56.125    326.705    589.865     930.42      0     1020.10
  510.569     36.773    461.059     920.65      0     1044.39
   15.890    353.851    345.385     655.05      1      932.65
  298.906    126.398    531.592    1093.24      0     1036.68
  280.401    105.089    497.296     727.87      0      867.00
  188.411    419.229    383.097     903.32      0     1114.49
   11.004    462.602    351.969     575.44      0      833.46
  357.141    324.151    507.147     686.02      0     1146.86
   74.029    403.535    372.881     520.79      0      836.79
  112.752    195.755    550.987    1048.71      0     1023.95
  189.496    273.100    400.458     550.31      0      834.02
  283.516    395.697    445.404     600.35      0     1064.84
  255.701    154.743    535.123    1078.51      0     1075.30
  ;
run ;
proc print data = a;
run;

proc iml ;

use a ;
read all var{prop} into y ;
read all var{income} into x2;
read all var{ratio} into x4 ;
read all var{stand} into x3;
read all var{double} into x5 ;
print y, x1;

n=nrow(y);
x = J(n,1,1)||x2||x3||x4||x5;
print x;

* Question a;
b = inv(x`*x)*x`*y;
print b;

* Question b;
ybar = sum(y)/n;
yhat = x*inv(x`*x)*x`*y;
ssr = (yhat-ybar)`*(yhat-ybar);
ssto = (y-ybar)`*(y-ybar);
R_2 = ssr/ssto;
print R_2;

* Question b;
sse = (y-yhat)`*(y-yhat);
adR_2 = 1-(29/25)*(sse/ssto);
print adR_2;

* Question c;
F_test_stat = (ssr/4)/(sse/25);
print F_test_stat;
crit_val = finv(.95,4,25);
print crit_val;

* Question d;
print sse;
mse = sse/25;
var_b = vecdiag(mse*inv(x`*x));
print var_b;
std_b3 = sqrt(var_b[3,]);
print std_b3;

t_value = 0.52021/std_b3;
print t_value;
t_crit_value = tinv(0.95,25);
print t_crit_value;

* Question e;
std_b = sqrt(var_b);
CLb_L = b-tinv(0.90,25)*std_b; 
CLb_U = b+tinv(0.90,25)*std_b;
print 'lower' CLb_L, 'upper' CLb_U;


new_x = x[,1]||x[,3]||x[,5];
print new_x;
yhat_n = new_x*inv(new_x`*new_x)*new_x`*y;
ssr_n = (yhat_n-ybar)`*(yhat_n-ybar);
ssto_n = (y-ybar)`*(y-ybar);
R_2_n = ssr_n/ssto_n;
print R_2_n;

F_n = (R_2-R_2_n)/((1-R_2)/27);
print F_n;


F_pvalue_n = 1-probf(F_n,2,27);
print F_pvalue_n;












quit ;
