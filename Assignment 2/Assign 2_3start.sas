options nocenter ;
data a;
input income1 income2 house stand double prop;
* create new variables;
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
  408.952    119.757    650.882     950.26      0     1044.39
  114.999    253.868    439.853     849.32      0      831.60
  200.932    141.234    400.907     571.64      0      773.70
  276.907    350.366    554.191     948.33      1     1273.24
  271.076    109.235    734.862     970.72      1     1124.66
  357.141    324.151    507.147     686.02      0     1146.86
   74.029    403.535    372.881     520.79      0      836.79
  112.752    195.755    550.987    1048.71      0     1023.95
  189.496    273.100    400.458     550.31      0      834.02
  283.516    395.697    445.404     600.35      0     1064.84
  255.701    154.743    535.123    1078.51      0     1075.30
  ;
run ;


* proc print data=a;
* run;



proc iml ;
use a;
read all var{prop} into y;
read all var{income} into x2;
read all var{stand} into x3;
read all var{ratio} into x4;
read all var{double} into x5;
n = nrow(y);
X = J(n,1,1)||x2||x3||x4||x5; * design matrix;

* Question A;
b = inv(X`*X)*X`*y;
print b;

* Question B;
ybar = (y`*J(n,1,1))/n;
yhat = X*b;
ssr = (yhat-ybar)`*(yhat-ybar); * uncorrected;
ssto = (y-ybar)`*(y-ybar); * uncorrected;
R2 = ssr/ssto;
sse = (y-yhat)`*(y-yhat);
p = ncol(X);
df_e = n-p;
df_t = n-1;
df_r = p-1;
R2_adj = 1-((sse/(df_e))/(ssto/(df_t)));
print R2, R2_adj;

print df_r, df_e, df_m;

* Question C;
msr = ssr/df_r;
mse = sse/df_e;
F = msr/mse;
F_crit = finv(0.95, df_r, df_e);
print F F_crit;


* Question D;
covb = mse*inv(X`*X);
var_b = vecdiag(covb);
std_b3 = sqrt(var_b[3]);
b3 = (b[3]);
t_value_b3 = b3/std_b3;
t_crit = tinv(0.975, df_r);
print t_value_b3 t_crit;


* Question E;
std_b = sqrt(var_b);
clb_l = b-tinv(0.95, df_e)*std_b;
clb_u = b+tinv(0.95, df_e)*std_b;
print 'lower: ' CLb_L, 'upper: ' CLb_U;

* Question E - Joint sigificance formula;
X_reduced = x[,1]||x[,3]||x[,5];
yhat_reduced = X_reduced*inv(X_reduced`*X_reduced)*X_reduced`*y;
ssr_reduced = (yhat_reduced - ybar)`*(yhat_reduced - ybar);
R2_reduced = ssr_reduced/ssto;
print R2_reduced;
F_reduced = ((R2-R2_reduced)/(1)) / (1-R2_reduced)/(27);
F_reduced_pval = 1-probF(F_reduced, 2, 27);
print F_reduced F_reduced_pval;
quit ;











* using proc reg;
* proc reg data=a plots=none alpha=0.10;
* 	model prop = income stand ratio double / clb;
* run;
	
