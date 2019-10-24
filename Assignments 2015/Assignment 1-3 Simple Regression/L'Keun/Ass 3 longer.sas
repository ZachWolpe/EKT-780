data a;
input income1 income2 house stand double prop;
income=income1+income2;
ratio=house/stand;
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

proc iml;
use a;
read all var{prop income stand ratio double} into YX;
tempresults=j(1,16,.);

start regress;

	n=nrow(YX);
	k=ncol(YX);
	Y=YX[,1];
	X=J(n,1,1)||YX[,2:k];

	JJ=j(n,n,1);
	Ybar=(sum(Y)/n)*j(n,1,1);
	dfr=k-1;
	dfe=n-k;
	dft=n-1;
	bhat=inv(X`*X)*X`*Y;

	SSR=bhat`*X`*Y-(1/n)*Y`*JJ*Y;
	SSE=Y`*Y-bhat`*X`*Y;
	SST=SSR+SSE;

	rsq=SSR/SST;
	Adjrsq=1-(SSE/dfe)/(SST/dft);

	Anova=j(3,5,.);
	cn={'DF','SSQ','MSQ','F','p-val'};
	rn={'Model','Error','Corrected Total'};
	Anova[1,1]=k-1;
	Anova[2,1]=n-k;
	Anova[3,1]=n-1;
	Anova[1,2]=SSR;
	Anova[2,2]=SSE;
	Anova[3,2]=SST;
	Anova[1,3]=SSR/dfr;
	Anova[2,3]=SSE/dfe;
	Anova[1,4]=(SSR/dfr)/(SSE/dfe);
	Anova[1,5]=1-probf(((SSR/dfr)/(SSE/dfe)),dfr,dfe);
	*print Anova[colname=cn rowname=rn];

	varcov=(SSE/(n-k))*inv(X`*X);
	stderr=sqrt(vecdiag(varcov));
	tval=bhat/stderr;
	pval=2*(1 - abs(probt(tval,n-k)));

finish;
print '******************a******************';
	run regress;
	print bhat;
	
print '******************b******************';
	holder=yx;
	do i=1 to 150 by 5;
	yx=holder;
	r=remove(YX,i:i+4);
	YYXX=shape(r,29,5);
	yx=yyxx;
	run regress;
	*print bhat tval pval;
	resl=bhat`||tval`||pval`||rsq;
	*print resl;
	tempresults=tempresults//resl;
	end;

colnm={'bh1','bh2','bh3','bh4','bh5','t1','t2','t3','t4','t5',
		 'pt1','pt2','pt3','pt4','pt5','rsq'};
cresults=tempresults[2:31,];
print cresults[colname=colnm];

print '******************c******************';
aresults=cresults[:,];
print aresults[colname=colnm];

print '******************d******************';
seed=1;
call randseed(seed);
sam = sample(1:nrow(YX), 30);  
s30 = YX[sam, ];               
YX= s30;
run regress;
print bhat;

print '******************e******************';
rep=100;
eresults=j(rep,16,777);
do c=1 to rep;
	call randseed(seed);
	sam = sample(1:nrow(YX), 30);  
	s30 = YX[sam, ];               
	YX= s30;
	run regress;
	resl=bhat`||tval`||pval`||rsq;
	eresults[c,]=resl;
	YX=holder;
end;

print eresults[colname=colnm];
aresults=eresults[:,];
print aresults[colname=colnm];
quit;