data a;
input X	Y;
cards;
10	44
11	34
12	19
13	14
14	24
15	29
16	26
17	22
18	16
19	18
20	19
21	23
22	31
23	80
24	38
25	39
26	41
27	42
28	36
29	28
30	24
31	18
32	18
33	25
34	29
35	30
36	30
37	34
38	32
39	36
;
run;
proc gplot data=a;
plot y*x;
run;
proc iml;
use a;
read all into xy;
n=nrow(xy);
span=0.7; 
n1=0.8*n;
m=round(span*n1); /*number of poins to include in the window*/

randsam=j(n,1,1);
randgen=rannor(randsam)||xy;
call sort (randgen, {1});
xy1=randgen[1:n1,2:3];
rows_xy1=nrow(xy1);
*nxy=nrow(xy1); 
*print xy1; /*checkpoint*/

minxy1=min(xy1[,1]);
maxxy1=max(xy1[,1]); /*can also use sort instead*/
**print minxy1 maxxy1;

/*kernel*/
do x_focal=1 to rows_xy1;
xc=xy1[,1]-xy1[x_focal,1];
xd=abs(xc);
mhood=xy1||xc||xd; /*neighborhood*/
call sort(mhood,{4}); /*sort according to abs differences*/
mhood1=mhood[1:m,];
h=0.5*(max(mhood1[,1])-min(mhood1[,1])); /*half the length of the window*/

z=mhood1[,4]/h; /*equation for abs(z) in the tricube weight*/
wt=((1-(z)##3)##3)#(z<1); /*condition matrix filled with 1's and 0's*//*'##'=element power*/
**print wt; /*checkpoint*/
y_hat=(wt`*mhood1[,2])/(j(1,m,1)*wt); /*calculation of weighted averages of the y values*/
**ker =ker//(x_focal||y0_hat);

**print y0_hat;

/**LPR with df=m-2**/
w=diag(wt);
X=J(nrow(mhood1),1,1)||mhood1[,3];
Y=mhood1[,2];
bhat_lpr=inv(x`*w*x)*x`*w*y;
bhat=bhat//bhat_lpr[1];

yhat=x*bhat_lpr;
**print yhat;
sse=(y-yhat)`*w*(y-yhat);
mse=sse/(m-2);
varb=mse*inv(x`*w*x);

k=tinv(0.975, m-2);
lowlim=lowlim//(bhat_lpr[1]-(k*(sqrt(varb[1,1]))));
uplim=uplim//(bhat_lpr[1]+(k*(sqrt(varb[1,1]))));
end;
lprmtrx=xy1||bhat||lowlim||uplim;

print lowlim uplim;
create assnglpr4 from lprmtrx[colname={'x' 'y' 'yhat' 'lowlim' 'uplim'}];
append from lprmtrx;


/*create assngkern4 from ker[colname={'x_focal' 'y_hat'}];*/
quit;


goptions reset = all;
symbol1 v=dot c=black h=1;
symbol2 c=red h=1;

/*proc gplot data=assngkern4;
plot y_hat*x_focal/overlay legend;
run;*/

proc gplot data=assnglpr4;
plot (y yhat lowlim uplim)*x/overlay legend /*vaxis=2 haxis=2*/;
run;

