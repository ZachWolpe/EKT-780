options ls=72 ps=1000 nodate pageno=1 ;
libname ekt "/folders/myfolders/sasuser.v94/EKT 720/Assignment 5" ;


* _____________________________ KNN1 2015memo _____________________________;
title "K-nearest neighbour modelling";

proc iml;
use ekt.yx;
read all into yx;

x = {1,4,6};
y = {2,4,5};
t = {9,9,9};


k=3;
step=1;
nn=nrow(yx);

rn = uniform(J(nn,1,0));

yx1 = rn || yx;
call sort(yx1, {1});
ssize = round(0.8*nn);
yx = J(ssize,1,1) // J(nn-ssize,1,2) || yx1[,2:4];

r=10;
x1min = round(min(yx[,3])) - r;
x1max = round(max(yx[,3])) + r;
x2min = round(min(yx[,4])) - r;
x2max = round(max(yx[,4])) + r;

print x1min x1max x2min x2max ;


do ii=x1min to x1max by step;
	do jj=x2min to x2max by step;
		yxg = yxg // (3 || 555 || ii || jj);
	end;
end;

yxm = yx // yxg;
*yxm contains all the data training, test and grid - 
the number 555 has no meaning - just a
placeholder since the actual outcome is not known ;

/*print "All data incl grid" ,  yxm ;*/



dist = J(nrow(yxm), ssize, 999999);

do i=1 to nrow(yxm);
	do j=1 to ssize;
		dist[i,j] = (yxm[i,3:4] - yx[j,3:4])*(yxm[i,3:4] - yx[j,3:4])`;
		if i <> j then dist[i,j] =
			sqrt(dist[i,j]);
	end;
end;
dist = dist`;
*the section above calculates a distance matrix with 
the rows the 200 training observations and
the columns all the observations that need 
to be classified ;


do m=1 to 5; * norw(yxm);
	dist_u = dist[,m] || yx[1:ssize,2];
	call sort(dist_u, {1});
	dist_u = dist_u[1:k,];
	teldt = J(1,k,1)*((dist_u[,2]=1) || 
		(dist_u[,2]=2) || (dist_u[,2]=3));
	maxt = max(teldt);
	
	** section below adds classification categories for 
  	the ties ;
  	ties=0;
  	ggt=0;
  	do ll=1 to 3;
  		if teld[1,ll] = maxt then do;
  			tie = tie + 1;
  			ggt = ggt + ll*10**(ll-1);
  		end;
  	end;
  	ntie = ntie // ggt;
  	
  	do bb=1 to 3;
  		if teldt[1,bb] = maxt then do;
  			gr = gr // (yxm[m,] || bb);
  			bb=200; * break;
  		end;
  	end;
end;


pacm = J(1,ssize, 1/ssize)*(gr[,2]=gr[,5])[1:ssize,];
test = nn-ssize;
pact = J(1,test,1/test)*(gr[,2]=gr[,5])[ssize+1:nn,];

print pacm pact;

gr = gr || ntie;
nm1 = {cat gr x1 x2 pgr ptie};

create plotd from gr[colname=nm1];
append from gr;
quit;








 data plotd ;
  set plotd ;
   grplot = pgr ;
   if cat=1 then grplot = gr +10;
   if cat=2 then grplot = gr+100 ;
   *1,2,3 grid -- 10,11,12 training data -- 101,102,103 test data ;

   grplott = ptie ;
   if cat=1 then grplott= gr+1000 ;
   if cat=2 then grplott= gr+2000 ;
   *1,20,21,300,301,320 grid outcomes including ties
   1001,1002,1003 training data
   2001,2002,2003 test data;
 run ;
   

 symbol1 value=dot 
        height=2;

 symbol2 value=star
        height=2;

proc gplot data=plotd ;
 plot x1*x2=grplot ;
run ;


proc gplot data=plotd ;
 plot x1*x2=grplott ;
run ;

proc gplot data=plotd ;
 plot x1*x2=grplott ;
 where cat=1 ;
run ;



proc freq data=plotd ;
 tables ptie * pgr ;
run ;

ods rtf close ;







* _____________________________ KNN1 2015memo2 _____________________________;














options ls=72 ps=1000 nodate /*pageno=1*/ ; 
**title "Data for KNN" ;
data ektyx;
input y x1 x2;
cards;
1	157.005	385.016
1	186.275	355.044
1	88.589	333.571
1	196.616	314.552
1	206.897	318.657
1	223.15	340.759
1	175.402	385.466
1	143.427	313.52
1	88.009	330.403
1	214.029	339.076
1	189.561	343.103
1	197.556	368.517
1	102.856	322.907
1	111.041	341.741
1	205.548	336.41
1	185.191	375.33
1	171.377	393.6
1	88.074	329.661
1	210.341	340.857
1	205.414	336.617
1	215.191	335.33
1	205.859	335.641
1	180.341	380.857
1	189.012	379.435
1	90.097	331.075
1	227.556	328.517
1	94.03	335.814
1	194.406	295.451
1	210.277	344.607
1	210.216	333.389
1	204.218	320.116
1	92.321	333.318
1	177.439	379.478
1	227.85	328.211
1	175.694	373.737
1	87.944	333.032
1	208.819	337.47
1	92.806	327.681
1	178.224	339.481
1	93.884	334.199
1	165.836	356.188
1	132.883	378.589
1	88.069	332.473
1	218.686	334.347
1	176.897	358.657
1	174.218	360.116
1	84.014	320.941
1	159.561	383.103
1	163.567	344.456
1	97.983	335.458
1	177.943	355.171
1	234.056	355.664
1	88.855	331.945
1	133.246	312.301
1	98.279	333.069
1	90.152	333.326
1	90.124	334.414
1	215.537	334.372
1	198.513	359.382
1	176.846	385.648
1	185.537	374.372
1	175.414	376.617
1	91.802	334.036
1	169.037	324.646
1	228.513	319.382
1	216.275	315.044
1	92.476	333.15
1	203.82	347.429
1	97.851	335.45
1	180.216	373.389
1	197.85	368.211
1	206.846	345.648
1	87.414	325.137
1	155.276	324.861
1	153.824	327.19
1	104.421	323.924
1	93.135	331.221
1	80.859	327.02
1	217.011	326.177
1	175.682	383.256
1	193.15	380.759
1	124.268	332.588
1	207.439	339.478
1	205.694	333.737
1	175.548	376.41
1	170.858	313.639
1	184.029	379.076
1	188.686	374.347
1	216.059	366.707
1	95.881	337.593
1	150.858	320.084
1	166.616	354.552
1	133.897	369.977
1	205.402	345.466
1	187.011	366.177
1	219.012	339.435
1	180.277	384.607
1	87.949	330.419
1	88.612	325.604
1	205.682	343.256
2	117.983	299.2
2	216.789	343.307
2	147.439	359.478
2	107.944	306.916
2	144.218	340.116
2	133.845	399.6
2	155.537	354.372
2	112.476	301.955
2	163.15	360.759
2	138.198	408.548
2	165.076	383.961
2	100.859	305.86
2	219.269	310.271
2	110.097	301.516
2	215.687	344.452
2	138.39	307.558
2	146.897	338.657
2	156.275	335.044
2	129.561	363.103
2	132.795	407.867
2	154.029	359.076
2	104.014	293.091
2	112.806	293.312
2	114.03	304.22
2	133.887	363.982
2	108.589	306.998
2	199.202	390.256
2	112.321	302.384
2	117.851	299.337
2	113.884	301.944
2	159.012	359.435
2	183.718	385.304
2	108.855	304.239
2	140.419	385.462
2	107.414	295.58
2	145.682	363.256
2	167.556	348.517
2	151.274	390.622
2	99.92	292.199
2	113.135	298.29
2	107.949	302.959
2	157.011	346.177
2	176.233	335.696
2	115.881	304.811
2	158.686	354.347
2	145.694	353.737
2	145.402	365.466
2	150.341	360.857
2	136.616	334.552
2	108.009	302.867
2	133.343	373.985
2	168.513	339.382
2	108.069	305.928
2	150.809	362.681
2	145.414	356.617
2	110.124	306.532
2	151.035	404.655
2	145.548	356.41
2	150.216	353.389
2	167.85	348.211
2	73.523	399.027
2	150.277	364.607
2	108.074	301.671
2	155.191	355.33
2	111.802	304.059
2	146.846	365.648
2	169.422	369.943
2	108.612	294.927
2	170.717	366.357
2	182.499	366.263
2	118.279	295.251
2	128.365	313.019
2	110.152	304.855
2	173.477	294.041
2	132.842	374.76
3	137.85	368.211
3	192.476	311.955
3	127.011	366.177
3	116.897	358.657
3	120.216	373.389
3	133.15	340.759
3	120.341	340.857
3	129.012	379.435
3	191.802	314.059
3	115.548	336.41
3	192.806	303.312
3	117.439	339.478
3	115.682	343.256
3	187.949	312.959
3	115.694	333.737
3	188.074	311.671
3	133.15	380.759
3	125.537	374.372
3	106.616	354.552
3	187.944	316.916
3	137.85	328.211
3	197.851	309.337
3	120.277	344.607
3	193.135	308.29
3	115.414	336.617
3	115.414	376.617
3	188.612	304.927
3	125.191	335.33
3	188.069	315.928
3	193.884	311.944
3	198.279	305.251
3	116.897	318.657
3	190.097	311.516
3	116.846	385.648
3	195.881	314.811
3	190.124	316.532
3	125.537	334.372
3	117.439	379.478
3	129.012	339.435
3	125.191	375.33
3	116.846	345.648
3	126.275	355.044
3	192.321	312.384
3	188.009	312.867
3	115.694	373.737
3	115.682	383.256
3	187.414	305.58
3	115.548	376.41
3	194.03	314.22
3	180.859	315.86
3	114.218	320.116
3	124.029	379.076
3	188.855	314.239
3	99.561	383.103
3	138.513	359.382
3	138.513	319.382
3	99.561	343.103
3	128.686	334.347
3	115.402	385.466
3	137.556	368.517
3	127.011	326.177
3	126.275	315.044
3	188.589	316.998
3	106.616	314.552
3	190.152	314.855
3	197.983	309.2
3	137.556	328.517
3	115.402	345.466
3	120.277	384.607
3	124.029	339.076
3	128.686	374.347
3	184.014	303.091
3	114.218	360.116
3	120.341	380.857
3	120.216	333.389
;
run;
proc iml ;
use ektyx ;
read all into yx ;


*input parameters used;
k=3;
step = 1 ; /*step size - determines the points of grid*/
*end input;
 
nn=nrow(yx) ;

print nn ;

rn = uniform(J(nn,1,0)) /*random generation*/;

yx1 = rn || yx ;
call sort(yx1,{1}) ;
ssize=round(0.8*nn) ;

print ssize ;

print yx1 ;
yx = (J(ssize,1,1) // J(nn-ssize,1,2)) || yx1[,2:4] ; /*'//' to vertically concatinate row of 1's and 2's*/
print yx ;

r=10 ; /*helps with grid design*/
x1min = round(min(yx[,3])) - r ;  /*r to extend the grid for more values*/
x1max = round(max(yx[,3])) + r ;
x2min = round(min(yx[,4])) - r ;
x2max = round(max(yx[,4])) + r ;
 
print x1min x1max x2min x2max ;


do ii = x1min to x1max by step /*step size*/ ;
 do jj = x2min to x2max by step /*step size*/;
  yxg = yxg // ( 3 || 555  ||   ii || jj) ; /*caterogry number||555(placeholder)||ii(from do loop)||jj(from do loop)*/
 end;
end;
/*double do loop to move vertically and horizontally -> full set of grid point*/

yxm = yx // yxg ;/*1st column=category 2nd=group 3rd=x1 4th=x2*/

*yxm contains all the data training(1), test(2) and grid(3) - 
the number 555 has no meaning - just a
placeholder since the actual outcome is not known ;

/*print "All data incl grid" ,  yxm ;*/


dist = J(nrow(yxm),ssize,.) ;

do i = 1 to nrow(yxm) ;
  do j = 1 to ssize ;
   dist[i,j] =  (yxm[i,3:4] - yx[j,3:4])*(yxm[i,3:4] - yx[j,3:4])` ; /*calculating the distance*/ 
   if i <> j then dist[i,j] = sqrt(dist[i,j]) ; /*'<>' - does not equal to (also use '^=')*/
  end;
end;
dist = dist`; /*transposing distance matrix*/

*print dist ;

*the section above calculates a distance matrix with 
the rows the 200 training observations and
the columns all the observations that need 
to be classified ;

do m = 1 to nrow(yxm) ;
  dist_u = dist[,m] || yx[1:ssize,2]; /*distance corresponds to group*/
  call sort(dist_u,{1}) ; /*sort distances*/
  *dist_u = dist_u[2:k+1,] ;  *excluding specific obs ;
  dist_u = dist_u[1:k,] ;	 *including specific obs ; /*pick k closest points*/
  *print dist_u ;
  teldt = J(1,k,1)*((dist_u[,2]=1)||(dist_u[,2]=2) || (dist_u[,2]=3))  ; /*testing columns for group values. If not, value=0 for that group ||*/
  maxt = max(teldt) ;


  ** section below adds classification categories for 
  the ties (must decide: smaller group, larger group or random assign); 
   tie=0;
   ggt=0;
   do ll = 1 to 3 ;
     if teldt[1,ll] = maxt then do ;
	   tie=tie+1 ;
/*Sollie code*/ggt = ggt + ll*10**(ll-1) ;/*shows tied groups with visual number*/
	 end;
   end;

	*print tie ggt;

    ntie = ntie // ggt ; /*keeps track*/
   ** end ties section ;
	   
   do bb = 1 to 3 ; /*assigns groups randomly in this section*/
     if teldt[1,bb] = maxt then do ;
       gr = gr // (yxm[m,] || bb) ; /*adding  a 5th column with new groups*/
	   bb=200 ; /*stop criteria*/
     end;
   end;
 end;

 * print  gr ;

  pacm = J(1,ssize,1/ssize)*(gr[,2]=gr[,5])[1:ssize,]; /*classification section that compares column 2 with 5*/
  test=nn-ssize ;
  pact = J(1,test,1/test)*(gr[,2]=gr[,5])[ssize+1:nn,];

  print pacm pact;

  gr = gr || ntie ; /*adding  a 6th column with ties data*/
  nm1 = {cat gr x1 x2 pgr ptie} ;

 * print gr[colname=nm1] ;  

  create plotd from gr[colname=nm1] ;
  append from gr ;

  *print ntie ;

quit ;

 data plotd ;
  set plotd ;
   grplot = pgr ;
   if cat=1 then grplot = gr +10;
   if cat=2 then grplot = gr+100 ;
   *1,2,3 grid -- 10,11,12 training data -- 101,102,103 test data ;

   grplott = ptie ;
   if cat=1 then grplott= gr+1000 ;
   if cat=2 then grplott= gr+2000 ;
   *1,20,21,300,301,320 grid outcomes including ties
   1001,1002,1003 training data
   2001,2002,2003 test data;
 run ;
   

 symbol1 value=dot 
        height=2;

 symbol2 value=star
        height=2;

proc gplot data=plotd ;
 plot x1*x2=grplot ;
run ;


proc gplot data=plotd ;
 plot x1*x2=grplott ;
run ;

proc gplot data=plotd ;
 plot x1*x2=grplott ;
 where cat=1 ;
run ;



proc freq data=plotd ;
 tables ptie * pgr ;
run ;

ods rtf close ;













