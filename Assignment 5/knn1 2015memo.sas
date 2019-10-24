options ls=72 ps=1000 nodate pageno=1 ;

libname ekt "c:\departement\ekt720\knearest" ;


ods rtf file="c:\departement\ekt720\knearest\knn.rtf" style=default;
ods noproctitle ;


title "K-nearest neighbour modelling" ;
proc iml ;
use ekt.yx ;
 read all into yx ;


*input parameters used;
k=3;
step = 1 ;
*end input;
 
nn=nrow(yx) ;

print nn ;

rn = uniform(J(nn,1,0)) ;

yx1 = rn || yx ;
call sort(yx1,{1}) ;
ssize=round(0.8*nn) ;

print ssize ;

print yx1 ;
yx = (J(ssize,1,1) // J(nn-ssize,1,2)) ||
    yx1[,2:4] ;
print yx ;

r=10 ;
x1min = round(min(yx[,3])) - r ;
x1max = round(max(yx[,3])) + r ;
x2min = round(min(yx[,4])) - r ;
x2max = round(max(yx[,4])) + r ;
 
print x1min x1max x2min x2max ;


do ii = x1min to x1max by step ;
 do jj = x2min to x2max by step ;
  yxg = yxg // ( 3 || 555  ||   ii || jj) ;
 end;
end;

yxm = yx // yxg ;

*yxm contains all the data training, test and grid - 
the number 555 has no meaning - just a
placeholder since the actual outcome is not known ;

/*print "All data incl grid" ,  yxm ;*/


dist = J(nrow(yxm),ssize,999999) ;

do i = 1 to nrow(yxm) ;
  do j = 1 to ssize ;
   dist[i,j] =  (yxm[i,3:4] - yx[j,3:4])*
      (yxm[i,3:4] - yx[j,3:4])` ; 
   if i <> j then dist[i,j] = 
       sqrt(dist[i,j]) ;
  end;
end;
dist = dist`; 

*print dist ;

*the section above calculates a distance matrix with 
the rows the 200 training observations and
the columns all the observations that need 
to be classified ;

do m = 1 to nrow(yxm) ;
  dist_u = dist[,m] || yx[1:ssize,2];
  call sort(dist_u,{1}) ;
  *dist_u = dist_u[2:k+1,] ;  *excluding specific obs ;
  dist_u = dist_u[1:k,] ;	 *including specific obs ;
  *print dist_u ;
  teldt = J(1,k,1)*((dist_u[,2]=1) || 
       (dist_u[,2]=2) || (dist_u[,2]=3))  ;
  maxt = max(teldt) ;


  ** section below adds classification categories for 
  the ties ;
   tie=0;
   ggt=0 ;
   do ll = 1 to 3 ;
     if teldt[1,ll] = maxt then do ;
	   tie=tie+1 ;
	   ggt = ggt + ll*10**(ll-1) ;
	 end;
   end;

	*print tie ggt;

    ntie = ntie // ggt ;
   ** end ties section ;
	   
   do bb = 1 to 3 ;
     if teldt[1,bb] = maxt then do ;
       gr = gr // (yxm[m,] || bb) ;
	   bb=200 ;
     end;
   end;
 end;

 * print  gr ;

  pacm = J(1,ssize,1/ssize)*(gr[,2]=gr[,5])[1:ssize,];
  test=nn-ssize ;
  pact = J(1,test,1/test)*(gr[,2]=gr[,5])[ssize+1:nn,];

  print pacm pact;

  gr = gr || ntie ;
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
