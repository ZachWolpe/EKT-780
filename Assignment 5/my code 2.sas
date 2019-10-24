options ls=72 ps=1000 nodate pageno=1 ;
libname ekt "/folders/myfolders/sasuser.v94/EKT 720/Assignment 5" ;


* New functions:
	Rand('uniform')
	uniform(J(n,1,1))
	loc()
	sort(x, 1)
;

proc iml;
use ekt.yx; read all into xy;

k=3;
split=0.3;
nn=nrow(xy);

start split_data;
	train = uniform(J(nn,1,1)) > split;
	xy = xy || train;
	xy_train = xy[loc(xy[,4]=1)`,];
	xy_test = xy[loc(xy[,4]=0)`,];
finish split_data;


call split_data;


start fit_model;
	y = xy_data[,1];
	x1 = xy_data[,2];
	x2 = xy_data[,3];
 	yt = xy_test[,1];
	x1t = xy_test[,2];
	x2t = xy_test[,3];
	n = nrow(x1);
	nt = nrow(x1t);
	
	dist = J(n,nt,999);
	do i=1 to n;
		* for each datapoint;
		do j=1 to nt;
			* for each xi,xj pair;
			dist[i,j] = sqrt((x1[i]-x1t[j])**2 + (x2[i]-x2t[j])**2);
		end;
	end;
	
	yy = J(nt,k,999);
	do i=1 to nt;
		* select K nearest neighbors;
		d = y||x1||x2||dist[,i];
		call sort(d,4);
		yy[i,] = d[1:k,1]`;
	end;
	
	count = J(nt,k,0);
	do i=1 to nt;
		* count appearences;
		do j=1 to k;
			if yy[i,j] = 1 then; count[i,1] = count[i,1] + 1;
			if yy[i,j] = 2 then; count[i,2] = count[i,2] + 1;
			if yy[i,j] = 3 then; count[i,3] = count[i,3] + 1;
		end;
	end;
	
	predict = J(nt,k,0);
	do i=1 to nt;
		* predict;
		predict[i,] = whichn(max(count[i,]), count[i,]);
	end;
	
	final_prediction = J(nt,1,0);
	do i=1 to nt;
		do j=1 to k;
			if predict[i,j] = 1 then; final_prediction[i] = j;
		end;
	end;	

finish fit_model;

xy_data = xy_train;
call fit_model;
test_data_pred = final_prediction || x1t || x2t;

accuracy = 0;
do i=1 to nt;
 	if yt[i]=final_prediction[i] then; 
		accuracy = accuracy + 1;
end;
accuracy = accuracy/nrow(yt);
print accuracy;



step=1;
y = xy[,1];
x1 = xy[,2];
x2 = xy[,3];

start gen_points;
	up1 = max(x1) + round(range(x1)/10);
	low1 = min(x1) - round(range(x1)/10);
	up2 = max(x2) + round(range(x2)/10);
	low2 = min(x2) - round(range(x2)/10);
	points = low1 || low2;
	
	xx1 = low1;
	xx2 = low2;
	do while(xx1 < up1);
		xx1 = xx1 + step;
		do while (xx2 < up2);
			xx2 = xx2 + step;
			points = points // (xx1 || xx2);
		end;
		xx2 = low2;
	end;
finish gen_points;


call gen_points;


* reset INPUTS;
* RERUN PREDICTORS;

xy_data = xy_train;
xy_test = J(nrow(points),1,1)||points;
call fit_model;
random_point_pred = final_prediction || x1t || x2t;


emp = J(nrow(random_point_pred),ncol(test_data_pred),.);
do i=1 to nrow(test_data_pred);
 	emp[i,] = test_data_pred[i,];
end;


dat = emp || random_point_pred;
create dat from dat[colname={'y' 'x1' 'x2' 'ypoints' 'x1points' 'x2points'}];
append from dat;
	
create points from random_point_pred[colname={'y' 'x1' 'x2'}];
append from random_point_pred;

create test from test_data_pred[colname={'y' 'x1' 'x2'}];
append from test_data_pred;
quit;


	


proc sgplot data=test;
	scatter x=x1 y=x2 / group=y markerattrs=(symbol=CircleFilled);



proc sgplot data=points;
	scatter x=x1 y=x2 / group=y;

proc sgplot data=dat;
	scatter x=x1 y=x2 / group=y markerattrs=(symbol=CircleFilled)
		markeroutlineattrs=(thickness=5);
	scatter x=x1points y=x2points / group=ypoints
		markeroutlineattrs=(thickness=0.02);

 




