 
* ___________________________________ Weighted KNN ___________________________________;
libname lib '/folders/myfolders/sasuser.v94/EKT 720/weighted_KNN/';

data xy;
set lib.yx;
run; 
proc print data=xy (obs=4);




title 'K Nearest Neighbours';
proc iml;
use xy; read all into xy;
n = nrow(xy);


* ------------------------- Train Test Split -------------------------;
start train_test_split;
	in = uniform(J(n,1,1))<train_percentage;
	xy_train = xy[loc(in=1),];
	xy_test = xy[loc(in=0),];
finish train_test_split;
* ------------------------- Train Test Split -------------------------;
train_percentage = 0.8;
call train_test_split;





* --------------------- Create Visualization Grid ---------------------;
x1_max = max(xy_train[,2]) + 10;
x1_min = min(xy_train[,2]) - 10;
x2_max = max(xy_train[,3]) + 10;
x2_min = min(xy_train[,3]) - 10;

do i=x1_min to x1_max by 5;
	do j=x2_min to x2_max by 5;
		grid = grid // (888||i||j);
	end;
end;
* --------------------- Create Visualization Grid ---------------------;






* ------------------------- K Nearest Neighbors Model -------------------------;
start K_nearest_neighbors(xy_train, xy_test, K);
	
	* compute distance matrix;
	dist = distance(xy_train[,2:3], xy_test[,2:3]);

	
	* sort data & collect K nearest;
	weightz = J(nrow(xy_test), K, 0);
	do i=1 to ncol(dist);
		* ---------- Weighting Function --------- ;
		* weight = 1/d;
		col = dist[,i] || (1/dist[,i]) || xy_train[,1];
		call sort(col, {1});
		* weight predictions;
		do j=1 to K;
			if col[j,3]=1 then; weightz[i,1] = weightz[i,1] + col[j,2];
			if col[j,3]=2 then; weightz[i,2] = weightz[i,2] + col[j,2];
			if col[j,3]=3 then; weightz[i,3] = weightz[i,3] + col[j,2];
		end;
	end;
	
	
	* predict;
	do i=1 to nrow(weightz);
		pred = pred // max(weightz[i,<:>]);
	end;
	
	return pred;
finish K_nearest_neighbors;
* ------------------------- K Nearest Neighbors Model -------------------------;
K=3;

* run on testing data;
pred = K_nearest_neighbors(xy_train, xy_test, K);


* ------------------------------ Confusion Matrix ------------------------------;
confusion_matrix = J(K,K,0);
do i=1 to nrow(pred);
	r=xy_test[i,1]; c=pred[i];
	confusion_matrix[r,c] = confusion_matrix[r,c] + 1;
end;

percent_confusion_matrix = confusion_matrix / (sum(confusion_matrix));
accuracy = sum(diag(confusion_matrix))/sum(confusion_matrix);
print confusion_matrix percent_confusion_matrix accuracy;
* ------------------------------ Confusion Matrix ------------------------------;


* run on grid;
grid_pred = K_nearest_neighbors(xy_train, grid, K);


grid_results = (grid_pred || grid[,2:3] || J(nrow(grid_pred),1,1));
res = (xy_train || J(nrow(xy_train),1,2)) // grid_results;

create res from res[colname={'yh' 'x1' 'x2' 'set'}];
append from res;
quit;





* ----------------------------------- Visualize Results -----------------------------------;
data rez;
set res;
if set=1 
	then do;
		grid_x1 = x1;
		grid_x2 = x2;
		grid_yh = yh;
	end;
if set=2 
	then do;
		true_x1 = x1;
		true_x2 = x2;
		true_yh = yh;
	end;
run;

proc sgplot data=rez;
	scatter x=true_x1 y=true_x2 / group=true_yh markerAttrs=(symbol=CircleFilled);
	scatter x=grid_x1 y=grid_x2 / group=grid_yh markerAttrs=(symbol=X);
title 'Grid Predictions';
run;
* ----------------------------------- Visualize Results -----------------------------------;

		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		


data a (keep = age chd )  ;
input nr age chd @@;
cards ;
1 20 0 2  23 0 3  24 0 4  25 0
5 25 1 6  26 0 7  26 0 8  28 0
9 28 0 10 29 0 11 30 0 12 30 0
13 30 0 14 30 0 15 30 0 16 30 1
17 32 0 18 32 0 19 33 0 20 33 0
21 34 0 22 34 0 23 34 1 24 34 0
25 34 0 26 35 0 27 35 0 28 36 0
29 36 1 30 36 0 31 37 0 32 37 1
33 37 0 34 38 0 35 38 0 36 39 0
37 39 1 38 40 0 39 40 1 40 41 0
41 41 0 42 42 0 43 42 0 44 42 0
45 42 1 46 43 0 47 43 0 48 43 1
49 44 0 50 44 0 51 44 1 52 44 1
53 45 0 54 45 1 55 46 0 56 46 1
57 47 0 58 47 0 59 47 1 60 48 0
61 48 1 62 48 1 63 49 0 64 49 0
65 49 1 66 50 0 67 50 1 68 51 0
69 52 0 70 52 1 71 53 1 72 53 1
73 54 1 74 55 0 75 55 1 76 55 1
77 56 1 78 56 1 79 56 1 80 57 0
81 57 0 82 57 1 83 57 1 84 57 1
85 57 1 86 58 0 87 58 1 88 58 1
89 59 1 90 59 1 91 60 0 92 60 1
93 61 1 94 62 1 95 62 1 96 63 1
97 64 0 98 64 1 99 65 1 100 69 1
;
run ;

proc means data=a;
		
		
		
		
		
		
		
		
		
		
		
proc iml;
use a;
read all into xy;
n=nrow(xy);
x = J(n,1,1)||xy[,1];
y = xy[,2];
x_orig=x;
y_orig=y;

do i =1 to 100;
in = sample((1:n), n, 'replace')`;
end;



*--Newton--*;
start logistic_regression;


	* initialize Beta=0 ==> p=1;
	bo = {0,0};
	diff=999;


	do i=1 to 20 while (diff>0.001);
		lo = x*bo;
		o = exp(lo);
		p = o/(1+o);
		w = diag(p#(1-p)); *think like binomial formula, this is step 3--the diagonal weighted matrix W;


		*--calculating the gradient matrix and hessian matrix--**;
		bn = bo + inv(x`*w*x)*x`*(y-p);
		diff = max(abs(bn-bo));
		bo = bn;
	end;
finish logistic_regression;


call logistic_regression;

* +++ BootStrap Pairs +++;
do i=1 to 100;
	in = sample((1:n), n, 'replace')`;
	y = y_orig[in];
	x = x_orig[in,];
		
	
	call logistic_regression;
	results = results // bn;
end;

print results;
quit ;



