
FILENAME REFFILE '/folders/myfolders/sasuser.v94/EKT 720/midterm/q2.csv';

PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=WORK.q2
	REPLACE;
	GETNAMES=YES; 
RUN;


proc sgplot data=q2;
	scatter x=x1 y=x2 / group=y;
run;



proc iml;
use q2;
read all into xy;
n=nrow(xy);



* ------------------------- K Nearest Neighbors Model: Y[1,2] -------------------------;
start K_nearest_neighbors(xy_train, xy_test, K);
	
	* compute distance matrix;
	dist = J(nrow(xy_train), nrow(xy_test), 99);
	do i=1 to nrow(xy_train);
		do j=1 to nrow(xy_test);
			dist[i,j] = sqrt((xy_train[i,1]-xy_test[j,1])**2 + (xy_train[i,2]-xy_test[j,2])**2);
		end;
	end;
	
	* sort data & collect K nearest;
	weightz = J(nrow(xy_test), 2, 0);
	do i=1 to ncol(dist);
		* ---------- Weighting Function --------- ;
		* weight = 1/d;
		col = dist[,i] || (1/dist[,i]) || xy_train[,3];
		call sort(col, {1});
		* weight predictions;
		do j=1 to K;
			if col[j,1]=1 then; weightz[i,1] = weightz[i,1] + col[j,2];
			if col[j,1]=2 then; weightz[i,2] = weightz[i,2] + col[j,2];
		end;
	end;
	
	* predict;
	do i=1 to nrow(weightz);
		pred = pred // max(weightz[i,<:>]);
	end;
	
	return pred;
finish K_nearest_neighbors;
* ------------------------- K Nearest Neighbors Model: Y[1,2] -------------------------;



* __________ metric: accuracy __________;
start acc(y,yh);
	acc = 0;
	do i=1 to nrow(y);
		if y[i]=yh[i] then; acc = acc + 1;
	end;
	acc = acc/nrow(yh);
	return acc;
finish acc;
* __________ metric: accuracy __________;




* --- --- --- --- --- --- K Fold Cross Validation --- --- --- --- --- --- ;
g=n/10; kk=100;
accuracy_scores = J(kk,10,999);
count = 0;

do i=1 to n by g;
	count = count + 1;
	* _____________ Moving Train Test Split _____________;
	in = (i:i+g-1)`;
	inv_in = (1:i-1)`//(i+g:n)`; 
	if i=1 then; inv_in = (i+g:n)`;
	if i=(n-(g-1)) then inv_in = (1:n-g)`;
	test = xy[in,];
	train = xy[inv_in,];
	* _____________ Moving Train Test Split _____________;
	
	
	do k=1 to kk;
	* _____________ for various values of K _____________;
		pred = K_nearest_neighbors(train,test,k);
		accuracy = acc(test[,3],pred);
		accuracy_scores[k,count] = accuracy;
	* _____________ For various values of K _____________;
	end;

end;
print accuracy_scores;



	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	