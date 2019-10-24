
************************ ASSIGNMENT 5 ************************;

/*
Zach Wolpe
u15094813
*/;

data knn_data;
set '/folders/myfolders/sasuser.v94/EKT 720/Assignment 5/yx.sas7bdat';
run;

proc print data=knn_data (obs=10);
run;

* plot the data;
proc sgplot data=knn_data;
	scatter x=x1 y=x2 /group=group;
run;

* train-test-split the data;
proc surveyselect data=knn_data out=knn_new rate=0.8 outall;
run;


* Compute all distances;
title 'Computational Results';
proc iml;
use knn_new;
read all into knn;
total_n = nrow(knn);

* train-test-split;
*********** *********** TRAIN-TEST-SPLIT *********** ***********;
do i=1 to total_n;
	if knn[i,1] = 1 then knn_train = knn_train // knn[i,];
	else knn_test = knn_test // knn[i,];
end;	

n = nrow(knn_train);
group = knn_train[,2];
x1 = knn_train[,3];
x2 = knn_train[,4];


*********** *********** SET K (Groups) *********** ***********;
*********** *********** SET K (Groups) *********** ***********;
k=3;



*********** *********** COMPUTE DISTANCES *********** ***********;
do i=1 to n;
* for each Point i;
	focal = x1[i];
	dist = j(n,1,.);
	
	do j=1 to n;
	* for each Point(i,j) combination;
		d = sqrt((x1[i] - x1[j])**2 +(x2[i]-x2[j])**2);
		dist[j]	= d;
	end;
	* distances = distances || dist;
	
	vec = group || dist;
	call sort(vec, 2);
	predict = predict // (vec[1:k,1])`;
	
end;


*********** *********** SUM GROUPS TO PREDICT *********** ***********;
scores = j(nrow(predict), k, 0);
do i=1 to (nrow(predict)); * for each row;
	
	do j=1 to (ncol(predict)); * for each col in each row;
	
		do kk=1 to k; * for each category;
			if predict[i,j] = kk then scores[i,kk] = scores[i,kk] + 1;
		end;
	end;
end;

**
**
**
;

*********** *********** MAKE PREDICTION *********** ***********;
** 1=predict category, 0=do not predict;
predictions = J(nrow(scores), ncol(scores), 0);
do r=1 to nrow(scores);
	maxx = max(scores[r,]);
	do c=1 to ncol(scores);
		do j=1 to ncol(predictions);
			if scores[r,c]=maxx then predictions[r,c]=1;
		end;
	end;
end;


*********** *********** Handle Tie Breaks *********** ***********;
do i=1 to nrow(scores);
	if sum(predictions[i,]) > 1 then;
		predict[i,k] = 0;
		scores[i,] = J(1,3,0);
		do p=1 to ncol(predict)-1;
			do kk=1 to k;
				if predict[i,p]=kk then scores[i,kk] =scores[i,kk]+1;
			end;
		end;
end;
**
**
**
;
*********** *********** MAKE PREDICTION *********** ***********;
** 1=predict category, 0=do not predict;
predictions = J(nrow(scores), ncol(scores), 0);
do r=1 to nrow(scores);
	maxx = max(scores[r,]);
	do c=1 to ncol(scores);
		do j=1 to ncol(predictions);
			if scores[r,c]=maxx then predictions[r,c]=1;
		end;
	end;
end;






*********** *********** Get Yhat *********** ***********;
yhat = J(nrow(predictions), 1, 0);
do i=1 to nrow(predictions);
	if predictions[i,1]=1 then yhat[i,] = 1;
	if predictions[i,2]=1 then yhat[i,] = 2;
	if predictions[i,3]=1 then yhat[i,] = 3;
end;

* print scores predictions yhat;
results = x1 || x2 || group || yhat;

cm = {'x1' 'x2' 'group' 'yhat'};



create results from results[colname=cm];
append from results;



*********** *********** ********* TESTING DATA ********* *********** ***********;
*********** *********** ********* TESTING DATA ********* *********** ***********;
*********** *********** ********* TESTING DATA ********* *********** ***********;
** COMPUTE DISTANCES W.R.T variables in the model; 

n_test = nrow(knn_test);
group_test = knn_test[,2];
x1_test = knn_test[,3];
x2_test = knn_test[,4];



*********** COMPUTE DISTANCES ***********;
do i=1 to n_test;
	dist = J(n,1,.);
	do j=1 to n; * slice from TESTING data & compare with TRAINING data;
		d = sqrt((x1_test[i] - x1[j])**2 +(x2_test[i]-x2[j])**2);
		dist[j]	= d;
	
	end;
	vec = group || dist;
	call sort(vec, 2);
	predict_test = predict_test // (vec[1:k,1])`;
	
end;

*********** SUM GROUPS TO PREDICT ***********;
scores_test = j(nrow(predict_test), k, 0);
do r=1 to (nrow(predict_test)); * for each row;
	do c=1 to (ncol(predict_test)); * for each col in each row;
		do kk=1 to k; * for each category;
			if predict_test[r,c] = kk then scores_test[r,kk] = scores_test[r,kk] + 1;
		end;
	end;
end;

*********** *********** MAKE PREDICTION *********** ***********;
** 1=predict category, 0=do not predict;
predictions = J(nrow(scores_test), ncol(scores_test), 0);
do r=1 to nrow(scores_test);
	maxx = max(scores_test[r,]);
	do c=1 to ncol(scores_test);
		do j=1 to ncol(predictions);
			if scores_test[r,c]=maxx then predictions[r,c]=1;
		end;
	end;
end;


*********** Handle Tie Breaks ***********;
do i=1 to nrow(scores_test);
	if sum(predictions[i,]) > 1 then;
		predict_test[i,k] = 0;
		scores_test[i,] = J(1,3,0);
		do p=1 to ncol(predict_test);
			do kk=1 to k;
				if predict_test[i,p]=kk then scores_test[i,kk] =scores_test[i,kk]+1;
			end;
		end;
end;
*********** MAKE PREDICTION ***********;
** 1=predict category, 0=do not predict;
predictions = J(nrow(scores_test), ncol(scores_test), 0);
do r=1 to nrow(scores_test);
	maxx = max(scores_test[r,]);
	do c=1 to ncol(scores_test);
		do j=1 to ncol(predictions);
			if scores_test[r,c]=maxx then predictions[r,c]=1;
		end;
	end;
end;

*********** Get Yhat ***********;
yhat_test = J(nrow(predictions), 1, 0);
do i=1 to nrow(predictions);
	if predictions[i,1]=1 then yhat_test[i,] = 1;
	if predictions[i,2]=1 then yhat_test[i,] = 2;
	if predictions[i,3]=1 then yhat_test[i,] = 3;
end;
* print scores predictions yhat;
results_test = x1_test || x2_test || group_test || yhat_test;

cm = {'x1' 'x2' 'group' 'yhat'};


create results_test from results_test[colname=cm];
append from results_test;








quit;


proc print data=results (obs=10);
title 'Modeled Training Data';
run;

proc print data=results_test (obs=10);
title 'Modeled Testing Data';
run;

*********** *********** ********* ASSESS RESULTS ********* *********** ***********;
title 'Confusion Matrix';
proc iml;
use results;
read all var{x1} into x1;
read all var{x2} into x2;
read all var{group} into group;
read all var{yhat} into yhat;

confusion_matrix = J(3,3,0);

do i=1 to nrow(group);
	cell = group[i] || yhat[i];
	confusion_matrix[cell[1], cell[2]] = confusion_matrix[cell[1], cell[2]] + 1;
end;

rn = {'A (y)' 'B (y)' 'C (y)'};
cn = {'A (yhat)' 'B (yhat)' 'C (yhat)'};
mattrib confusion_matrix rowname=rn colname=cn;

do i=1 to ncol(confusion_matrix);
	sumz = sumz || sum(confusion_matrix[,i]);
end;

correct_percentage_per_group = vecdiag(confusion_matrix)` / sumz;

accuracy = sum(vecdiag(confusion_matrix))/sum(confusion_matrix);
print confusion_matrix ,, accuracy ,, correct_percentage_per_group;
quit;

*********** *********** ********* ASSESS RESULTS ********* *********** ***********;
title 'Graphical Results';

proc sgplot data=results;
	scatter x=x1 y=x2 /group=yhat     GROUPDISPLAY=CLUSTER;
run;






*********** *********** ******* TEST DATA RESULTS ******* *********** ***********;
*********** *********** ********* ASSESS RESULTS ********* *********** **********;
title 'Testing Data Confusion Matrix';
proc iml;
use results_test;
read all var{x1} into x1;
read all var{x2} into x2;
read all var{group} into group;
read all var{yhat} into yhat;

confusion_matrix = J(3,3,0);

do i=1 to nrow(group);
	cell = group[i] || yhat[i];
	confusion_matrix[cell[1], cell[2]] = confusion_matrix[cell[1], cell[2]] + 1;
end;

rn = {'A (y)' 'B (y)' 'C (y)'};
cn = {'A (yhat)' 'B (yhat)' 'C (yhat)'};
mattrib confusion_matrix rowname=rn colname=cn;

do i=1 to ncol(confusion_matrix);
	sumz = sumz || sum(confusion_matrix[,i]);
end;

correct_percentage_per_group = vecdiag(confusion_matrix)` / sumz;

accuracy = sum(vecdiag(confusion_matrix))/sum(confusion_matrix);
print confusion_matrix ,, accuracy ,, correct_percentage_per_group;
quit;

*********** *********** ********* ASSESS RESULTS ********* *********** ***********;
title 'Graphical Results';

proc sgplot data=results_test;
	scatter x=x1 y=x2 /group=yhat     GROUPDISPLAY=CLUSTER;
run;



*********** *************** Using K-Fold Cross Validation  *************** ***********;
*********** *************** Using K-Fold Cross Validation  *************** ***********;
*********** *************** Using K-Fold Cross Validation  *************** ***********;
*********** *************** Using K-Fold Cross Validation  *************** ***********;
proc iml;
title 'Using K-Fold Cross Validation';
print 'Using K-Fold Cross Validation';
use knn_data;
read all into d;
K_folds = 10;


folds = round(nrow(d) / K_folds);
do i=1 to K_folds;
	class = class // repeat(i, folds);
end;
d = d || class;


do w=1 to K_folds;
	print 'iteration' w;
	remove x1_test x2_test y_test 
	x1_train x2_train y_train ;
	
	x1_test = J(folds,1,.);
	x2_test = J(folds,1,.);
	y_test = J(folds,1,.);
	x1_train = J((nrow(d)-folds),1,.);
	x2_train = J((nrow(d)-folds),1,.);
	y_train = J((nrow(d)-folds),1,.);

	count_test = 1;
	count_train = 1;
	do j=1 to nrow(d);
		if d[j,ncol(d)] = w then do;
			y_test[count_test,1] = d[j,1];
			x1_test[count_test,1] = d[j,2];
			x2_test[count_test,1] = d[j,3];
			count_test = count_test + 1;
		end;
		if d[j,ncol(d)] ^= w then do;
			y_train[count_train,1] = d[j,1];
			x1_train[count_train,1] = d[j,2];
			x2_train[count_train,1] =  d[j,3];
			count_train = count_train + 1;
		end;
	end;
		
	** ARANGE DATA;
	x1 = x1_train;
	x2 = x2_train;
	group = y_train;
	
	
	K = 3;
	n = nrow(x1_train);
	remove predict_test;
	
	*********** *********** ********* TESTING DATA ********* *********** ***********;
	** COMPUTE DISTANCES W.R.T variables in the model; 
	
	n_test = nrow(y_test);
	group_test = y_test;
	

	
	*********** COMPUTE DISTANCES ***********;
	remove dist;
	remove vec;
	remove predict_test;
	do i=1 to n_test;
		dist = J(n,1,.);
		do j=1 to n; * slice from TESTING data & compare with TRAINING data;
			d_T = sqrt((x1_test[i] - x1[j])**2 +(x2_test[i]-x2[j])**2);
			dist[j]	= d_T;
		end;
		vec = group || dist;
		call sort(vec, 2);
		predict_test = predict_test // (vec[1:k,1])`;
		
	end;


	*********** SUM GROUPS TO PREDICT ***********;
	scores_test = j(nrow(predict_test), k, 0);
	do r=1 to (nrow(predict_test)); * for each row;
		do c=1 to (ncol(predict_test)); * for each col in each row;
			do kk=1 to k; * for each category;
				if predict_test[r,c] = kk then scores_test[r,kk] = scores_test[r,kk] + 1;
			end;
		end;
	end;
	
	*********** *********** MAKE PREDICTION *********** ***********;
	** 1=predict category, 0=do not predict;
	remove predictions;
	predictions = J(nrow(scores_test), ncol(scores_test), 0);
	do r=1 to nrow(scores_test);
		maxx = max(scores_test[r,]);
		do c=1 to ncol(scores_test);
			do j=1 to ncol(predictions);
				if scores_test[r,c]=maxx then predictions[r,c]=1;
			end;
		end;
	end;


	*********** Handle Tie Breaks ***********;
	do i=1 to nrow(scores_test);
		if sum(predictions[i,]) > 1 then;
			predict_test[i,k] = 0;
			scores_test[i,] = J(1,3,0);
			do p=1 to ncol(predict_test);
				do kk=1 to k;
					if predict_test[i,p]=kk then scores_test[i,kk] =scores_test[i,kk]+1;
				end;
			end;
	end;
	

	*********** MAKE PREDICTION ***********;
	** 1=predict category, 0=do not predict;
	remove predictions;
	predictions = J(nrow(scores_test), ncol(scores_test), 0);
	do r=1 to nrow(scores_test);
	remove maxx;
		maxx = max(scores_test[r,]);
		do c=1 to ncol(scores_test);
			do j=1 to ncol(predictions);
				if scores_test[r,c]=maxx then predictions[r,c]=1;
			end;
		end;
	end;

	
	*********** Get Yhat ***********;
	
	remove yhat_Test;
	yhat_test = J(nrow(predictions), 1, 0);
	do i=1 to nrow(predictions);
		if predictions[i,1]=1 then yhat_test[i,] = 1;
		if predictions[i,2]=1 then yhat_test[i,] = 2;
		if predictions[i,3]=1 then yhat_test[i,] = 3;
	end;

	* print scores predictions yhat;
	
	results_test = x1_test || x2_test || group_test || yhat_test;
	
	end;	

	

	
	
	
	x1 = results_test[,1];
	x2 = results_test[,2];
	group = results_test[,3];
	yhat = results_test[,4];
	
	
	title 'Testing Data Confusion Matrix';
	
	confusion_matrix = J(3,3,0);
	
	do i=1 to nrow(group);
		cell = group[i] || yhat[i];
		confusion_matrix[cell[1], cell[2]] = confusion_matrix[cell[1], cell[2]] + 1;
	end;
	
	remove sumz;
	do i=1 to ncol(confusion_matrix);
		sumz = sumz || sum(confusion_matrix[,i]);
	end;
	
	correct_percentage_per_group = vecdiag(confusion_matrix)` / sumz;
	
	accuracy = sum(vecdiag(confusion_matrix))/sum(confusion_matrix);
	
	
	
	
	
	tt = accuracy || correct_percentage_per_group;
	
	final_results = final_results // tt;
	
	print final_results;


end;


tm = {'accuracy' 'Percent Correct A' 'Percent Correct B' 'Percent Correct C'};

create final_data from final_results[colname=tm];
append from final_results;
run;



proc print data = final_data;
run;







	































