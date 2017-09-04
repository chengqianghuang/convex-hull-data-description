% This is README file for the experiments of Convex Hull Data Description

%% 1. The description of the folders and files
% Files in folder "datasets" are data files used for the experiments;
% Files in folder "methods" contains methods to be compared with CHDD;
% Files in folder "tools" contains the framework that CHDD built upon;
% Files in folder "results" maintains the results of the experiments;
% Files in folder "figures" has all the source code for reproducing the
% figures in the paper.

% "DD_Method_*.m" are step-by-step introductions of CHDD;
% "chdd.m" is the file that integrates CHDD into dd_tools framework;
% "chdc_train.m" and "chdc_test.m" are for convex hull clustering;

% "Func_KKmeans.m" is a simple wrapper for lmkkmeans; (Please check ref)
% "Func_AMI.m" is for calculating AMI of clustering results; (Please check ref)

% "Comparison_*.m" are files for the comparison of CHDD with other method
% in one-class classification or clustering tasks.


%% 2. A simple case study to show how to use the programs
% use separate programs to show the whole process of using CHDD
% Run: DD_Method_XW_1Tra      % training
% Run: DD_Method_XW_2Thr      % thresholding
% Run: DD_Method_XW_3Tst      % testing


%% 3. Test the validity of chdd function under the dd_tools framework
% please run all the files in step 2 at first
t = T(1:2,:)';

Q = chdd(Ring2',0.15,0.5,1.0e-4);
res = +(t*Q);

figure(3)
hold on
scatter(Q.data.sv(:,1), Q.data.sv(:,2), 'bo')
scatter(T(1,:),T(2,:),'kx')
scatter(T(1,abs(res(:,1))>Q.data.threshold),T(2,abs(res(:,1))>Q.data.threshold),'go')
hold off


%% 4. References, please refer to the following for more details
% dd_tools: http://prlab.tudelft.nl/david-tax/dd_tools.html
% prtools: http://prtools.org/software/
% dbscan: https://uk.mathworks.com/matlabcentral/fileexchange/52905-dbscan-clustering-algorithm
% kkmeans: https://github.com/mehmetgonen/lmkkmeans