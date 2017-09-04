% This file is to compare the data description methods in terms of AUC:
% 1. CHDD   -	chdd
% 2. SVDD   - 	svdd              dd_tools
% 3. GMM    -	mog_dd            dd_tools
% 4. Parzen -   parzen_dd         dd_tools
% 5. PCA    -	pca_dd            dd_tools
% 6. Kmeans -   kmeans_dd         dd_tools
% 7. KNN    -	knndd             dd_tools
% 8. LOF    -   lofdd             dd_tools
% 9. LPDD   -   lp_dd             dd_tools

% To prepare the datasets used in this file, please add the folder
% prtools, dd_tools and datasets into the working environment.
% Note that to gain prdatasets, the following commands could be used:
% >> prdatasets             % list the datasets available
% >> prdatasets ALL         % download all the datasets

% TO HAVE A RIGHT RESULT, THIS FILE HAS MODIFIED THE LIBRARY FILE:
% dd_gridsearch.m TO TRAIN MODELS ONLY WITH TARGET DATA INSTANCES. PLEASE
% REFER TO THE LIBRARY FILE FOR MORE INFORMATION.

% For selecting hyperparameters, 5-fold cross-validation is used.
% The paramters for the methods are as follows:
% 1. chdd:          extfrac, sigma, tol
% 2. svdd:          rejfrac, sigma
% 3. mog_dd:        rejfrac, n
% 4. parzen_dd:     rejfrac, H
% 5. pca_dd:        rejfrac, npc
% 6. kmeans_dd:     rejfrac, n, tol (for convergence)
% 7. knndd:         rejfrac, neighbor, method
% 8. lofdd:         rejfrac, k
% 9. lp_dd:         rejfrac, sigma, distance type, distance par

%% 0. basic preparation (data / parameter)
M = 9;       % number of models
N = 15;      % number of datasets

D = cell(N,1);      % data set
L = cell(N,1);      % lable set

% D here contains all the data instances with a target class and outliers
% L here contains all the labels with 1 for target class and 2 for outliers

% easy datasets
[D{1},L{1}] = oc_set(iris,1);                   % 50 50 50
[D{2},L{2}] = oc_set(wine,1);                   % 59 71 48
[D{3},L{3}] = oc_set(breast,1);                 % 458 241
[D{4},L{4}] = oc_set(car,1);                    % 384 69 1210 65
[D{5},L{5}] = oc_set(biomed,1);                 % 67 127
[D{6},L{6}] = oc_set(diabetes,1);               % 268 500
[D{7},L{7}] = oc_set(sonar,1);                  % 111 97
[D{8},L{8}] = oc_set(breastdiag,1);             % 357 212

% harder datasets
[D{9},L{9}] = oc_set(glass,1);                  % 70 76 17 51
[D{10},L{10}] = oc_set(liver,1);                % 145 200
[D{11},L{11}] = oc_set(ionosphere,1);           % 225 126
[D{12},L{12}] = oc_set(imox,1);                 % 48 48 48 48
[D{13},L{13}] = oc_set(auto_mpg,1);             % 229 169
[D{14},L{14}] = oc_set(chromo,1);               % sum = 1143
[D{15},L{15}] = oc_set(ecoli,1);                % 143 2 77 2 35 20 5 52

w = cell(N,M);              % all the M mapping for N datasets
w_auc = zeros(N,M);         % all the M*N aucs (for one specific run)
w_auc_avg = zeros(N,M);
w_auc_std = zeros(N,M);

bestarg = cell(N,M);        % all the M*N best parameters


%% Loop over each dataset for performance analysis
for i=1:N   % for all the datasets

    A = D{i};
    I = L{i};

    %% 1. For each dataset, run 10 times to get average performances
    w_stat = zeros(M,10);
    for t=1:10              % test each method in each dataset 10 times
                            % to get the mean and variance of performance.
        i
        t
        
        %% 1.1 preprocess the dataset using a scale mapping
        A = A * scalem(A,'variance');

        target = target_class(A,'target');
        target = target(randperm(length(target)),:);
        outlier = target_class(A,'outlier');
        outlier = outlier(randperm(length(outlier)),:);

        %% 1.2 prepare the training dataset (target class only)
        frac = 0.9;                              % the fraction for training
        lengtar = floor(length(target)*frac);    % the num of target for training
        lengout = floor(length(outlier)*frac);   % the num of target for training    

        a = gendatoc(target(1:lengtar,:),outlier(1:lengout,:));                       % model training data
        b = gendatoc(target(lengtar+1:end,:),outlier(:,:));                           % model testing data

        %% 1.3 prepare the candidates of the parameters
        % 1. chdd:          etrfrac, sigma, tol
        % 2. svdd:          rejfrac, sigma
        % 3. mog_dd:        rejfrac, n
        % 4. parzen_dd:     rejfrac, H
        % 5. pca_dd:        rejfrac, npc
        % 6. kmeans_dd:     rejfrac, n, tol (for convergence)
        % 7. knndd:         rejfrac, k, method
        % 8. lofdd:         rejfrac, k
        % 9. lpdd:          rejfrac, sigma, distance type, distance par

        % In this version of CHDD, the accuracy of the method is determine by
        % simply tuning the number of extreme points.
        % To further improve the usability of the method, it would be better to
        % link the reconstruction error with the number of extreme points. In
        % other words, the selected extreme points should be a minimum set of
        % points that describe all the other data with minimum error. This will
        % be added into later versions. (modify thresholding process)

        etrfrac = [0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5];
        rejfrac = 0;

        sigma = sqrt(0.1:0.2:2) * sqrt(size(+a,2));
        H = sigma.^2;

        tol = 1.0e-5;

        n = 1:10;
        k = 1:10;
        
        npc = 1:floor(sqrt(size(target,2)));

        %% 1.4 5-fold cross-validataion
        % I. chdd
        [w{i,1},bestarg{i,1}] = dd_gridsearch(a,'chdd',5,etrfrac,sigma,tol);

        % II. svdd
        [w{i,2},bestarg{i,2}] = dd_gridsearch(a,'svdd',5,rejfrac,sigma);

        % III. mog_dd
        [w{i,3},bestarg{i,3}] = dd_gridsearch(a,'mog_dd',5,rejfrac,n);

        % IV. parzen_dd
        [w{i,4},bestarg{i,4}] = dd_gridsearch(a,'parzen_dd',5,rejfrac,H);

        % V. pca_dd
        [w{i,5},bestarg{i,5}] = dd_gridsearch(a,'pca_dd',5,rejfrac,npc);

        % VI. kmeans_dd
        [w{i,6},bestarg{i,6}] = dd_gridsearch(a,'kmeans_dd',5,rejfrac,n,tol);

        % VII. knndd
        [w{i,7},bestarg{i,7}] = dd_gridsearch(a,'knndd',5,rejfrac,k);

        % VIII. lofdd
        [w{i,8},bestarg{i,8}] = dd_gridsearch(a,'lofdd',5,rejfrac,k);

        % IVV. lpdd
        [w{i,9},bestarg{i,9}] = dd_gridsearch(a,'lpdd',5,rejfrac,sigma);

        %% 1.5 plot the results
        w_col = ['r' 'b' 'y' 'm' 'k' 'g' 'c' '--' '-.' '-'];

        for j=1:M    
            if isempty(w{i,j}) == 1
                continue;
            end
            e = dd_roc(b*w{i,j});
            w_auc(i,j) = dd_auc(e);
            w_stat(j,t) = w_auc(i,j);
        end
        
        figure          % ROC
        clf
        hold on
        for j=1:M    
            if isempty(w{i,j}) == 1
                continue;
            end    
            e = dd_roc(b*w{i,j});
            w_auc(i,j) = dd_auc(e);
            plotroc(e,w_col(j));
        end
        hold off

    end
    
    %% 2. Get the average performances of all the methods in a dataset
    for j=1:M
        w_auc_avg(i,j) = mean(w_stat(j,:),2);
        w_auc_std(i,j) = sqrt(var(w_stat(j,:),0,2));
    end
end