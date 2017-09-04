% This file is to compare the data description methods in terms of AMI:
% 1. CHDC - chdc
% 2. Kernel Kmeans
% 3. DBSCAN
% 4. Linkage Clustering

% To prepare the datasets used in this file, please add the folder
% prtool, dd_tools and datasets into the working environment.
% Note that to gain prdatasets, the following commands could be used:
% >> prdatasets             % list the datasets available
% >> prdatasets ALL         % download all the datasets

% For selecting hyperparameters, 10-fold cross-validation is used.
% The hyperparamters for the methods are as follows:
% 1. CHDC - chdc                rejfrac, sigma, tol, reltol
% 2. Kernel Kmeans              c, sigma
% 3. DBSCAN                     epsilon, minptns
% 4. Linkage Clustering         c


%% 0. basic preparation (data / parameter)
load('all-datasets.mat');

M = 4;      % number of models
N = 8;      % number of datasets

D = cell(N,1);      % data set
L = cell(N,1);      % lable set

% D here contains all the data instances with a target class and outliers
% L here contains all the labels with 1 for target class and 2 for outliers

D{1} = bootmotor(:,1:end-1);
L{1} = bootmotor(:,end);

T = unique(carPrestigeLabel);
D{2} = carPrestige;
L{2} = strcmp(carPrestigeLabel,T(1));
L{2} = L{2} + 2*strcmp(carPrestigeLabel,T(2));

T = unique(HistDataOldMapsLabel);
D{3} = HistDataOldMaps;
L{3} = zeros(length(HistDataOldMapsLabel),1);
for i=1:length(T)
    L{3} = L{3} + i*strcmp(HistDataOldMapsLabel,T(i));
end

T = unique(imagesegmentationLabel);
D{4} = imagesegmentation;
L{4} = zeros(length(imagesegmentationLabel),1);
for i=1:length(T)
    L{4} = L{4} + i*strcmp(imagesegmentationLabel,T(i));
end

D{5} = penbasedrecognitionhandwrittendigits(1:1000,1:end-1);
L{5} = penbasedrecognitionhandwrittendigits(1:1000,end);

T = unique(DataUserModelingDatasetHamdiTolgaKAHRAMANSLabel);
D{6} = DataUserModelingDatasetHamdiTolgaKAHRAMANS;
L{6} = zeros(length(DataUserModelingDatasetHamdiTolgaKAHRAMANSLabel),1);
for i=1:length(T)
    L{6} = L{6} + i*strcmp(DataUserModelingDatasetHamdiTolgaKAHRAMANSLabel,T(i));
end

D{7} = drivFaceD.data;
L{7} = drivFaceD.nlab;

D{8} = movementlibras(:,1:end-1);
L{8} = movementlibras(:,end);

w = cell(N,M);              % all the M mapping for N datasets
w_auc = zeros(N,M+1);       % all the M*N aucs
bestarg = cell(N,M);        % all the M*N best parameters


%% 1. loop over each dataset for performance analysis
for i=1:N	% for all the datasets
    i

    A = D{i};
    I = L{i} + 1;

    %% 1.1 preprocess the dataset using a scale mapping
    A = A * scalem(A,'variance');

    %% 1.2 prepare the candidates of the parameters
    % 1. CHDC - chdc                rejfrac, sigma, tol, reltol
    % 2. Kernel Kmeans              c, sigma
    % 3. DBSCAN                     epsilon, minptns
    % 4. Linkage Clustering         c
    
    etrfrac = [0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8 0.9 1];
    sigma1 = sqrt([0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8 0.9 1]);
    sigma2 = sqrt(0.1:0.2:2) * sqrt(size(+A,2));
    sigma = unique([sigma1, sigma2]);
    tol = 1.0e-4;
    reltol = 1e-2;

    C = 1:15;
    epsilon = [0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8 0.9 1];
    minptns = [2 3 5 7 10 15 20 25 50 75];
    
    %% 1.3 calculate best AMI for each method
    % I. chdc
    AMI_chdc = [];
    AMI_chdc_arg = [];
    for j = 1:length(etrfrac)
        for k = 1:length(sigma)
            for l=1:length(reltol)
                clust = chdc_train(A,etrfrac(j),sigma(k),tol,reltol(l));
                AMI_chdc_arg = [AMI_chdc_arg; etrfrac(j) sigma(k) tol reltol(l)];
                AMI_chdc = [AMI_chdc Func_AMI(I,clust+1)];
            end
        end
    end
    [val,idx] = max(AMI_chdc);    
    w_auc(i,1) = val;

    % II. kernel kmeans
    % This method comes from the Matlab and R implementations of the clustering algorithms 
    % in "Localized Data Fusion for Kernel k-Means Clustering with Application to Cancer Biology",
    % which is appearing in Advances in Neural Information Processing Systems 27 (NIPS 2014).
    % https://github.com/mehmetgonen/lmkkmeans
    AMI_KKmeas = [];
    for j = 1:length(C)
        for k = 1:length(sigma)
            clust = Func_KKmeans(A,C(j),sigma(k));
            AMI_KKmeas = [AMI_KKmeas Func_AMI(I,clust+1)];
        end
    end
    w_auc(i,2) = max(AMI_KKmeas);
    
    % III. dbscan
    AMI_DBSCAN = [];
    for j = 1:length(epsilon)
        for k = 1:length(minptns)
            clust = DBSCAN(A,epsilon(j),minptns(k));
            AMI_DBSCAN = [AMI_DBSCAN Func_AMI(I,clust+1)];
        end
    end
    w_auc(i,3) = max(AMI_DBSCAN);
    
    % IV. linkage clustering
    AMI_Linkage = [];
    for k = 1:length(C)
        clust = cluster(linkage(A),'maxclust',C(k));
        AMI_Linkage = [AMI_Linkage Func_AMI(I,clust+1)];
    end
    w_auc(i,4) = max(AMI_Linkage);
    
    % V. benchmarking
    w_auc(i,5) = Func_AMI(I,[1:length(I)]);
end