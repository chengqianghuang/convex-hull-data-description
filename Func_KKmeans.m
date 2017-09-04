function [clust, state] = Func_KKmeans(X,C,sigma)
% https://github.com/mehmetgonen/lmkkmeans

%set the number of clusters
parameters = struct();
parameters.cluster_count = C;

%initialize the kernel
%should be an N x N matrix containing similarity values between samples
K = exp(-distm(X,X)/(sigma*sigma));
K = (K+K')/2;

%perform training
state = kkmeans_train(K, parameters);

%display the clustering
clust = state.clustering;
