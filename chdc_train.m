% CHDC Convex Hull Data Clustering
% 
%       Cluster = CHDC(A,FRACREJ,SIGMA,TOL,RELATIONTOL)
%
% INPUT
%   A           One-class dataset
%   FRACREJ     Fraction of extreme points (fixed 0.1)
%   SIGMA       Width parameter in the RBF kernel (default = 1)
%   TOL         Tolerance of updating (convergence)
%   RELATIONTOL Tolerance of the relation (use for clustering)
%
% OUTPUT
%   Cluster     Cluster labels for each data instance
%   V           The clustering model
% 

% This file is to clustering data instance according to their weights in
% other data instances, especially the extreme points.

%
% Author: Chengqiang Huang, ch544@exeter.ac.uk
% PhD Candidate in Computer Science, University of Exeter
% 

function [clust,V] = chdc_train(a,fracrej,sigma,tol,reltol)

W = chdd(a,fracrej,sigma,tol);

thres = reltol;         % default = 1.0e-2

X = W.data.X;
N = size(X,2);
D = size(X,1);

label = 1;
clust = zeros(1,N);

% search through data instances to assign label
for i=1:D
    % find related data
    indx = X(i,:) > thres;
    
    indc = max(clust(indx));
    
    if indc == 0                    % not assigned
        clust(indx) = label;
        label = label + 1;
    else                            % assigned (single/multiple clusters)
        indc_set = unique(clust(indx));
        % run over the label and assign a same one
        for j=1:length(indc_set)
            if indc_set(j) == 0
                continue;
            else
                clust(clust(1,:) == indc_set(j)) = indc;
            end
        end
        
        clust(indx) = indc;
    end
end

% immerge cluster (few data) with larger cluster
while 1
    flag = 0;
    
    clusters = unique(clust);
    for i=1:length(clusters)
        indx = (clust == clusters(i));
        if sum(indx) <= N/100
            flag = 1;
            
            Y = X(:,indx);
            Y(Y > thres) = 0;
            
            [~,I] = max(max(Y,[],2));
            
            Y(I,:) = 0;
            X(:,indx) = Y;
            
            [~,J] = max(X(I,:));
            clust(indx) = clust(J);
        end
    end
    
    if flag == 0
        break;
    end
end

V.W = W;
V.cluster = clust;
V.clusterthres = thres;
