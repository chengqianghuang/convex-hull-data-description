% CHDC Convex Hull Data Clustering
% 
%       Cluster = CHDC_TEST(W,A,FRACREJ,SIGMA,TOL,RELATIONTOL)
%
% INPUT
%   V           The mapping for data clustering based on chdd
%   a           The data, whose cluster are to be assigned (row vectors)
%
% OUTPUT
%   Cluster     Cluster labels for each data instance
%   Rank        The dependencies of data instances with the extreme points
%   List        The list of possible clusters
% 

% This file is to clustering data instance according to their weights in
% other data instances, especially the extreme points.

%
% Author: Chengqiang Huang, ch544@exeter.ac.uk
% PhD Candidate in Computer Science, University of Exeter
% 

function [cluster,rank,list] = chdc_test(V,a)

W = V.W.data;

% Note that, a' = S'*X
% and here we go:

%% 1. setup var
D = size(W.sv,1);       % dimension of X (column vector)
N = size(a,1);          % number of column vector in X
t = [+a W.w*ones(N,1)];

atT = exp(-distm(W.sv,t)/(W.s*W.s));
TtT = exp(-distm(t,t)/(W.s*W.s));

%% 2 init X and get concised A
X = ones(D,N) ./ D;
sqrteps = sqrt(eps);

%% 3 enter main loop, FIND X that minimize the difference
for i=1:W.maxiter
    X0 = X;

    X = max(0,X .* (atT ./ (W.a*X)));
    X(X<1/(D*100)) = 0;

    % Check for convergence
    d = trace(TtT - 2*atT'*X + X'*W.a*X);
    dnorm = d/N;
    delta = max(max(abs(X-X0) / (sqrteps+max(max(abs(X0))))));

    if i>1
        if delta <= W.tolx
            break;
        elseif dnorm0-dnorm <= W.tolfun*max(1,dnorm0)
            break;
        elseif i == W.maxiter
            break
        end
    end

    dnorm0 = dnorm;
end

% This code is for one-class classification
% out = diag(TtT) - 2*diag(atT'*X) + diag(X'*W.a*X);

%% The results
rank = X;               % column vectors
list = V.cluster;       % column vector

[~,idx] = max(X);
cluster = V.cluster(idx);   % row vector