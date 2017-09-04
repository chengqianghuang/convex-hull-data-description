% This file contains the code for describing a dataset using CHDD
% with semi-NMF. (Ring Dataset)

% This file is responsible for comparing different algorithms in solving
% the problem of nonnegative least sqaure problem (NNLS)
% - Quadratic Programming (too slow and inaccurate)
% - Multiplicative Updating (Lee, NMF)
% - Multiplicative Updating (Ding, Semi-NMF)
% - Rank-one Residual Iteration (RRI)

% A = [A; 1']
% min ||A - AX||^2
% s.t. X >= 0


%% Comment
% The different choice of W will effect the result:
% 1) a high W will cause the algorithm to run slower, because it needs to
% find a better result satisfying the convex constraint;
% 2) a high W will cause the algorithm to find more accurate results, which
% is a better set of extreme points.
% To conclude with, W exchanges time for better results.


%% A. Prepare A Dataset for Training
% A.1 a training/testing dataset
load('Dat_Ring.mat')
A = Ring2;

N = size(A,2);      % num of data
D = size(A,1);      % dim of data
param = 0.5;        % parameter for gaussian kernel
W = sqrt(D);        % weights for convex constraint

maxiter = 100;
tolx = 1.0e-4;
tolfun = 1.0e-4;

fracrej = 0.1;


%% B. Find the Convex Hull (Support Vectors)
% B.1 argment training dataset
A = [A; W*ones(1,N)];

% B.2 init X
X = ones(N) ./ N;

% B.3 enter main loop, FIND X that minimize the difference
AtA = exp(-distm(A',A')/(param*param));
ata = diag(AtA);

Kpos = max(0,AtA);
Kneg = max(0,-AtA);

sqrteps = sqrt(eps);
nm = numel(AtA);
timer1 = 0;
timer2 = 0;
dnormrec = [];

for i=1:maxiter
    X0 = X;
    
    %% Denote the follow commented parts to see the effect of different algorithms in extreme point identification
    
    % Method 1, Rank-one Residual (RRI): the original way of updating the coefficient matrix X (W)
    %{
    tic;
    for j=1:N
        X(j,:) = max(0,(AtA(j,:) - AtA(j,:) * X + ata(j) * X(j,:))) / ata(j);
    end
    timer1 = timer1 + toc;
    %}
    
    % Method 2, Multiplicative Updating (Semi-NMF version):
    %{
    tic;
    G = X';
    G = G .* sqrt((Kpos + G*Kneg) ./ (Kneg + G*Kpos));
    X = G';
    timer1 = timer1 + toc
    %}
    
    % Method 3, Multiplicative Updating (NMF version):
    % The following code is partially from nnmf
    tic;
    X = max(0,X .* (AtA ./ (AtA*X + eps(AtA))));
    X(X<1/(N*100)) = 0;
    timer1 = timer1 + toc;
    
    %% Check for convergence
    tic;
    d = AtA*X;
    d = trace(AtA) - 2*trace(d) + trace(X'*d); 
    dnorm = d/N;
    timer2 = timer2 + toc;
    
    delta = max(max(abs(X-X0) / (sqrteps+max(max(abs(X0))))));
    
    if i>1
        if delta <= tolx
            break;
        elseif dnorm0-dnorm <= tolfun*max(1,dnorm0)
            break;
        elseif i == maxiter
            break
        end
    end
    
    dnorm0 = dnorm;
    
    dnormrec = [dnormrec; dnorm0];
end

%% C. Change fracrej here to see different effects (extreme points)
% fracrej = 0.3
I = diag(X) >= dd_threshold(diag(X),1-fracrej);

% plot basic contour
figure(1)
clf
hold on
scatter(A(1,:),A(2,:),'b+')
scatter(A(1,I),A(2,I),'ro')
hold off

figure(2)
clf
plot(dnormrec(1:end))

timer1
timer2