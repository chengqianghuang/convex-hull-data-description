% This file contains the code for convex hull approximation, i.e.,
% identifying the most possible extreme points of a given dataset to
% describe the dataset using the concept of convex hull.

% This file is responsible for comparing different algorithms in solving
% the formulation of the proposed convex hull approximation.

% 1) The formulation:
% A = [A; 1']
% min ||A - AX||^2
% s.t. X >= 0

% 2) The tested algorithms:
% - Multiplicative Updating (MU, Lee & Seung)
% - Rank-one Residual Iteration (RRI, Ho)


%% A. Prepare A Dataset for Training
% A.1 a training/testing dataset
load('Dat_Ring.mat')    % data instances are column vectors
load('Dat_moons.mat');
load('Dat_blobs.mat');
load('Dat_Multi.mat');

for data_index = 1:4
    
    if data_index ==1
        A = [Ring2 Ring3];
    elseif data_index==2
        A = moons';
    elseif data_index==3
        A = blobs';
    else
        A = Multi;
    end

    A = (A' * scalem(A','variance'))';

    N = size(A,2);      % num of data
    D = size(A,1);      % dim of data

    param = 0.3;        % parameter for gaussian kernel
    W = sqrt(D);        % weights for convex constraint

    maxiter = 100;
    tolx = 1.0e-4;
    tolfun = 1.0e-4;

    fracrej = 0.1;

    %% B. Find the Convex Hull (Support Vectors)
    % B.1 argment training dataset
    A = [A; W*ones(1,N)];

    % B.2 enter main loop, FIND X that minimize the difference
    if param ~= 0
        % 1) Gaussian Kernel
        AtA = exp(-distm(A',A')/(param*param));
    else
        % 2) Linear Kernel
        AtA = A'*A;
    end
    ata = diag(AtA);

    Kpos = max(0,AtA);
    Kneg = max(0,-AtA);

    sqrteps = sqrt(eps);
    nm = numel(AtA);

    %% ======= Method 1, Modified Multiplicative Updating (MU): =========
    % The multiplicative updating methdo from chris ding is applicable to
    % semi-NMF problems, which is:
    %                   X_+- = F_+- * G_+-^T.
    % But the sqrt function in the updating process is neglected.
    % ===================================================================
    X = ones(N) ./ N;
    timer1 = 0;
    timer2 = 0;
    dnormrecMU2 = [];

    for i=1:maxiter
        X0 = X;

        % The following code is partially from nnmf
        tic;
        % X = max(0, X .* (AtA ./ (AtA*X + eps(AtA))));
        X = max(0,X .* ((Kpos + Kneg*X) ./ (Kneg + Kpos*X + eps(AtA))));              % Chris Ding, Modified Semi-NMF (quick and broadly applicable)
        X(X<1/(N*100)) = 0;
        timer1 = timer1 + toc;

        % Check for convergence
        tic;    
        % d = trace(AtA) - 2 * trace(AtA*X) + trace(X'*AtA*X);
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
                break;
            end
        end

        dnorm0 = dnorm;

        dnormrecMU2 = [dnormrecMU2; dnorm0];
    end

    I1 = diag(X) >= dd_threshold(diag(X),1-fracrej);            % c_ii is among top 20%
    I2 = diag(X) >= dd_threshold(diag(X),1-2*fracrej);          % c_ii is among top 40%
    I3 = diag(X) >= dd_threshold(diag(X),1-4*fracrej);          % c_ii is among top 60%
    I4 = diag(X) >= dd_threshold(diag(X),1-8*fracrej);          % c_ii is among top 80%

    % plot basic contour
    figure(1)
    subplot(2,4,data_index)
    hold on
    scatter(A(1,:),A(2,:),'yo')
    scatter(A(1,I4),A(2,I4),'bo')
    scatter(A(1,I3),A(2,I3),'go')
    scatter(A(1,I2),A(2,I2),'go')
    scatter(A(1,I1),A(2,I1),'ro')
    hold off

    timer1
    timer2

    %% ========== Method 2, Rank-one Residue Iteration (RRI): ===========
    % RRI is applicable to:
    %                   X_+- = F_+- * G_+-^T.
    % ===================================================================
    X = ones(N) ./ N;
    timer1 = 0;
    timer2 = 0;
    dnormrecRRI = [];

    for i=1:maxiter
        X0 = X;

        % The following code is partially from nnmf
        tic;    
        for j=1:N
            % kernelized method
            X(j,:) = max(0,(AtA(j,:) - AtA(j,:) * X + ata(j) * X(j,:))) / ata(j);
        end
        X(X<1/(N*100)) = 0;
        timer1 = timer1 + toc;

        % Check for convergence
        tic;    
        % d = trace(AtA) - 2 * trace(AtA*X) + trace(X'*AtA*X);
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
                break;
            end
        end

        dnorm0 = dnorm;

        dnormrecRRI = [dnormrecRRI; dnorm0];
    end

    I1 = diag(X) >= dd_threshold(diag(X),1-fracrej);            % c_ii is among top 20%
    I2 = diag(X) >= dd_threshold(diag(X),1-2*fracrej);          % c_ii is among top 40%
    I3 = diag(X) >= dd_threshold(diag(X),1-4*fracrej);          % c_ii is among top 60%
    I4 = diag(X) >= dd_threshold(diag(X),1-8*fracrej);          % c_ii is among top 80%

    % plot basic contour
    figure(1)
    subplot(2,4,4+data_index)
    hold on
    scatter(A(1,:),A(2,:),'yo')
    scatter(A(1,I4),A(2,I4),'bo')
    scatter(A(1,I3),A(2,I3),'go')
    scatter(A(1,I2),A(2,I2),'go')
    scatter(A(1,I1),A(2,I1),'ro')
    hold off

    timer1
    timer2

end