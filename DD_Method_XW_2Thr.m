%% D Thresholding
% find the threshold for one-class classification

%% D.1 init new X
N = N;
D = sum(I);

X = ones(D,N) ./ D;

K = AtA;
K1 = K(I,:);
K2 = K(I,I);
k = diag(K2);

dnormrec = [];

%% D.2 enter main loop, FIND X that minimize the difference
for i=1:maxiter
    X0 = X;

    %% Multiplicative Updating (NMF version)
    X = max(0,X .* (K1 ./ (K2*X)));
    X(X<1/(D*100)) = 0;
    
    %% Check for convergence
    d = trace(K - 2*K1'*X + X'*K2*X);
    dnorm = d/N;

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

%% D.3 Get the threshold
[thres,indx] = max(diag(K) - 2*diag(K1'*X) + diag(X'*K2*X));

figure(1)
hold on
scatter(A(1,indx),A(2,indx),'gx')
hold off

figure(2)
hold on
plot(dnormrec(1:end))
hold off