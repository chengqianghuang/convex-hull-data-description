%% E. Testing
% test the effect of CHDD in one-class classification

%% E.1 testing dataset
T = normrnd(0,2,1,500);
T = [T; normrnd(0,2,1,500)];
T = [T; W*ones(1,500)];

N = size(T,2);
D = sum(I);

atT = exp(-distm(A(:,I)',T')/(param*param));
TtT = exp(-distm(T',T')/(param*param));

X = ones(D,N) ./ D;

dnormrec = [];

%% E.2 enter main loop, FIND X that minimize the difference
for i=1:maxiter
    X0 = X;
    
    %% Multiplicative Updating (Lee version)
    X = max(0,X .* (atT ./ (K2*X)));
    X(X<1/(D*100)) = 0;

    %% Check for convergence
    d = trace(TtT - 2*atT'*X + X'*K2*X);
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

%% E.3 Get the anomalies
error = diag(TtT) - 2*diag(atT'*X) + diag(X'*K2*X);

% plot the results
figure(1)
hold on
scatter(T(1,:),T(2,:),'kx')
scatter(T(1,error>thres),T(2,error>thres),'go')
hold off

figure(2)
hold on
plot(dnormrec(1:end))
hold off