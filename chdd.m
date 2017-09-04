% CHDD Convex Hull Data Description
% 
%       W = CHDD(A,FRACREJ,SIGMA,TOL)
%       W = A*CHDD([],FRACREJ,SIGMA,TOL)
%       W = A*CHDD(FRACREJ,SIGMA,TOL)
%
% INPUT
%   A           One-class dataset
%   FRACREJ     Fraction of extreme points (default = 0.1)
%   SIGMA       Width parameter in the RBF kernel (default = 1)
%   TOL         Tolerance of updating (for early stop, default = 1.0e-5)
%
% OUTPUT
%   W         Convex hull data description
% 
% DESCRIPTION
% CHDD finds a set of extreme points for approximating the convex hull of a
% given dataset. With the convex hull, CHDD is able to describe the dataset
% and identify possible anomalies. The data description uses the 
% Gaussian kernel by default. FRACREJ gives the fraction of the extreme
% points out of all the points, i.e., control the number of extreme points
% and how well the convex hull is approximated.

% Note: This work is built based on the framework of dd_tools.
% http://prlab.tudelft.nl/david-tax/dd_tools.html

% Author: Chengqiang Huang, ch544@exeter.ac.uk
% PhD Candidate in Computer Science, University of Exeter


%function W = chdd(a,fracrej,sigma)
function W = chdd(varargin)

argin = shiftargin(varargin,'scalar');
argin = setdefaults(argin,[],0.1,1,1.0e-5);

if mapping_task(argin,'definition')
   W = define_mapping(argin,'untrained','CHDD');

elseif mapping_task(argin,'training')
   [a,fracrej,sigma,tol] = deal(argin{:});
   
	if isempty(sigma)
		error('This versions needs a sigma.');
	end
	% introduce outlier label for outlier class if it is available.
	if isocset(a)
        a = target_class(a);
		signlab = getoclab(a);
		if all(signlab<0), error('CHDD needs target objects!'); end
    else
        signlab = ones(size(a,1),1);
    end

	% Standard optimization procedure: A ~~ A*X
    %% A. Set up the parameters
    N = size(+a,1);      % num of data
    D = size(+a,2);      % dim of data
    w = sqrt(D);         % weights for convex constraint

    maxiter = 100;
    tolx = tol;
    tolfun = tol;

    %% B. Find the Convex Hull (Support Vectors)
    % B.1 argment training dataset
    t = [+a w*ones(N,1)];

    % B.2 init X
    X = ones(N) ./ N;

    % B.3 enter main loop, FIND X that minimize the difference
    if sigma == 0
        K = t*t';
    else
        K = exp(-distm(t,t)/(sigma*sigma));
    end
    Kpos = max(0,K);
    Kneg = max(0,-K);

    sqrteps = sqrt(eps);
    dnormrec = nan(maxiter,1);
    
    for i=1:maxiter
        X0 = X;
        
        %% Update the solution
        X = max(0,X .* ((Kpos+Kneg*X) ./ (Kneg + Kpos*X + eps(K))));
        X(X<1/(N*100)) = 0;

        %% Check for convergence
        d = K*X;
        d = trace(K) - 2*trace(d) + trace(X'*d);
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

        dnormrec(i) =  dnorm0;
    end
    
    % B.4 get the extreme points:
    I = diag(X) >= dd_threshold(diag(X),1-fracrej);
    alf = t(I,:);
    
    
    %% C Thresholding
    % C.1 init new X
    N;
    D = sum(I);
    
    X = ones(D,N) ./ D;
    
    K1 = K(I,:);
    K2 = K(I,I);
    K1pos = max(0,K1);
    K1neg = max(0,-K1);
    K2pos = max(0,K2);
    K2neg = max(0,-K2);

    % C.2 enter main loop, FIND X that minimize the difference
    for i=1:maxiter
        X0 = X;
        
        %% Update the solution
        X = max(0,X .* ((K1pos + K2neg*X) ./ (K1neg + K2pos*X)));
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
    end
    
    thres = max(diag(K) - 2*diag(K1'*X) + diag(X'*K2*X));
    % thres = mean(diag(K) - 2*diag(K1'*X) + diag(X'*K2*X));
    % thres = median(diag(K) - 2*diag(K1'*X) + diag(X'*K2*X));
    
	%% Compute the offset (not important, but now gives the possibility to
	% interpret the output as the distance to the center of the sphere)
	offs = 1;

	% store the results
    W.maxiter = maxiter;
    W.tolx = tolx;
    W.tolfun = tolfun;
	W.s = sigma;                % sigma for gaussian kernel
    W.w = w;                    % weight for convex constraint
    W.si = I;                   % index of support vectors in original data
	W.sv = alf;                 % argmented support vectors [sv 1]
	W.a = K2;                   % kernel matrix of support vectors
    W.X = X;                    % parameter matrix A = AX (clustering)
	W.threshold = thres;
	W.offs = offs;
    W.p = [];
    W.exitflag = 1;
    W.output = dnormrec;    
	W = prmapping(mfilename,'trained',W,char('target','outlier'),size(a,2),2);
	W = setname(W,'CHDD');
    
elseif mapping_task(argin,'trained execution') %testing

   [a,fracrej] = deal(argin{1:2});
	W = getdata(fracrej);
	m = size(a,1);

    % check if alpha's are OK
    if isempty(W.a)
        warning('dd_tools:OptimFailed','The CHDD is empty or not well defined');
        out = zeros(m,1);
    else
        % and here we go:
        %% 1. setup var
        D = size(W.sv,1);
        N = size(a,1);
        t = [+a W.w*ones(N,1)];
        
        if W.s == 0
            atT = W.sv * t';
            TtT = t * t';
        else
            atT = exp(-distm(W.sv,t)/(W.s*W.s));
            TtT = exp(-distm(t,t)/(W.s*W.s));
        end        
        atTpos = max(0,atT);
        atTneg = max(0,-atT);
        apos = max(0,W.a);
        aneg = max(0,-W.a);

        %% 2. init X and get concised A
        X = ones(D,N) ./ D;
        sqrteps = sqrt(eps);

        %% 3. enter main loop, FIND X that minimize the difference
        for i=1:W.maxiter
            X0 = X;
            
            % Update X
            X = max(0,X .* ((atTpos + aneg*X) ./ (atTneg + apos*X)));
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

        %% 4. Output
        out = diag(TtT) - 2*diag(atT'*X) + diag(X'*W.a*X);
    end
    
    newout = [out repmat(W.threshold,m,1)];
        
	% Store the distance as output:
	W = setdat(a,-newout,fracrej);
else
   error('Illegal call to SVDD.');
end
return


