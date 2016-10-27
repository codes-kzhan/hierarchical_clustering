%% function hiCluster for RECOMB
% Author: Xiaoqian Wang
% \min_{W,Qi,F\in Ind} ||X'W+1b'-F||_F^2 + gamma*\sum_i^g||WQ_i||_{2,1}^2 + lambda*Tr(F'11'F)
function [outF, outQ, outW, outObj, outNumIter] = hiCluster(X, numC, numG, gamma, lambda, maxIter, inThresh)
% n: number of samples
% d: number of SNPs
% c: number of clusters
% X: d*n SNP data
% numC: number of clusters
% numG: number of groups of the clusters

% outF: n*c cluster indicator
% outQ: c*numG group indicator
% outW: d*c weight matrix

%% Initialization
%
if nargin < 7
    inThresh = 10^-6;
end
if nargin < 6
%    maxIter = 10^4;
    maxIter = 500;
end
if nargin < 5
    lambda = 1;
end
if nargin < 4
    gamma = 1;
end
delta = 10^-8;

%
c = numC;
[d, n] = size(X);

%
s = RandStream.create('mt19937ar','seed',7);  %% seed, "457" can be changed
RandStream.setGlobalStream(s);
A = rand(c, numG);
Q = A./(sum(A,2)*ones(1,numG));

%
tmp = randi(c,n,1);
F = dummyvar(tmp);

%
Xc = X - repmat(mean(X,2),1,n); % time: O(dn)
%Xc = X;
XX = Xc*Xc'; % time: O(d^2n)
for l = 1 : c
    XF(:,l) = sum(Xc(:,F(:,l)==1),2); % time: O(dc)
end
%XF = Xc*Fc; % time: O(dnc)
W = (XX + gamma*eye(d))\XF; % time: O(d^2c)
b = mean(F',2) - W'*mean(X,2);

%
D = getQL21Square(W, Q, delta); % time: O(dc^3)

%% Main code
obj = zeros(maxIter, 1);

for iter = 1: maxIter

    % fix D, W, update Q
    for g = 1 : numG
        A(:, g) = diag(W'*D{g}*W); % time: O(d^2c)
    end
    for t = 1 : c
        a = A(t, :)+eps;
        Q(t,:) = 1/sum(1./a)*1./a;
    end
    
    % fix W, b, update F
%     [~,tmp] = max(W'*X + b*ones(1,n));
%     %[~,tmp] = max(W'*X);
%     F = dummyvar(tmp);
%     if max(tmp) < c
%         F(:,max(tmp)+1:c) = zeros(n,c-max(tmp));
%     end
%     Fc = F - repmat(mean(F,1),n,1);
%     %Fc = F;
%     for l = 1 : c
%         XF(:,l) = sum(Xc(:,F(:,l)==1),2); % time: O(dc)
%     end
%     %XF = Xc*Fc; % time: O(dnc)
% 
%     F = solveF(X'*W + ones(n,1)*b');
%     F = solveF_svd(X'*W + ones(n,1)*b');
    [F, Fsum] = solveF_balance(X'*W + ones(n,1)*b', lambda);
    Fc = F - repmat(mean(F,1),n,1);
    XF = Xc*Fc; % time: O(dnc)
    
    % fix Q, D, F, update W and b
    for t = 1 : c
        tmpD = zeros(d, d);
        for g = 1: numG
            tmpD = tmpD + D{g}*(Q(t, g)^2);
        end
        W(:,t) = (XX + gamma*tmpD)\XF(:,t); % time: O(d^2n)
    end
    b = mean(F',2) - W'*mean(X,2);
    
    % fix W, Q, update D
    [D, Q21] = getQL21Square(W, Q, delta);

    % calculate obj
    Loss  =  norm(Fc - Xc'*W, 'fro')^2; % time: O(dnc)
    obj(iter) = Loss + gamma*Q21 + lambda*Fsum;

    if(iter > 1)
        if((obj(iter-1) - obj(iter))/obj(iter-1) < inThresh)
            break;
        end
    end

    if mod(iter, 10) == 0
        fprintf('process iteration %d, the obj is %d ...\n', iter, obj(iter));
    end
    
end

%% Outputs
%
[~, qMat] = max(Q, [], 2);
outQ = zeros(c,numG);
for t = 1:c     
    outQ(t, qMat(t)) = 1;
end

%
outNumIter = iter;

%
outObj = obj(1:iter);

%
[~, outF] = max(F');

%
outW = W;

end


% Di = diag(||WQ_i||_{2,1} / ||w^kQ_i||_2)
function [D, Q21] = getQL21Square(W, Q, delta)
% W: d*c weight matrix
% Q: c*numG group indicator

%% Initialization
%
numG = size(Q, 2);
Q21 = 0;

%% Main code
for g = 1:numG

    WQ = W*diag(Q(:,g));
    s = sum(WQ.*WQ,2).^0.5+delta;
    t = sum(s);
    D{g} = diag(t./s);
    Q21 = Q21 + t^2;
    
end

end