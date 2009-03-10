function p = mirhmm(x,varargin)
%   p = mirhmm(m) describes the signal x as the most probable sequence of
%       states of a Hidden Markov Model (HMM).
%   Optional argument: 
%       mirhmm(...,'States',Q) specifies the number of hidden states.
%           default value: Q = 3
%       mirhmm(...,'Distrib',M) specifies the number of distributions
%           in the Gaussian mixture associated with each state.
%           default value: M = 10
%       mirhmm(...,'MaxIter',mi) specifies the maximal number of EM
%           iterations.
%           default value: mi = 20
%       mirhmm(...,c) where c = 'full', 'diag' (default value) or 'spherical'
%           Covering type
%
% From "How to use the HMM toolbox"
% http://www.cs.ubc.ca/~murphyk/Software/HMM/hmm_usage.html

[Q,M,c,mi] = scanargin(varargin);
cc = mirmfcc(mirframe(x,'WinLength',.03,'Hop',.5),10);
obs = get(cc,'Data');
for k = 1:length(obs)
    data = obs{k};
    O = size(data,1);
    % First, let us fit a mixture of M Gaussians
    % for each of the Q states using K-means.
    try
        prior0 = normalise(rand(Q,1));
    catch
        error('Please install HMM toolbox.');
    end
    transmat0 = mk_stochastic(rand(Q,Q));
    [mu0, Sigma0] = mixgauss_init(Q*M,data,c);
    mu0 = reshape(mu0, [O Q M]);
    Sigma0 = reshape(Sigma0, [O O Q M]);
    mixmat0 = mk_stochastic(rand(Q,M));

    % Then let us improve these parameter estimates using EM.
    [LL, prior1, transmat1, mu1, Sigma1, mixmat1] = ...
     mhmm_em(data, prior0, transmat0, mu0, Sigma0, mixmat0, ...
                                                'max_iter', mi);
    % Since EM only finds a local optimum, good initialisation is crucial. 
    % The initialisation procedure illustrated above is very crude,
    % and is probably not adequate for real applications... 
    % Use trainHMM for a real-world example of EM 
    % with mixtures of Gaussians using BNT.

    % Computing the most probable sequence (Viterbi)
    % First we need to evaluate B(i,t) = P(y_t | Q_t=i) for all t,i
    B = mixgauss_prob(data, mu1, Sigma1, mixmat1);
    % Finally,
    path{k} = viterbi_path(prior1, transmat1, B);
end    
p = mirscalar(cc,'Data',path,'Title','HMM path');    


function [Q,M,c,mi] = scanargin(v)
Q = 3;
M = 10;
c = 'diag';
mi = 20;
if exist('v')
    i = 1;
    while i <= length(v)
        arg = v{i};
        %if isnumeric(arg)
        %else
        if strcmpi(arg,'States')
            i = i+1;
            Q = v{i};
        elseif strcmpi(arg,'Distrib')
            i = i+1;
            M = v{i};
        elseif strcmpi(arg,'MaxIter')
            i = i+1;
            mi = v{i};
        elseif strcmpi(arg,'full') || strcmpi(arg,'diag') || ...
                                     strcmpi(arg,'spherical')
            c = arg;
        else
            error('ERROR IN HMM: Syntax error. See help hmm.');
        end
        i = i+1;
    end
end


function h = trainHMM(O,N,M)
% Based on the Short Tutorial "Creating a Gaussian Mixture Model using BNT"
% by Richard W. DeVaul, Version 1.0
% http://www.media.mit.edu/wearables/mithril/BNT/mixtureBNT.txt
% PROBLEM: IT NEEDS A SUPERVISED TRAINING!

maxiter=10;     %% The number of iterations of EM (max)
epsilon=1e-100; %% A very small stopping criterion
dag = [ 0 1 1 ; 0 0 1 ; 0 0 0 ];
discrete_nodes = [1 2];
nodes = [1 : 3];
for k = 1:length(O)
    Ok = O{k};
    node_sizes=[ N M size(Ok,1) ];
    % The model is a N-class, M component mixture model
    bnet = mk_bnet(dag, node_sizes, 'discrete', discrete_nodes);
    bnet.CPD{1} = tabular_CPD(bnet,1);
    bnet.CPD{2} = tabular_CPD(bnet,2);
    bnet.CPD{3} = gaussian_CPD(bnet, 3);
    
    % First, we create the (exact) inference engine 
    % that will allow EM to estimate the model parameters.
    engine = jtree_inf_engine(bnet);
    
    % Next, we fit the model.
    [bnet2, ll, engine2] = learn_params_em(engine,training,maxiter,epsilon);
end