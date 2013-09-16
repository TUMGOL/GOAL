% GeOmetric Analysis operator Learning (GOAL)
% (c) Simon Hawe
% Institute for Data Processing, Technische Universitaet Muenchen, 2012
% contact: simon.hawe@tum.de
% http://www.gol.ei.tum.de

addpath('.\utilis\'); close all;

%% set parameters for learning
% initialize with default parameters
para            = init_learning_parameters();   % load defaults

% general
precision       = 'single';     % floating point precision 'single'/'double'
para.zmean      = 1;            % zero mean training patches and operator

% training data
data_path       = '.\training\';% directory containing training images
para.p_sz       = 8;            % width of square patches
total_patches   = 50000;        % total number of training patches

% operator
operator_type   = 'RAND';       % initialization method {'TV', 'SVD', 'SVDN','DCT','RAND','ELAD','LOAD','NONE'}
lift            = 2;            % operator overcompleteness

% objective function
para.Sp_type    = 'LogSquare';  % type of sparsifying function
para.p          = 1;            % p,q parameters of the sparsifying function
para.q          = 0;
para.nu         = 5e3;          % smoothing parameter of sparsifying function
para.kappa      = 5*1e3;        % balancing parameter for rank term
para.mu         = 3*1e3;        % balancing parameter for mutual coherence term

% optimization
para.max_iter   = 400;          % number of conjugate gradient iterations

% para processing
if para.q   == 0                % if q > 0 --> adjust kappa, gamma
    para.kappa    = para.kappa  * 4e-3;
    para.mu       = para.mu     * 2e-3;
end


%% collect and pre-process training data
S = complete_training_set(data_path, total_patches, para.p_sz);

% remove mean
if para.zmean
    AT  = eye(size(S,1)) - 1/size(S,1)*ones(size(S,1));
    para.Mul_mat = AT;
    S   = AT*S;
end

% normalize training data to unit length
S = bsxfun(@times, S, 1./(sqrt(sum(S.^2))));

S(:,isnan(1./sum(S))) = [];     % remove constant patches
% S = unique(S','rows')';       % remove duplicates from training data


%% initialize operator with specified method
para.Omega    = initOperator(operator_type, S, para.p_sz, lift, para.zmean);


%% adjust floating point precision
if strcmp(precision, 'single')
    S                   = single(S);
    para.Omega    = single(para.Omega);
end


%% start learning
para = goal(S, para);
