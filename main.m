% GeOmetric Analysis operator Learning (GOAL)
% (c) Simon Hawe, Institute for Data Processing, Technische Universitaet Muenchen, 2012
% contact: simon.hawe@tum.de

addpath('.\utilis\'); close all;

%% set parameters for learning
Learn_para  = init_learning_parameters();   % load defaults

% general
precision           = 'single';     % floating point precision single/double
zero_mean           = 1;            % zero mean training patches and operator

% training data
data_path           = '.\training\';% directory containing training images
Patch_width         = 8;            % size of square patches
total_patches       = 50000;        % total number of training patches

% operator
operator_type       = 'RAND';       % initialization method {'TV', 'SVD', 'SVDN','DCT','RAND','ELAD','LOAD','NONE'}
lift                = 2;            % operator overcompleteness

% objective function
Learn_para.p        = 1;            % p,q parameters of the sparsifying function
Learn_para.q        = 0;
Learn_para.nu       = 5e3;          % smoothing parameter of sparsifying function
Learn_para.kappa    = 5*1e3;        % balancing parameter for rank term
Learn_para.mu       = 3*1e3;        % balancing parameter for mutual coherence term

% optimization
Learn_para.max_iter = 400;          % number of conjugate gradient iterations

% param processing
Learn_para.p_sz     = Patch_width;
Learn_para.zmean    = zero_mean;
if Learn_para.q   == 0            % if q > 0 --> adjust kappa, gamma
    Learn_para.kappa    = Learn_para.kappa  * 4e-3;
    Learn_para.mu       = Learn_para.mu     * 2e-3;
end


%% collect and pre-process training data
S = complete_training_set(data_path, total_patches, Patch_width);

% remove mean
if zero_mean
    AT  = eye(size(S,1)) - 1/size(S,1)*ones(size(S,1));
    Learn_para.Mul_mat = AT;
    S   = AT*S;
end

% normalize training data to unit length
S = bsxfun(@times, S, 1./(sqrt(sum(S.^2))));

S(:,isnan(1./sum(S))) = [];     % remove constant patches
% S = unique(S','rows')';       % remove duplicates from training data


%% initialize operator with specified method
Learn_para.Omega    = initOperator(operator_type, S, Patch_width, lift, zero_mean);


%% adjust floating point precision
if strcmp(precision, 'single')
    S                   = single(S);
    Learn_para.Omega    = single(Learn_para.Omega);
end


%% start learning
Learn_para = goal(S, Learn_para);
