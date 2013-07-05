%(c) Simon Hawe, Lehrstuhl fuer Datenverarbeitung Technische Universitaet
%Muenchen, 2012. Contact: simon.hawe@tum.de
addpath('.\utilis\'); close all;
% Currently possible selections
Operators    = {'TV', 'SVD', 'SVDN','DCT','RAND','ELAD','LOAD','NONE'};
IName = 'couple.png';
%IName = 'barbara.png';
Patch_width     = 8;

%Rot_M = create_circle(Patch_width);

%q=fspecial('gaussian',Patch_width,(Patch_width)/3);
%Rot_M =diag(q(:));
% Rot_M = Rot_M -mean(Rot_M(:));
% Rot_M = Rot_M./max(Rot_M(:));
% Do we have to learn?
do_learning     = 1;
normalizer      = 1;

% Name of initial analysis operator, can be extended arbitarilly
% TV = 1; SVD = 2; SVDN = 3; DCT = 4; RAND = 5; ELAD = 6; LOAD = 7; NONE = 8
Operator_type = {Operators{5}};%,Operators{4}};
lift                = 2;
%S                   = [];
IdxSet              = {};
Omega = [];
Rot_M = 1;
%% First we start here with the learning phase if required
if do_learning || ~exist('Learn_para','var')
    
    % Get the images from which we learn
    clear Images;
    
    % Number of patches per image
    n_patches    = 35000;
    Images{1}    = '.\training\';
    Images{2}    = IName;
    
    % S containes the entire trainings set
    [S] =  complete_training_set(Images, n_patches, Patch_width);
    
    AT = eye(Patch_width^2)-1/Patch_width^2*ones(Patch_width^2);
    S =AT*S;
    S = bsxfun(@times,S,1./(sqrt(sum(S.^2))));
    S(:,isnan(1./sum(S)))=[];
    S = unique(S','rows')';
    
    Omega = [];
    for i=1:numel(Operator_type)
        % Selecting some type of initial analysis operator
        switch Operator_type{i}
            case 'TV'
                O = create_tvmat(Patch_width);
            case 'SVD'
                [O,~]=eig(S*S');
                O = create_lifting(lift,O);
            case 'SVDN'
                [O,vals]=eig(S*S');
                O = bsxfun(@times,O,1./diag(vals));
                O = bsxfun(@times,O,1./sqrt(sum(O.^2,1)));
                O = create_lifting(lift,O);
            case 'DCT'
                O = kron(dctmtx(Patch_width),dctmtx(Patch_width));
                O = create_lifting(lift,O);
            case 'RAND'
                randn('state',0);
                O = randn(round(lift*size(S,1)),size(S,1));
                %[O,~,~]=svd(O);
                %O = create_lifting(lift,O);
            case 'ELAD'
                O = initOmega(S,round(lift*size(S,1)));
            case 'LOAD'
                O = load(sprintf('%dx%dSquared.mat',Patch_width,Patch_width));
                O = create_lifting(lift,O.Omega);
            case 'NONE'
                O = Learn_para.Omega;
        end
        Omega = [Omega;O];
    end
    Learn_para  = init_learning_parameters();
    Omega       = bsxfun(@times,Omega,1./sqrt(sum(Omega.^2,2)));
    
    % Initialize the parameters required for learning things
    Learn_para.Omega    = Omega;
    Learn_para.max_iter = 400;
    
    Learn_para.kappa    = 1*1e4;
    Learn_para.mu       = 3*1e0;
    
    Learn_para.p        = 1;
    Learn_para.q        = 0;
    
    
    Learn_para.p_sz    = Patch_width;
    
    if Learn_para.q   == 0
        Learn_para.kappa    = Learn_para.kappa*4e-3;
        Learn_para.mu    = Learn_para.mu*2e-3;
    end
    
    Learn_para.nu       = 1e3;
    
    Learn_para.Mul_mat = AT;
    
    Learn_para = goal(S, Learn_para);
end












