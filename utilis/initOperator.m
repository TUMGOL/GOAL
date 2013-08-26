function [ operator ] = initOperator( Operator_type, S, p_sz, lift, z_mean )
%INITOPERATOR Initializes operator according to specific method
%   Detailed explanation goes here

% Currently possible selections
% Operators    = {'TV', 'SVD', 'SVDN','DCT','RAND','ELAD','LOAD','NONE'};

% Selecting some type of initial analysis operator
switch Operator_type
    case 'TV'
        operator = create_tvmat(p_sz);
    case 'SVD'
        [operator,~]=eig(S*S');
        operator = create_lifting(lift,operator);
    case 'SVDN'
        [operator,vals]=eig(S*S');
        operator = bsxfun(@times,operator,1./diag(vals));
        operator = bsxfun(@times,operator,1./sqrt(sum(operator.^2,1)));
        operator = create_lifting(lift,operator);
    case 'DCT'
        operator = kron(dctmtx(Patch_width),dctmtx(p_sz));
        operator = create_lifting(lift,operator);
    case 'RAND'
        randn('state',0);
        operator = randn(round(lift*size(S,1)),size(S,1));
    case 'ELAD'
        operator = initOmega(S,round(lift*size(S,1)));
    case 'LOAD'
        operator = load(sprintf('%dx%dSquared.mat',p_sz,p_sz));
        operator = create_lifting(lift,operator.Omega);
%         load('C:\Users\eMKay\Desktop\Omega.mat');
%         operator = Omega;
    case 'NONE'
        operator = Learn_para.Omega;
end

if z_mean
    % project analysis operator atoms onto 1-normal hyperplane
    operator = bsxfun(@minus, operator, mean(operator,2));
end

% normalize rows of operator to unit length
operator = bsxfun(@times, operator, 1./sqrt(sum(operator.^2, 2)));

if z_mean
    % project analysis operator atoms onto 1-normal hyperplane
    operator = bsxfun(@minus, operator, mean(operator,2));
end

end

