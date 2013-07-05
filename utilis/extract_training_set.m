%% This function extracts n training samples of size sz from Signal S.
% The input S can be either a 1-D signal or 2-D signal (image). For a 1-D
% signal n training sequences having sz samples will be extracted from
% random positions. For a 2-D Signal n training patches of size
% [sz(1),sz(2)] or [sz,sz] will be extracted from n random positions.
function [X, IdxSet] = extract_training_set(S, n, sz, IdxSet, op)

% if isscalar(sz) && min(size(S)) > 1
%     sz = [sz,sz];
% end
    sz_h = floor(sz./2);
    % Extraction window. awesomely coded :)
    Extractor = bsxfun(@plus,[-sz_h(1):(-sz_h(1)+sz(1)-1)]',(-sz_h(2):(-sz_h(2)+sz(2)-1))*size(S,1));
    if numel(sz) == 3
        Extractor = KronSum([0:numel(S)/sz(end):numel(S)-1],Extractor);
    end
    Extractor = Extractor(:);
    
    if ~isempty(op) && ~isscalar(op)
        op  = FullOp(op, sz(1), 1,1, size(S));
        res = sum(log(1+(op*double(S)).^2),3);
        res(1:sz_h(1),:)=0;
        res(end-sz_h(1)+1:end,:)=0;
        res(:,1:sz_h(2))=0;
        res(:,end-sz_h(2)+1:end)=0;
        res(IdxSet) = 0;
        [~,idx]=sort(res(:),'descend');
    else
        pos = reshape(1:numel(S(:,:,1)),size(S,1),size(S,2));
        pos(IdxSet)= 0; 
        pos(1:sz_h(1),:)=[];
        pos(end-sz_h(1)+1:end,:)=[];
        pos(:,1:sz_h(2))=[];
        pos(:,end-sz_h(2)+1:end)=[];
        
        if ~isempty(op) && op == 1
            idx = pos(1:sz(1):end,1:sz(2):end);
            idx = idx(:);
            n = numel(idx);
        else
            pos = pos(:);
            pos(pos==0)=[];
            sel = randperm(numel(pos));
            idx = pos(sel);
        end
    end
    C_pos = bsxfun(@plus,idx(1:n)',Extractor);
    X = S(C_pos);
    IdxSet = [IdxSet,idx(1:n)'];

end


% if ~isempty(op)
%     
%     sz_h = floor(sz./2);
%     
%     % Extraction window. awesomely coded :)
%     Extractor = bsxfun(@plus,[-sz_h(2):(-sz_h(2)+sz(2)-1)]',(-sz_h(1):(-sz_h(1)+sz(1)-1))*size(S,1));
%     Extractor = Extractor(:);
%     
%     op  = FullOp(op, sz(1), 1,1, size(S));
%     res = sum(abs(op*double(S)).^.2,3);
%     res(1:sz_h(1),:)=0;
%     res(end-sz_h(1)+1:end,:)=0;
%     res(:,1:sz_h(2))=0;
%     res(:,end-sz_h(2)+1:end)=0;
%     
%     [~,idx]=sort(res(:),'descend');
%     
%     C_pos = bsxfun(@plus,idx(1:n)',Extractor);
%     X = S(C_pos);
%    
% else
%     [i1,i2] = reggrid(size(S)-sz+1, n);
%     X = sampgrid(S, sz, i1, i2);
% end


%X = X(:,randperm(n));
%X = bsxfun(@minus,X,mean(X));
