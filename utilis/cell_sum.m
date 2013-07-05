function Out = cell_sum(SC)
Out = zeros(size(SC{1}));
for i = 1:numel(SC)
    Out = Out + SC{i};
end
