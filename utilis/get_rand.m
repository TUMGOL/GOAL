function series = get_rand(range,nbr)
    series = randperm(range);
    if nargin == 2
        series = series(1:nbr)';
    end
end