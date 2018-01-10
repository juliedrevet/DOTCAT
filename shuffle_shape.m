function shape = shuffle_shape(expe)

sub_shape = shuffle_subshape;
shape_order = findmirror(expe);

% 2 last training shapes: not in the first two test blocks
pres = shuffle_subshape;
s = sub_shape(shape_order(5:end),:);
while sum(sum(ismember(pres(3:4,:),s(1:2,:))))~=0
    pres = shuffle_subshape;
end

% allocate shapes index
shape = [pres; s];

end

function subshape = shuffle_subshape

subshape = zeros(4,2);

k = 1;
shape_k = randperm(8,2);
while ~any(subshape(4,:))
    while any(ismember(shape_k,subshape))
        shape_k = randperm(8,2);
    end
    subshape(k,:) = shape_k;
    k = k+1;
end

end

