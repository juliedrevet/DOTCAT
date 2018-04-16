function shape = shuffle_shape(expe)
% only works for 2 practice + 8 testing blocks


shape = zeros(length(expe.blck),2);
nprac = expe.cfg.nprac;
sub_shape = shuffle_subshape;
% condtn always alternating 1 2 1 2 1 2 1 2
% apply sub_shape to all condtn(1)
shape((nprac+1):2:end,:) = sub_shape;
% apply sub_shape to all mirror
% mirror from 1:4 is mod(((1:4)+1),4)+1
shape((nprac+2):2:end,:) = sub_shape(mod(((1:4)+1),4)+1,:);

% apply shape for practice blocks (should differ from 2 first blocks)
prac_shape = shape((nprac+3):(nprac+4),:);

shape(1:nprac,:) = prac_shape';
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

