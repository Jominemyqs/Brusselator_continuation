function out = objective_evaluation(feature1, feature2)
% Outputs distance between feature collections

out = 0;

type = checkType(feature1);

switch type
    case 'scalar'
        out = abs(feature1 - feature2);

    case 'vector'
        out = ws_distance(feature1, feature2, 2);

    case 'cell'
        path(path,'mexEMD/');
        Rep = length(feature1);
        dist_mat = zeros(Rep, Rep);
        for i = 1:Rep
            for j = 1:Rep
                dist_mat(i,j) = ws_distance(feature1{i}, feature2{j}, 2);
            end
        end
        out = mexEMD(ones(Rep,1)/Rep, ones(Rep,1)/Rep, dist_mat);

    case 'matrix'
        % Not expected in your pipeline; fail loudly so you notice immediately.
        error('objective_evaluation:MatrixInput', ...
              'Received a matrix feature; expected scalar/vector/cell.');

    otherwise
        error('objective_evaluation:UnknownType', 'Unknown feature type.');
end

end


function type = checkType(input)

if iscell(input)
    type = 'cell';
    return
end

if ~isnumeric(input)
    type = 'unknown';
    return
end

if isscalar(input)
    type = 'scalar';
elseif isvector(input)
    type = 'vector';
else
    type = 'matrix';
end

end

