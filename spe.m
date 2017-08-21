function Y = spe(X, no_dims, varargin)
%SPE Perform the Stochastic Proximity Embedding algorithm
%
%   Y = spe(X, no_dims, 'Global')
%   Y = spe(X, no_dims, 'Local', k)

if nargin <= 2
    variant = 'Global';
else
    variant = varargin{1};
end
if strcmp(variant, 'Local')
    if nargin > 3, k = varargin{2};
    else k = 12; end
end
if ~strcmp(variant, 'Global') && ~strcmp(variant, 'Local')
    error('Unknown parameter.');
end
if exist('k', 'var') && ischar(k)
    error('Adaptive neighborhood selection is not yet supported in SPE.');
end

% Initialize parameters
lambda = 1;                                         % initial learning parameter
s = 100;                                            % number of updates per iteration
%     max_iter = 20000 + round(.04 * size(X, 1) ^ 2);     % number of iterations
max_iter = 20000 + round(.0004 * size(X, 1) ^ 2);     % number of iterations
tol = 1e-5;                                         % regularlization parameter
n = size(X, 1);                                     % number of datapoints
if strcmp(variant, 'Local')
    max_iter = max_iter * 3;
end

% Compute proximity matrix in original space
if strcmp(variant, 'Global')
    R = pdist(X', X');
    R = R / max(max(R)) * sqrt(2);
else
    [R, n_ind] = find_nn(X, k);
end

% Initialize datapoints randomly
Y = rand(n, no_dims);

% Perform SPE
for i=1:max_iter
    if rem(i, 10000) == 0
        disp(['Iteration ' num2str(i) ' of ' num2str(max_iter) '...']);
    end
    
    % Select points that should be updated
    J = randperm(n);
    ind1 = J(1:s);
    if strcmp(variant, 'Global')
        ind2 = J(s+1:2*s);
    else
        ind2 = double(n_ind(ind1,:))';
        J = round(rand(1, size(ind2, 2)) * (k - 1)) + 1;
        J = J + ((0:length(J) - 1) * k);
        ind2 = ind2(J);
    end
    
    % Compute distances between points in embedded space
    D = sqrt(sum((Y(ind1,:) - Y(ind2,:)) .^ 2, 2));
    
    % Get corresponding distances in real space
    Rt = R((ind1 - 1) * size(R, 1) + ind2)';
    
    % Update locations of points
    Y(ind1,:) = Y(ind1,:) + lambda * (1/2) * repmat(((Rt - D) ./ (D + tol)), 1, no_dims) .* (Y(ind1,:) - Y(ind2,:));
    Y(ind2,:) = Y(ind2,:) + lambda * (1/2) * repmat(((Rt - D) ./ (D + tol)), 1, no_dims) .* (Y(ind2,:) - Y(ind1,:));
    
    % Update learning parameter
    lambda = lambda - (lambda / max_iter);
end



function [D, ni] = find_nn(X, k)
%FIND_NN Finds k nearest neigbors for all datapoints in the dataset
%
%	[D, ni] = find_nn(X, k)
%

if ~exist('k', 'var') || isempty(k)
    k = 12;
end

% Perform adaptive neighborhood selection if desired
if ischar(k)
    [D, max_k] = find_nn_adaptive(X);
    ni = zeros(size(X, 1), max_k);
    for i=1:size(X, 1)
        tmp = find(D(i,:) ~= 0);
        tmp(tmp == i) = [];
        tmp = [tmp(2:end) zeros(1, max_k - length(tmp) + 1)];
        ni(i,:) = tmp;
    end
    
    % Perform normal neighborhood selection
else
    
    % Compute distances in batches
    n = size(X, 1);
    sum_X = sum(X .^ 2, 2);
    batch_size = round(2e7 ./ n);
    D = zeros(n, k);
    ni = zeros(n, k);
    for i=1:batch_size:n
        batch_ind = i:min(i + batch_size - 1, n);
        DD = bsxfun(@plus, sum_X', bsxfun(@plus, sum_X(batch_ind), ...
            -2 * (X(batch_ind,:) * X')));
        [DD, ind] = sort(abs(DD), 2, 'ascend');
        D(batch_ind,:) = sqrt(DD(:,2:k + 1));
        ni(batch_ind,:) = ind(:,2:k + 1);
    end
    D(D == 0) = 1e-9;
    D = sparse([repmat(1:n, [1 k])'; ni(:)], [ni(:); repmat(1:n, [1 k])'], [D(:); D(:)], n, n);
end