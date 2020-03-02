function [x_pr, val] = phase_retrieval(y, A, n_iter, step_size)
%PHASE_RETRIEVAL solves a phase retrieval problem
% The forward model is y = abs(A * x).^2
% From y and A, it returns the result of the PR algorithm

% General parameters
[n, d] = size(A);
alpha = n / d;

%% Spectral initialization
% Normalization of the intensity
y = y / sqrt(var(y));
t = (y-1) ./ (y+sqrt(2*alpha)-1+1e-6);  % optimal weights
Z = A' * diag(t) * A;
[V, D] = eig(Z);
first_eigenvector = V(:, 1);

%% Gradient descent iterations
curr_estimate = first_eigenvector;
for i_iter = 1:n_iter
    % Compute current estimation of y
    field = A * curr_estimate;
    curr_y = abs(field).^2;
    if i_iter == 1
        zero_err = sum((curr_y-y).^2);
    end
    val(i_iter) = sum((curr_y-y).^2)/zero_err;
    % Compute gradient
    grad = A' * ((curr_y-y) .* field);
    % Gradient descent update
    curr_estimate = curr_estimate - step_size * grad;
end
x_pr = curr_estimate;

end

