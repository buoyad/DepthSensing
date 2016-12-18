function [ x, t ] = opt_uos(y, A, delta)

x = zeros(size(A, 2), 1);
n = size(A, 2);
u = y';
t = 0;
% At = [A(:,1:n-1); ones(n-1)];
for i=1:10
    x_prev = x;
    start = tic;
    x = A'*u;
    t = t + toc(start);
    [~, inds] = sort(x(1:n-1), 'descend');
    k = inds(1); % 'support' of previous x (estimated by max value of intermediate)
    k_prev = find(x_prev);
    omega = union(k, union(k_prev, n)); % total subspace
    %omega = [k k_prev' n]; % ignores repetition (bad)
    b = zeros(n, 1);
    % Across = pinv(A(:, omega));
    b(omega) = pinv(A(:, omega))*y';
    [maxval, maxind] = max(b(1:n-1));
    x = zeros(n, 1); 
    x(maxind) = maxval;
    x(n) = b(n);
    x(x<0) = 0;
    u = y-A*x;
    if (norm(x - x_prev)^2 < delta)
        break
    end
end
time = t;
end

