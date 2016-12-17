function [ x ] = opt_uos(y, A, delta)

x = zeros(size(A, 2), 1);
u = y';
for i=1:10
    x_prev = x;
    x = A'*u;
    n = size(A, 2);
    [~, inds] = sort(x(1:n-1), 'descend');
    k = inds(1); % 'support' of previous x (only one nonzero entry besides B)
    k_prev = find(x_prev);
    omega = union(k, union(k_prev, n)); % total subspace
    b = zeros(n, 1);
    % Across = pinv(A(:, omega));
    b(omega) = pinv(A(:, omega))*y';
    [maxval, maxind] = max(b);
    x = zeros(n, 1); 
    x(maxind) = maxval;
    x(n) = b(n);
    u = y-A*x;
    if (norm(x - x_prev)^2 < delta)
        break
    end
end
end

