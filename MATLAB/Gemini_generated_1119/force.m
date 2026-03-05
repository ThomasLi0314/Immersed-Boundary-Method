function F = force(X, K, index_forward, index_backward, dtheta)
    F = K * (X(index_forward, :) + X(index_backward, :) - 2 * X) / (dtheta ^2);
end