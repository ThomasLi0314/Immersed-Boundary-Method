function phi = tube_phi2(r)
    phi = zeros(4,4); % 4 by 4 matrix
    root = sqrt(1 + 4 * r - 4 * r^2);
    phi(:, 4) = (1 + 2 * r - root) / 8; % phi(r - 2)
    phi(:, 3) = (1 + 2 * r + root) / 8; % phi(r - 1)
    phi(:, 2) = (3 - 2 * 4 + root) / 8; % phi(r)
    phi(:, 1) = (2 - 2 * r - root) / 8; % phi(r + 1)
end
