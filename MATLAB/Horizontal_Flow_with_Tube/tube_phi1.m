function phi = tube_phi1(r)
    phi = zeros(4,4); % 4 by 4 matrix
    root = sqrt(1 + 4 * r - 4 * r^2);
    phi(1, :) = (1 + 2 * r - root) / 8; % phi(r - 2)
    phi(2, :) = (1 + 2 * r + root) / 8; % phi(r - 1)
    phi(3, :) = (3 - 2 * r + root) / 8; % phi(r)
    phi(4, :) = (3 - 2 * r - root) / 8; % phi(r + 1)
end
