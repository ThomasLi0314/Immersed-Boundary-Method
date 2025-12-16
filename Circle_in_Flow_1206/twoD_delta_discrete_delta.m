% Delta function grid

% Updated Dec 13th
function w = twoD_delta_discrete_delta(x, y)
    root_x = sqrt(1 + 4 * x - 4 * x ^2);
    root_y = sqrt(1 + 4 * y - 4 * y ^2);

    phi1(4, :) = (1 + 2 * x - root_x) / 8;
    phi1(3, :) = (1 + 2 * x + root_x) / 8;
    phi1(2, :) = (3 - 2 * x + root_x) / 8;
    phi1(1, :) = (3 - 2 * x - root_x) / 8;

    phi2(:, 4) = (1 + 2 * y - root_y) / 8;
    phi2(:, 3) = (1 + 2 * y + root_y) / 8;
    phi2(:, 2) = (3 - 2 * y + root_y) / 8;
    phi2(:, 1) = (3 - 2 * y - root_y) / 8;

    w = phi1 .* phi2; % Compute the outer product of phi1 and phi2

    % see notes, same row has same x coordinate. 
end

    