% Navier Stokes Solver for pIB

function [u_final, u_intermediate] = navier_stokes_solver_pIB(u, Force, a, el_fx, el_bx, el_fy, el_by, dx, dy, dt, mu, rho)
    %% Prelimnayr step
    w_pre = u - dt / 2 * skew_symmetric_advection(u, el_fx, el_bx, el_fy, el_by, dx, dy) + dt / (2 * rho) * Force;

    % Discrete Fourier transform
    w_pre = fft(w_pre, [], 1);
    w_pre = fft(w_pre, [], 2);

    % Update the fourier transform of u
    u_temp(:,:,1)=a(:,:,1,1).*w_pre(:,:,1)+a(:,:,1,2).*w_pre(:,:,2);
    u_temp(:,:,2)=a(:,:,2,1).*w_pre(:,:,1)+a(:,:,2,2).*w_pre(:,:,2);
    u_real_temp=ifft(u_temp,[],2);
    u_intermediate = real(ifft(u_real_temp, [], 1));

    %% FInal Step
    w_final = u - dt * skew_symmetric_advection(u_intermediate, el_fx, el_bx, el_fy, el_by, dx, dy) + (dt / 2) * (mu / rho) * laplacian_tube(u, el_fx, el_bx, el_fy, el_by, dx) + dt / rho * Force;
    w_final = fft(w_final, [], 1);
    w_final = fft(w_final, [], 2);

    % Update the velocity
    u_final_temp(:, :, 1) = a(:,:,1,1).*w_final(:,:,1)+a(:,:,1,2).*w_final(:,:,2);
    u_final_temp(:, :, 2) = a(:,:,2,1).*w_final(:,:,1)+a(:,:,2,2).*w_final(:,:,2);
    u_real_final = ifft(u_final_temp, [], 1);
    u_final = real(ifft(u_real_final, [], 2));
end