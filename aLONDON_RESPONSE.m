%% Parameters
t = .5; % Hopping parameter
Delta = 0.5; % Superconducting gap
mu = -1; % Chemical potential
J = 0.4;
e = 1.6e-19; % Elementary charge
hbar = 1; % Reduced Planck's constant
c = 1; % Speed of light
tolerance = 10^(-1.5); % Tolerance for selecting Fermi surface points, Debye cutoff
k_points = linspace(-pi, pi, 100); % Discretize k-space for integration (100x100 grid)

% Create a 2D k-space grid
[kx, ky] = meshgrid(k_points, k_points);

% Compute normal state energy xi(k) for each (kx, ky)
J_k = J .* sin(kx) .* sin(ky);
xi_kup = -2 * t * (cos(kx) + cos(ky)) - mu + J_k;
xi_kdown = -2 * t * (cos(kx) + cos(ky)) - mu - J_k;


% Select points on the Fermi surface (where xi_k is close to 0)
figure;hold on
[fermi_surface_maskup, kx_fsup, ky_fsup] = FSmasking(xi_kup, tolerance, kx, ky);
[fermi_surface_maskdown, kx_fsdown, ky_fsdown] = FSmasking(xi_kdown, tolerance, kx, ky);
%% 

% Initialize response tensor components
R_xx = 0;
R_yy = 0;
R_xy = 0;

jjj = 21; % gapless FF state
jjj = 11; % gapped s-wave

% I2의 온도 의존성 계산하기. self-consistent Delta(T) > Ekup, Ekdown 스펙트럼 > I2(T) 계산 
for j = 1: length(T_ran)
    j
    T = T_ran(j);
    J = J_ran(jjj);
    qx = q_self(jjj);
    for i = 1: length(Dinput)
        Delta = Dinput(i);

        J = J_ran(jjj);

        [Xi, La] = XiLaMaker(kx, ky, qx, 0, t, J, mu);
        [Ea1, Da1, ~, ~] = FFgapsolver_finiteTemp(Xi, La, Delta, V, T);

        Dout1(i) = Da1;
        E0out1(i) = Ea1;

    end

    [~ , I] = min(abs(Dout1 ./ (Dinput + epsilon) - 1));
    if (Dout1(2) - 0) / (Dinput(2) - Dinput(1)) < 1
        I = 1;
        Dsq1 = 0;
        E0sq1 = E0out1(I);
    else
        Dsq1 = Dinput(I);
        E0sq1 = E0out1(I);
    end
    Fq(1) = E0sq1;
    E_kup = sqrt(xi_kup.^2 + Delta.^2) + J_k;
    E_kdown = sqrt(xi_kdown.^2 + Delta.^2) - J_k;

end

% Loop over each (kx, ky) on the Fermi surface
for idx = 1:length(kx_fsup)
    
    [Rxx_out, Ryy_out, Rxy_out] = LONDON_RESPONSE_TENSOR_singlet(kx_fs, ky_fs, idx, t, mu, J, Delta, s, Rxx, Ryy, Rxy, T);
end

% Normalize by the number of points on the Fermi surface for averaging
num_fs_points = length(kx_fs);
R_xx = R_xx / num_fs_points;
R_yy = R_yy / num_fs_points;
R_xy = R_xy / num_fs_points;

% Display the results
fprintf('The computed response tensor component R_xx is: %.4e\n', R_xx);
fprintf('The computed response tensor component R_yy is: %.4e\n', R_yy);
fprintf('The computed response tensor component R_xy is: %.4e\n', R_xy);

%% function
function [fermi_surface_mask, kx_fs, ky_fs] = FSmasking(xi_k, wDebye, kx, ky)
fermi_surface_mask = abs(xi_k) < wDebye;
kx_fs = kx(fermi_surface_mask);
ky_fs = ky(fermi_surface_mask);
scatter(kx_fs,ky_fs)
end


function [Rxx_out, Ryy_out, Rxy_out] = LONDON_RESPONSE_TENSOR_singlet(kx_fs, ky_fs, idx, t, mu, J, Delta, s, Rxx, Ryy, Rxy, T)
    % Compute normal state energy for this point
    xi_k_point = -2 * t * (cos(kx_fs(idx)) + cos(ky_fs(idx))) - mu;
    J_k_point = J * sin(kx_fs(idx)) * sin(ky_fs(idx));
    % Compute quasiparticle energy for superconducting state
    E_k = sqrt(xi_k_point^2 + Delta^2) + s * J_k_point;
    
    % Fermi velocity components for the normal state
    v_Fx = 2 * t * sin(kx_fs(idx)) + s * J * cos(kx_fs(idx)) * sin(ky_fs(idx));
    v_Fy = 2 * t * sin(ky_fs(idx)) + s * J * sin(kx_fs(idx)) * cos(ky_fs(idx));
    
    % Energy integral part
    f_prime = -(1/T)*exp(-E_k/T) ./ (1 + exp(-E_k/T)).^2; % Approximation for dF/dE at low T
    N_E = E_values ./ sqrt(E_values.^2 - Delta^2); % Density of states approximation
    
    % Contribution to R_xx, R_yy, and R_xy
    Rxx_out = Rxx + (e^2 / (4 * pi^3 * hbar * c)) * (v_Fx^2 / abs(v_Fx)) * energy_integral;
    Ryy_out = Ryy + (e^2 / (4 * pi^3 * hbar * c)) * (v_Fy^2 / abs(v_Fy)) * energy_integral;
    Rxy_out = Rxy + (e^2 / (4 * pi^3 * hbar * c)) * (v_Fx * v_Fy / sqrt(v_Fx^2 + v_Fy^2)) * energy_integral;

end