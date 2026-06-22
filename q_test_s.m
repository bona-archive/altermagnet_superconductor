
%% q_test: s-wave
clc
% clear
% close all
st = 's_wave';

%% Parameters
Delta0 = 1e-5;
Nk = 70;
Nq = 201;
t = 1;

U = 3.5;
V = 1;
mu = - 1.6;
n_iter = 1e2;

%% Variables
qx_ran = linspace(-pi/2, pi/2, Nk+1);
qx_ran(end) = [];
qy_ran = qx_ran;
[qx_ran, qy_ran] = meshgrid(qx_ran, qy_ran);
qx_ran = reshape(qx_ran, [numel(qx_ran), 1]);
qy_ran = reshape(qy_ran, [numel(qy_ran), 1]);
x = qx_ran;
y = qy_ran;

%% Plotting options
ms = 60;
xl = 'qx';
yl = 'qy';
zl1 = 'E_0(q)';
zl2 = 'Delta(q)';

%%
% J_ran = 0: .2: 2;
% parfor p = 1: length(J_ran)
J_ran = 1;
for p = 1
    J = J_ran(p);
    dn = fullfile('q_test', st, sprintf('J%.3f', J));
    [~, ~] = mkdir(dn);
    tl = sprintf('%s: t=%.3f, J=%.3f, U=%.3f, V=%.3f, mu=%.3f', st, t, J, U, V, mu);
    dest1 = fullfile(dn, 'E_0.png');
    dest2 = fullfile(dn, 'Delta.png');
%     if isfile(dest1)
%         continue
%     end

    %% k-meshes
    kx = linspace(-pi, pi, Nk+1);
    kx(end) = [];
    ky = linspace(-pi, pi, Nk+1);
    ky(end) = [];
    [kx, ky] = meshgrid(kx, ky);
    kx = reshape(kx, [numel(kx), 1]);
    ky = reshape(ky, [numel(ky), 1]);
    xik = - 2 * t * cos(kx) - 2 * t * cos(ky);              % xi(k)
    Jk = J * sin(kx) .* sin(ky);                            % J(k)

    %% Square harmonics
    g_s_k = ones(size(kx));
    g_es_k = cos(kx) + cos(ky);
    g_pp_k = sin(kx) + 1i * sin(ky);
    g_pm_k = sin(kx) - 1i * sin(ky);
    g_d_k = cos(kx) - cos(ky);

    %% Specify SC state
    gamma = U;
    gk = g_s_k;

    tic
    %% Solve gap eq
    Delta_lst = zeros(size(qx_ran));
    E_0_lst = zeros(size(qx_ran));
    for i = 1: length(qx_ran)
        kx1 = kx + qx_ran(i);
        kx2 = kx - qx_ran(i);
        ky1 = ky + qy_ran(i);
        ky2 = ky - qy_ran(i);

        xik1 = - 2 * t * cos(kx1) - 2 * t * cos(ky1) - mu;  % xi(k+q)
        xik2 = - 2 * t * cos(kx2) - 2 * t * cos(ky2) - mu;  % xi(k-q)
        Jk1 = J * sin(kx1) .* sin(ky1);                     % J(k+q)
        Jk2 = J * sin(kx2) .* sin(ky2);                     % J(k-q)

        Xi = .5 * (xik1 + Jk1 + xik2 - Jk2);                % Xi(k;q)
        La = .5 * (xik1 + Jk1 - xik2 + Jk2);                % Lambda(k;q)

        Delta = Delta0;

        % solve gap equation
        for iter = 1: n_iter
            Ek = sqrt(Xi.^2 + (Delta * gk).^2);             % E_eta(k;q)
            Ek_u = Ek + La;                                 % E_u^eta(k;q)
            Ek_d = Ek - La;                                 % E_d^eta(k;q)
            ind_u = Ek_u > 0;
            ind_d = Ek_d > 0;
            Delta_k = gamma * (Delta * gk) / 2 ./ Ek .* (ind_u - (~ind_d));
            Delta = mean(Delta_k, 'omitnan');               % Delta_eta(q)
        end
        % compute ground state energy
        E_0_k = Xi - Ek + (Delta * gk).^2 ./ Ek / 2 .* (ind_u - (~ind_d))...
            + Ek_u .* (~ind_u) + Ek_d .* (~ind_d);
        E_0 = mean(E_0_k);                          % E_0^eta(q)

        % archive
        Delta_lst(i) = Delta;
        E_0_lst(i) = E_0;
    end
    toc

    %% Plotting
    z1 = E_0_lst;
    z2 = Delta_lst;
    figure;
    scatter3(x, y, z1, ms, z1, 'filled')
    colormap(gca);
    colorbar
    xlabel(xl)
    ylabel(yl)
    zlabel(zl1)
    view(2)
    title(tl)
    saveas(gcf, dest1)

    figure;
    scatter3(x, y, z2, ms, z2, 'filled')
    colormap(gca);
    colorbar
    xlabel(xl)
    ylabel(yl)
    zlabel(zl2)
    view(2)
    title(tl)
    saveas(gcf, dest2)
end