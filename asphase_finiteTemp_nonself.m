% clear all

ParamN = 50;
t = 0.5;
mu = -1;
V = 1.5;
J_ran = linspace(0, 1, ParamN);
T_ran = linspace(10^(-4), 0.01, 100);
% 0 ~ 0.02 구간: 밀도가 높은 부분
DN = 5000;
D1end = 0.01;
D2end = 0.025;
D1 = linspace(0, D1end, DN);  % 예: 50개의 샘플

% 0.02 ~ 0.3 구간: 밀도가 낮은 부분
D2 = linspace(D1end, D2end, DN+1); % 예: 50개의 샘플

% 두 구간 병합
Dinput = [D1, D2(2:end)]; % 중복되는 값 제거


kx = linspace(-pi, pi, N+1);
ky = linspace(-pi, pi, N+1);
kx(end) = [];
ky(end) = [];
qx_scan = kx;

[kx, ky] = meshgrid(kx, ky);
kx = reshape(kx, [numel(kx), 1]);
ky = reshape(ky, [numel(ky), 1]);

epsilon = 1e-6;

gpk = sin(kx) + 1i * sin(ky);
gpka = abs(gpk) .^ 2;
gesk = cos(kx) + cos(ky);
geska = abs(gesk) .^ 2;
gdk = cos(kx) - cos(ky);
gdka = abs(gdk) .^ 2;
ek = - 2* t * cos(kx) - 2 * t * cos(ky);

Dout1 = zeros(1, length(Dinput));
E0out1 = zeros(1, length(Dinput));
Dout2 = zeros(1, length(Dinput));
E0out2 = zeros(1, length(Dinput));
Dout3 = zeros(1, length(Dinput));
E0out3 = zeros(1, length(Dinput));
Fq = zeros(1, 3);
rho = zeros(1, length(T_ran));

%% s-wave
disp('start s-wave')

dqx = qx_scan(2)-qx_scan(1);
jjj=11;
% jjj = jjff;

for j = 1: length(T_ran)
    j
    T = T_ran(j);
    J = J_ran(jjj);
    % Delta = Dff_self(jjj);
    qx = q_self(jjj);



        Delta = 0.22;



        [Xi, La] = XiLaMaker(kx, ky, qx-dqx, 0, t, J, mu);
        [~, Da1, ~, Ea1] = FFgapsolver_finiteTemp(Xi, La, Delta, V, T);

    Fq(1) = Ea1;

    % [Xi, La] = XiLaMaker(kx, ky, qx-dqx, 0, t, J, mu);
    % [Ea1, ~, ~] = FFgapsolver_finiteTemp(Xi, La, Delta, V, T);
    % Fq(1) = Ea1;

        [Xi, La] = XiLaMaker(kx, ky, qx, 0, t, J, mu);
        [~, Da2, ~ ,Ea2] = FFgapsolver_finiteTemp(Xi, La, Delta, V, T);

    Fq(2) = Ea2;

    % [Xi, La] = XiLaMaker(kx, ky, qx, 0, t, J, mu);
    % [Ea2, Da2, ~] = FFgapsolver_finiteTemp(Xi, La, Delta, V, T);
    % Fq(2) = Ea2;

        [Xi, La] = XiLaMaker(kx, ky, qx+dqx, 0, t, J, mu);
        [~, Da3, ~, Ea3] = FFgapsolver_finiteTemp(Xi, La, Delta, V, T);

    Fq(3) = Ea3;
    % [Xi, La] = XiLaMaker(kx, ky, qx+dqx, 0, t, J, mu);
    % [Ea3, ~, ~] = FFgapsolver_finiteTemp(Xi, La, Delta, V, T);
    % Fq(3) = Ea3;

    rho(j) = (Fq(3) - 2 * Fq(2) + Fq(1)) / dqx / dqx;
end

lam = 1 ./ sqrt(abs(rho));


figure
plot(T_ran, rho/rho(1));xlabel('T [K]');ylabel('K(0,T)/K(0, 0)')


figure
plot(T_ran, lam / lam(1));xlabel('T [K]');ylabel('\lambda(T)/\lambda(0)')
%% function
function [Xi,La] = XiLaMaker(kx, ky, qx, qy, t, J, mu)

kx1 = kx + qx;
kx2 = -kx + qx;
ky1 = ky + qy;
ky2 = -ky + qy;

Jk1 = J * sin(kx1) .* sin(ky1);
Jk2 = J * sin(kx2) .* sin(ky2);

xik1 = - 2 * t * cos(kx1) - 2 * t * cos(ky1) - mu;
xik2 = - 2 * t * cos(kx2) - 2 * t * cos(ky2) - mu;

Xi = .5 * (xik1 + Jk1 +xik2 - Jk2);                            
La = .5 * (xik1 + Jk1 -xik2 + Jk2);

end

function [E_0, Dffout,Ek, F] = FFgapsolver_finiteTemp(Xi, La, Delta, U, T)
Ek = sqrt(Xi .^2 + Delta^2);
Ek_u = Ek + La;
Ek_d = Ek - La;

Ek_merged = [Ek_u;Ek_d;-Ek_u;-Ek_d];
% ind_u = Ek_u > 0;
% ind_d = Ek_d > 0;
Delta_k = U * Delta /4 ./ Ek .* (tanh(Ek_u / 2 / T) + tanh(Ek_d / 2 / T));

Dffout = mean(Delta_k);

if T == 0
    F = Ek_u.*1./(1+exp(Ek_u/T))+Ek_d.*1./(1+exp(Ek_d/T))+(-Ek_u).*1./(1+exp(-Ek_u/T))+(-Ek_d).*1./(1+exp(-Ek_d/T));
    F = mean(F, 'omitnan');
else
    c = max(-Ek_merged/T);
    logZ = c + log(sum(exp(-Ek_merged/T-c)));
    F = -T * logZ;
end
% E_0_k = Xi - Ek + Delta^2 ./ Ek /4 .* (tanh(Ek_u / 2 / T) + tanh(Ek_d / 2 /T))...
%     + Ek_u.*1./(1+exp(Ek_u/T))+Ek_d.*1./(1+exp(Ek_d/T))+(-Ek_u).*1./(1+exp(-Ek_u/T))+(-Ek_d).*1./(1+exp(-Ek_d/T));
E_0_k = Ek_u.*1./(1+exp(Ek_u/T))+Ek_d.*1./(1+exp(Ek_d/T))+(-Ek_u).*1./(1+exp(-Ek_u/T))+(-Ek_d).*1./(1+exp(-Ek_d/T));
E_0 = mean(E_0_k, 'omitnan');

end