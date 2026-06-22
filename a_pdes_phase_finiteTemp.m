clear all

ParamN = 50;
N=1000;
t = 0.5;
mu = -1;
V_ran = linspace(0.6, 1.5, ParamN);
V = V_ran(23);
V =1;
J_ran = linspace(0, 1, ParamN);
epsilon = 1e-6;
T_ran = linspace(10^(-3), 0.009, 20);
T_ran = linspace(10^(-3), 4, 5);
% 0 ~ 0.02 구간: 밀도가 높은 부분
DN = 500;
D1end = 0.012;
D2end = 0.02;
D1 = linspace(0, D1end, DN);  % 예: 50개의 샘플

% 0.02 ~ 0.3 구간: 밀도가 낮은 부분
D2 = linspace(D1end, D2end, DN + 1); % 예: 50개의 샘플

% 두 구간 병합
Dinput = [D1, D2(2:end)]; % 중복되는 값 제거

D2end = 5;
Dinput = linspace(0,D2end, 100);


kx = linspace(-pi, pi, N+1);
ky = linspace(-pi, pi, N+1);
kx(end) = [];
ky(end) = [];
qx_scan = kx;

[kx, ky] = meshgrid(kx, ky);
kx = reshape(kx, [numel(kx), 1]);
ky = reshape(ky, [numel(ky), 1]);

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

Dpout_u1 = zeros(1, length(Dinput));
E0pout1 = zeros(1, length(Dinput));
Dpout_u2 = zeros(1, length(Dinput));
E0pout2 = zeros(1, length(Dinput));
Dpout_u3 = zeros(1, length(Dinput));
E0pout3 = zeros(1, length(Dinput));
DeltaT = zeros(1, length(T_ran));
%% p-, es-, d- wave

disp('start p-, es-, d-wave')

dqx = qx_scan(2)-qx_scan(1);
jjp = 40;

% figure
% hold on

for j = 1: length(T_ran)
    j
    T = T_ran(j);
    J = J_ran(jjp);

    % qx = - dqx 계산
    for i = 1: length(Dinput)
            Delta = Dinput(i);

            J = J_ran(jjp);

            [Xi, La, gpka] = tripletXiLaMaker(kx, ky, -dqx, 0, t, J, mu);
            [Ea1, Dua1, ~, ~] = tripletFFgapsolver_finiteTemp(Xi, La, gpka, Delta, V, T);

            Dpout_u1(i) = Dua1;
            E0pout1(i) = Ea1;

    end

    [~ , I] = min(abs(Dpout_u1 ./ (Dinput + epsilon) - 1));
    if (Dpout_u1(2) - 0) / (Dinput(2) - Dinput(1)) < 1
        I = 1;
        Dp1 = 0;
        E0p1 = E0pout1(I);
    else
        Dp1 = Dinput(I);
        E0p1 = E0pout1(I);
    end
    Fq(1) = E0p1;

    % qx = 0 계산
    for i = 1: length(Dinput)
            Delta = Dinput(i);

            J = J_ran(jjp);

            [Xi, La, gpka] = tripletXiLaMaker(kx, ky, 0, 0, t, J, mu);
            [Ea2, Dua2, ~, ~] = tripletFFgapsolver_finiteTemp(Xi, La, gpka, Delta, V, T);

            Dpout_u2(i) = Dua2;
            E0pout2(i) = Ea2;

    end
    % plot(Dinput,Dout1)

    [~ , I] = min(abs(Dpout_u2 ./ (Dinput + epsilon) - 1));
    if (Dpout_u2(2) - 0) / (Dinput(2) - Dinput(1)) < 1
        I = 1;
        Dp2 = 0;
        E0p2 = E0pout2(I);
    else
        Dp2 = Dinput(I);
        E0p2 = E0pout2(I);
    end
    Fq(2) = E0p2;

    % qx = dqx 계산
    for i = 1: length(Dinput)
            Delta = Dinput(i);

            J = J_ran(jjp);

            [Xi, La, gpka] = tripletXiLaMaker(kx, ky, dqx, 0, t, J, mu);
            [Ea3, Dua3, ~, ~] = tripletFFgapsolver_finiteTemp(Xi, La, gpka, Delta, V, T);

            Dpout_u3(i) = Dua3;
            E0pout3(i) = Ea3;

    end

    [~ , I] = min(abs(Dpout_u3 ./ (Dinput + epsilon) - 1));
    if (Dpout_u3(2) - 0) / (Dinput(2) - Dinput(1)) < 1
        I = 1;
        Dp3 = 0;
        E0p3 = E0pout3(I);
    else
        Dp3 = Dinput(I);
        E0p3 = E0pout3(I);
    end
    Fq(3) = E0p3;

    rho(j) = (Fq(3) - 2 * Fq(2) + Fq(1)) / dqx / dqx;
    DeltaT(j) = Dp2;
end

lam = 1 ./ sqrt(abs(rho));

Tcindex = find(DeltaT == 0);
Tcindex = Tcindex(1);
Tc = T_ran(Tcindex);
T_norm = T_ran(1:Tcindex-1)/Tc;

Tlowindex = find(T_norm < 0.3);
Tlowindex = Tlowindex(end);
T_low = T_norm(1:Tlowindex);


figure
plot(T_ran/Tc, DeltaT);xlabel('T/T_{c} [K]');ylabel('\Delta(T)')

figure
plot(T_ran(1:Tcindex-1)/Tc, rho(1:Tcindex-1)/rho(1));xlabel('T/T_{c} [K]');ylabel('K(0,T)/K(0, 0)')


figure
plot(T_ran(1:Tcindex-1)/Tc, lam(1:Tcindex-1) / lam(1));xlabel('T/T_{c} [K]');ylabel('\lambda(T)/\lambda(0)')


%% function

function [Xi,La, gpka] = tripletXiLaMaker(kx, ky, qx, qy, t, J, mu)

kx1 = kx + qx;
kx2 = -kx + qx;
ky1 = ky + qy;
ky2 = -ky + qy;

Jk1 = J * sin(kx1) .* sin(ky1);
Jk2 = J * sin(kx2) .* sin(ky2);

xik1 = - 2 * t * cos(kx1) - 2 * t * cos(ky1) - mu;
xik2 = - 2 * t * cos(kx2) - 2 * t * cos(ky2) - mu;
gpk = sin(kx1) + 1i * sin(ky1);

Xi = .5 * (xik1 + Jk1 +xik2 + Jk2);                            
La = .5 * (xik1 + Jk1 -xik2 - Jk2);
gpka = abs(gpk) .^ 2;

end

function [E_0, Dffout,Ek, F] = tripletFFgapsolver_finiteTemp(Xi, La, gpka, Delta, U, T)
Ek = sqrt(Xi .^2 + Delta^2 .* gpka);
Ek_u = Ek + La;
Ek_d = Ek - La;

Ek_merged = [Ek_u;Ek_d;-Ek_u;-Ek_d];
ind_u = Ek_u > 0;
ind_d = Ek_d > 0;
Delta_k = U * Delta * gpka /4 ./ Ek .* (tanh(Ek_u / 2 / T) + tanh(Ek_d / 2 / T));

Dffout = mean(Delta_k);

% E_0_k = Xi - Ek + Delta^2 * gpka ./ Ek /4 .* (tanh(Ek_u / 2 / T) + tanh(Ek_d / 2 /T))...
    % + Ek_u .* (~ind_u) + Ek_d .* (~ind_d);

if T == 0
    F = Ek_u.*1./(1+exp(Ek_u/T))+Ek_d.*1./(1+exp(Ek_d/T))+(-Ek_u).*1./(1+exp(-Ek_u/T))+(-Ek_d).*1./(1+exp(-Ek_d/T));
    F = mean(F, 'omitnan');
else
    c = max(-Ek_merged/T);
    logZ = c + log(sum(exp(-Ek_merged/T-c)));
    F = -T * logZ;
end
E_0_k = Ek_u.*1./(1+exp(Ek_u/T))+Ek_d.*1./(1+exp(Ek_d/T))+(-Ek_u).*1./(1+exp(-Ek_u/T))+(-Ek_d).*1./(1+exp(-Ek_d/T));
E_0 = mean(E_0_k,'omitnan');

% E_0 = mean(E_0_k, 'omitnan')

end