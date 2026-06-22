% clear all

N = 1000;
ParamN = 50;
t = 0.5;
mu = -1;
V=0.7;
J_ran = linspace(0, 1, ParamN);


kx = linspace(-pi, pi, N+1);
ky = linspace(-pi, pi, N+1);
kx(end) = [];
ky(end) = [];
dqx = kx(2)-kx(1);

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

Fq = zeros(1, 3);

%% s-wave
disp('start s-wave')

figure
hold on

for j = 1: 26
    j
    J = J_ran(j);
    Delta = Dff_self(j);
    qx = q_self(j);
    J = J_ran(j);
    Jk = J .* sin(kx) .* sin(ky);

    [Xi, La] = XiLaMaker(kx, ky, qx-dqx, 0, t, J, mu);
    [Ea1, ~, ~] = FFgapsolver_finiteTemp(Xi, La, Delta, V, 0);
    Fq(1) = Ea1;
    [Xi, La] = XiLaMaker(kx, ky, qx, 0, t, J, mu);
    [Ea2, ~, ~] = FFgapsolver_finiteTemp(Xi, La, Delta, V, 0);
    Fq(2) = Ea2;
    [Xi, La] = XiLaMaker(kx, ky, qx+dqx, 0, t, J, mu);
    [Ea3, ~, ~] = FFgapsolver_finiteTemp(Xi, La, Delta, V, 0);
    Fq(3) = Ea3;


    plot(Fq)
    rho(j) = (Fq(3) - 2 * Fq(2) + Fq(1)) / dqx / dqx;

end

figure
plot(J_ran(1:26), rho);xlabel('J/2t');ylabel('n_s')
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

function [E_0, Dffout,Ek] = FFgapsolver_finiteTemp(Xi, La, Delta, U, T)
Ek = sqrt(Xi .^2 + Delta^2);
Ek_u = Ek + La;
Ek_d = Ek - La;

% ind_u = Ek_u > 0;
% ind_d = Ek_d > 0;
Delta_k = U * Delta /4 ./ Ek .* (tanh(Ek_u / 2 / T) + tanh(Ek_d / 2 / T));

Dffout = mean(Delta_k);

% E_0_k = Xi - Ek + Delta^2 ./ Ek /4 .* (tanh(Ek_u / 2 / T) + tanh(Ek_d / 2 /T))...
%     + Ek_u.*1./(1+exp(Ek_u/T))+Ek_d.*1./(1+exp(Ek_d/T))+(-Ek_u).*1./(1+exp(-Ek_u/T))+(-Ek_d).*1./(1+exp(-Ek_d/T));
E_0_k = Ek_u.*1./(1+exp(Ek_u/T))+Ek_d.*1./(1+exp(Ek_d/T))+(-Ek_u).*1./(1+exp(-Ek_u/T))+(-Ek_d).*1./(1+exp(-Ek_d/T));
E_0 = mean(E_0_k, 'omitnan');

end
















function [Ek_u, Ek_d] = signeltEkJkMaker(kx, ky, ek, mu, Delta, J,gka)
Jk = J * sin(kx) .* sin(ky);
Ek = sqrt((ek-mu).^2 + Delta^2*gka);
Ek_u = Ek + Jk;
Ek_d = Ek - Jk;
end

function [E0out, Dout,correl,minEk_d] = singletgapsolver(ek, mu, Delta, Jk, U)
Ek = sqrt((ek-mu).^2 + Delta^2);
Ek_u = Ek + Jk;
Ek_d = Ek - Jk;
ind_u = Ek_u > 0;
ind_d = Ek_d > 0;
Delta_k = U * Delta / 2 ./ Ek .* (ind_d - (~ind_u));
E_0_s_k = (ek - mu) - Ek + Delta^2 ./ Ek / 2 .* (ind_u - (~ind_d))...
    + Ek_u .* (~ind_u) + Ek_d .* (~ind_d);
correl = Delta ./ Ek / 2 .* (ind_u - (~ind_d));
Dout = mean(Delta_k, 'omitnan');
E0out = mean(E_0_s_k, 'omitnan');
minEk_d = min(Ek_d);
end

function [E0out, Dout,correl,minEk_d] = singletgapsolver_finiteTemp(ek, mu, Delta, Jk, U, T)
Ek = sqrt((ek-mu).^2 + Delta^2);
Ek_u = Ek + Jk;
Ek_d = Ek - Jk;
ind_u = Ek_u > 0;
ind_d = Ek_d > 0;
Delta_k = U * Delta / 4 ./ Ek .* (tanh(Ek_u/ 2 / T) + tanh(Ek_d / 2/ T));
E_0_s_k = (ek - mu) - Ek + Delta^2 ./ Ek / 4 .* (tanh(Ek_u / 2 / T) + tanh(Ek_d / 2/ T))...
    + Ek_u .* (~ind_u) + Ek_d .* (~ind_d);
correl = Delta ./ Ek / 2 .* (tanh(Ek_u/ 2 / T) + tanh(Ek_d / 2/ T));
Dout = mean(Delta_k, 'omitnan');
E0out = mean(E_0_s_k, 'omitnan');
minEk_d = min(Ek_d);
end

function [E_0, Dffout,Ek] = FFgapsolver(Xi, La, Delta, U)
Ek = sqrt(Xi .^2 + Delta^2);
Ek_u = Ek + La;
Ek_d = Ek - La;

ind_u = Ek_u > 0;
ind_d = Ek_d > 0;
Delta_k = U * Delta /2 ./ Ek .* (ind_u - (~ind_d));

Dffout = mean(Delta_k);

E_0_k = Xi - Ek + Delta^2 ./ Ek /2 .* (ind_u -(~ind_d))...
    + Ek_u .* (~ind_u) + Ek_d .* (~ind_d);

E_0 = mean(E_0_k, 'omitnan');


% E1ff=mean(Xi - Ek - Delta^2 ./ Ek /2 .* (ind_u -(~ind_d)));
% E2ff=mean(Ek_u .* (~ind_u) + Ek_d .* (~ind_d));
end