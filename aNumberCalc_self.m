% clear all

N = 1000;
ParamN = 50;
t = 0.5;
mu = -1;
mu_ran = ones(1,1)*(mu);
U_ran = ones(1,1) * 1.5;
U = U_ran(1);
J_ran = linspace(0, 1, ParamN);
% J_ran = ones(1,1) * 0.4;
V_ran = linspace(0.6, 1.5, ParamN);
V_ran = linspace(0, 1.5, ParamN);

Dend = 0.3;
Dinput = linspace(0, Dend, 100);


kx = linspace(-pi, pi, N+1);
ky = linspace(-pi, pi, N+1);
kx(end) = [];
ky(end) = [];
% [kmx, kmy] = meshgrid(kx, ky);
% kkx = reshape(kmx, [numel(kmx), 1]);
% kky = reshape(kmy, [numel(kmy), 1]);

qx = linspace(-pi, pi, N+1);qx(end)=[];
qy = qx;
[qx, qy] = meshgrid(qx, qy);
qx = reshape(qx, [numel(qx), 1]);
qy = reshape(qy, [numel(qy), 1]);
qx0line = qx(N * (N / 2) + 1 + N/2 : N * (N / 2 + 1));
qyline = qy(N * (N / 2) + 1  + N/2 : N * (N / 2 + 1));


epsilon = 1e-6;

gpk = sin(kx) + 1i * sin(ky);
gpka = abs(gpk) .^ 2;
gesk = cos(kx) + cos(ky);
geska = abs(gesk) .^ 2;
gdk = cos(kx) - cos(ky);
gdka = abs(gdk) .^ 2;
%% Delta-N 함수 그리기
% s-wave 케이스
beta = 5000; % zero temperature
Ns = zeros(1,length(Ds_self));
Nff = zeros(1,length(Ds_self));
Np = zeros(1,length(Ds_self));
Nd = zeros(1,length(Ds_self));
Nes = zeros(1,length(Ds_self));


for i = 1: length(Ds_self)
    Deltas = Ds_self(i);
    J=J_ran(i);
    [correl,~,~,~,~,~]=Hbdg_solver('s',N^2,t,mu,J,beta,1,kx,ky,Deltas,false,0,0);
    Ns(i) = sum(sum(correl(1,1,:,:))) + sum(sum(correl(2,2,:,:)));
end

% FF state 케이스

for i = 1: length(Ds_self)
    Deltas = Dinput(i);
    q = q_self(i);
    J=J_ran(i);
    [correl,~,~,~,~,~]=Hbdg_solver('s',N^2,t,mu,J,beta,1,kx,ky,Deltas,true,q,0);
    Nff(i) = sum(sum(correl(1,1,:,:))) + sum(sum(correl(2,2,:,:)));
end

% p+ip-wave 케이스

for i = 1: length(Ds_self)
    Deltap = Dp_self(i, 23);
    J=J_ran(i);
    [correl,~,~,~,~,~]=Hbdg_solver('pp11',N^2,t,mu,J,beta,1,kx,ky,Deltas,false,q,0);
    Np(i) = sum(sum(correl(1,1,:,:))) + sum(sum(correl(2,2,:,:)));
end

% d-wave 케이스
for i = 1: length(Ds_self)
    Deltad = Dd_self(i, 23);
    J=J_ran(i);
    [correl,~,~,~,~,~]=Hbdg_solver('d',N^2,t,mu,J,beta,1,kx,ky,Deltas,false,q,0);
    Nd(i) = sum(sum(correl(1,1,:,:))) + sum(sum(correl(2,2,:,:)));
end

% es-wave 케이스
for i = 1: length(Ds_self)
    Deltad = Des_self(i, 23);
    J=J_ran(i);
    [correl,~,~,~,~,~]=Hbdg_solver('es',N^2,t,mu,J,beta,1,kx,ky,Deltas,false,q,0);
    Nes(i) = sum(sum(correl(1,1,:,:))) + sum(sum(correl(2,2,:,:)));
end

figure
plot(Ns); title('s-wave')
figure
plot(Nff); title('FF state')
figure
plot(real(Np)); title('p+ip-wave')
figure
plot(Nd); title('d-wave')
figure
plot(Nes); title('es-wave')
%% function
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