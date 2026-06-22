clear all

N = 100;
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

Dend = 1;
Dinput = linspace(0, Dend, 100);


kx = linspace(-pi, pi, N+1);
ky = linspace(-pi, pi, N+1);
kx(end) = [];
ky(end) = [];
[kmx, kmy] = meshgrid(kx, ky);
kkx = reshape(kmx, [numel(kmx), 1]);
kky = reshape(kmy, [numel(kmy), 1]);

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
J=0.2;
for i = 1: length(Dinput)
    Deltas = Dinput(i);
    for j = 1: length(kx)
        for l = 1: length(ky)
            kkx = kx(j);
            kky = ky(j);
            ek = - 2* t * cos(kkx) - 2 * t * cos(kky);
            Jk = J * sin(kkx) * sin(kky);
            Hbdg = [(ek-mu)+Jk 0 0 Deltas;
                0 (ek-mu)-Jk -Deltas 0;
                0 -Deltas -(ek-mu)+Jk 0;
                Deltas 0 0 -(ek-mu)-Jk];
            % Hbdg = [(ek-mu)-Jk -Deltas;
            %         -Deltas -(ek-mu)+Jk];
            [u, d] = eig(Hbdg);
            F=diag(1./(exp(diag(d)*beta)+1));
            correl(:,:, j, l) = u * F * u';

        end
    end
    Ns3(i) = sum(sum(correl(3,3,:,:)));
    Ns1(i) = sum(sum(correl(1,1,:,:)));
    Ns(i) = sum(sum(correl(3,3,:,:))) + sum(sum(correl(1,1,:,:)));
    % Ns(i) = sum(sum(correl(1,1,:,:)));
end

% FF state 케이스

for i = 1: length(Dinput)
    Deltas = Dinput(i);
    for j = 1: length(kx)
        for l = 1: length(ky)
            kkx = kx(j);
            kky = ky(j);
            ek = - 2* t * cos(kkx) - 2 * t * cos(kky);
            Jk = J * sin(kkx) * sin(kky);
            Hbdg = [(ek-mu)+Jk 0 0 Deltas;
                0 (ek-mu)-Jk -Deltas 0;
                0 -Deltas -(ek-mu)+Jk 0;
                Deltas 0 0 -(ek-mu)-Jk];
            [u, d] = eig(Hbdg);
            F=diag(1./(exp(diag(d)*beta)+1));
            correl(:,:, j, l) = u * F * u';

        end
    end
    Nff(i) = sum(sum(correl(1,1,:,:)))+ sum(sum(correl(3,3,:,:)));
end

% p+ip-wave 케이스

for i = 1: length(Dinput)
    Deltap = Dinput(i);
    for j = 1: length(kx)
        for l = 1: length(ky)
            kkx = kx(j);
            kky = ky(j);
            gk = sin(kkx) + 1i * sin(kky);
            ek = - 2* t * cos(kkx) - 2 * t * cos(kky);
            Jk = J * sin(kkx) * sin(kky);
            Hbdg = [(ek-mu)+Jk 0 Deltap*gk 0;
                0 (ek-mu)-Jk 0 Deltap*gk;
                Deltap*gk' 0 -(ek-mu)+Jk 0;
                0 Deltap*gk' 0 -(ek-mu)-Jk];
            [u, d] = eig(Hbdg);
            F=diag(1./(exp(diag(d)*beta)+1));
            correl(:,:, j, l) = u * F * u';

        end
    end
    Np(i) = real(sum(sum(correl(1,1,:,:))))+ real(sum(sum(correl(3,3,:,:))));
end

% % d-wave 케이스
% for i = 1: length(Dinput)
%     Deltad = Dinput(i);
%     for j = 1: length(kx)
%         for l = 1: length(ky)
%             kkx = kx(j);
%             kky = ky(j);
%             gk = cos(kkx) - cos(kky);
%             ek = - 2* t * cos(kkx) - 2 * t * cos(kky);
%             Jk = J * sin(kkx) * sin(kky);
%             Hbdg = [(ek-mu)+Jk 0 0 Deltad*gk;
%                 0 (ek-mu)-Jk -Deltad*gk 0;
%                 0 -Deltad*gk -(ek-mu)+Jk 0;
%                 Deltad*gk 0 0 -(ek-mu)-Jk];
%             [u, d] = eig(Hbdg);
%             F=diag(1./(exp(diag(d)*beta)+1));
%             correl(:,:, j, l) = u * F * u';
% 
%         end
%     end
%     Nd(i) = sum(sum(correl(1,1,:,:)))+ sum(sum(correl(3,3,:,:)));
% end
% 
% % es-wave 케이스
% for i = 1: length(Dinput)
%     Deltad = Dinput(i);
%     for j = 1: length(kx)
%         for l = 1: length(ky)
%             kkx = kx(j);
%             kky = ky(j);
%             gk = cos(kkx) + cos(kky);
%             ek = - 2* t * cos(kkx) - 2 * t * cos(kky);
%             Jk = J * sin(kkx) * sin(kky);
%             Hbdg = [(ek-mu)+Jk 0 0 Deltad*gk;
%                 0 (ek-mu)-Jk -Deltad*gk 0;
%                 0 -Deltad*gk -(ek-mu)+Jk 0;
%                 Deltad*gk 0 0 -(ek-mu)-Jk];
%             [u, d] = eig(Hbdg);
%             F=diag(1./(exp(diag(d)*beta)+1));
%             correl(:,:, j, l) = u * F * u';
% 
%         end
%     end
%     Nes(i) = sum(sum(correl(1,1,:,:)))+ sum(sum(correl(3,3,:,:)));
% end

figure
plot(Dinput, Ns/N/N/2); title('s-wave')
figure
plot(Dinput, Nff/N/N/2 - Nff(1)/N/N/2); title('FF state')
figure
plot(Dinput, Np/N/N/2 - Np(1)/N/N/2); title('p+ip-wave')
% figure
% plot(Dinput, Nd/N/N/2 - Nd(1)/N/N/2); title('d-wave')
% figure
% plot(Dinput, Nes/N/N/2 - Nes(1)/N/N/2); title('es-wave')
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

function lambda = calculate_PenetrationDepth_Rho(t, en)
    lambda = zeros(1, length(t));
    % 자유에너지 계산
    en = en(en > 0);
    for i=1:length(t)
        lambda(i)=(1+2*sum(exp(-1/t(i)*en)./((exp(-1/t(i)*en)+1).^2)));
    end
    plot(t,lambda)
    hold on
    xlabel('Temperature, T')
    ylabel('Penetration depth, \lambda(T)')
    fontsize(16,"points")
end