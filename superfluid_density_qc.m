%% 
clear all
%% 
t=0.5;
J_ran = linspace(0, 1, 50);
mu=-1;
m=1/40;
U=-1.5;
% Vs=-0.1;
Delta = 0.018;
Delta = 0.2305;
str='J/2t=%.2f, \\mu/2t=%.2f, U/2t=%.2f';
% str='t=%.2f, J=%.2f, \\mu=%.2f, V_s=%.2f, \\Delta=%.2f';
beta=5000;
NOGx=1001; % Number Of Grid x
NOGy=1001; % Number Of Grid y
Np=NOGx-1;
NF=1/((NOGx-1)*(NOGy-1)); % Normalization Factor
kx = linspace(-pi, pi, N+1);
ky = linspace(-pi, pi, N+1);
kx(end) = [];
ky(end) = [];

qqx = kx(kx<pi/8 & kx>-0.0001);
qqy = ky(ky<pi/8 & ky>-0.0001);

[kx, ky] = meshgrid(kx, ky);
kx = reshape(kx, [numel(kx), 1]);
ky = reshape(ky, [numel(ky), 1]);
superfluid_density_qcc = zeros(1,22);
index_test = zeros(1,22);
dq = qqx(2) - qqx(1);
for jj = 1: 22
J = J_ran(jj);
    jj
    Delta = Dff_self(jj);
for i=1:length(qqx)
    % for j=1:length(qqy)
        qx = qqx(i);
        [Xi,La] = XiLaMaker(kx, ky, qx, 0, t, J, mu);
        [E_0, Dffout,Ek] = FFgapsolver(Xi, La, Delta, U);
        EE(i) = E_0;

    % end
end
    superfluid_density = gradient(gradient(EE))/dq/dq;
    index_qc = abs(qqx - q_self(jj)) < 0.0001;
    index_qc = find(index_qc == true);
    index_test(jj) = index_qc;
    % index_qc = find(EE == min(EE));
    % index_test(jj) = index_qc(1);
    superfluid_density_qcc(jj) = superfluid_density(index_qc(1));
end
% mode=floor(mod(EE,10^3));
% E=round(EE-mode+10^3-1,3);
E=fix((EE*10^(-3))*10^(11))/10^11;
%% 1D plot

figure
plot(E)
title(sprintf(str, t, J, mu))
[xx, yy]=find(E==min(E));
qqx(yy(1))
%% finding 2D minimum q points

minE=min(min(E));
[minx, miny]=find(E==minE);
qxmin=qqx(minx);
qymin=qqy(miny);
valid_indices = find(qxmin >= -pi/2 & qxmin <= pi/2 & qymin >= -pi/2 & qymin <= pi/2);

% scatter를 이용하여 [-pi/2, pi/2] 범위 내에 속하는 점만 찍기
% scatter(qqx(valid_indices), qqy(valid_indices), 'o');
%% 2D q surf
figure
surf(qqx,qqy,E)
xlabel('q_x');ylabel('q_y');
% xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
% yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
xticks([-pi/2 0 pi/2]);xticklabels({'-\pi/2', '0', '\pi/2'});
yticks([-pi/2 0 pi/2]);yticklabels({'-\pi/2', '0', '\pi/2'});
%% 2D q case contour and scatter
figure
contour(qqx,qqy,E);hold on;scatter(qxmin,qymin,'k','filled')
xlabel('q_x');ylabel('q_y');
% xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
% yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
xlim([-pi/2 pi/2]);ylim([-pi/2 pi/2])
xticks([-pi/2 0 pi/2]);xticklabels({'-\pi/2', '0', '\pi/2'});
yticks([-pi/2 0 pi/2]);yticklabels({'-\pi/2', '0', '\pi/2'});
title(sprintf(str, J, mu, Vs))
axis square
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
%% Superfluid density (second derivative of free energy)

[delEx,delEy]=gradient(E);
[delExx,delExy]=gradient(delEx);
[delEyx,delEyy]=gradient(delEy);
Superfluid_Density_tensor = zeros(Np, Np, 2, 2);
Superfluid_Density_tensor(:, :, 1, 1) = delExx;  % 텐서의 (1,1) 위치에 delExx
Superfluid_Density_tensor(:, :, 1, 2) = delExy;  % 텐서의 (1,2) 위치에 delExy
Superfluid_Density_tensor(:, :, 2, 1) = delEyx;  % 텐서의 (2,1) 위치에 delEyx
Superfluid_Density_tensor(:, :, 2, 2) = delEyy;  % 텐서의 (2,2) 위치에 delEyy

%% function
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
