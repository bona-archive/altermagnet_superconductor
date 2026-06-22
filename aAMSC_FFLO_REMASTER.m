%% 
clear all

%% parameter
t=1;
J=0.2;
mu=-1;
m=1/40;
U=1.5;
% Vs=-0.1;
% Del=0.018;
Del = 0.2;
Del = 0.2305;
Del = 0.15; % J=0.4 mu=-1 FF state magnify
Del = 0.2305; % J=0.2 mu=-1 s-wave state
% Del = 0.025; % J=0.4 mu=-1 FF state

str='J/2t=%.2f, \\mu/2t=%.2f, U/2t=%.2f';
% str='t=%.2f, J=%.2f, \\mu=%.2f, V_s=%.2f, \\Delta=%.2f';
beta=5000;
NOGx=501; % Number Of Grid x
NOGy=501; % Number Of Grid y
Np=NOGx-1;
NF=1/((NOGx-1)*(NOGy-1)); % Normalization Factor
kkx=linspace(-pi,pi,NOGx);kkx(end)=[];
kky=linspace(-pi,pi,NOGy);kky(end)=[];
[kkx, kky] = meshgrid(kkx, kky);

qend=pi;
Npselect = round(Np/16);
Npendp = Np/2+1+Npselect;
Npendm = Np/2+1-Npselect;

% qqx=kkx(end/2+1:end);
% qqy=kky(end/2+1:end);
test=linspace(-pi,pi,NOGx);test(end)=[];
qqx=linspace(-qend/8,qend/8,NOGx);qqx(end)=[];
qqy=linspace(-qend/8,qend/8,NOGy);qqy(end)=[];

qqx = test(Np/2+1-Npselect:Np/2+1+Npselect);

qqy = test(Np/2+1-Npselect:Np/2+1+Npselect);
%% colormap

% Define the color shifting range to avoid pure white
shift_start = 0.1;  % Starting point of the colormap to avoid the lightest colors
shift_end = 0.9;    % Ending point of the colormap to avoid the lightest colors

% 'hot' 컬러맵을 불러옴
numColors = 256;  % Increase the number of colors for a smoother gradient
original_hot = hot(numColors);
original_hot = original_hot(round(shift_start * numColors):round(shift_end * numColors), :);

% 각 색상 채널을 반전시켜 보색 구하기
complementary_hot = 1 - original_hot;

% 원래 컬러맵과 보색 컬러맵을 이어붙임
custom_colormap = [complementary_hot; original_hot];
%% calculation

for i=1:length(qqx)
    for j=1:length(qqy)

        [Xi,La] = XiLaMaker(kkx, kky, qqx(i), qqy(j), t, J, mu);

        Delta = Del;

        [EEq, ~, ~, ~, ~, ~] = FFgapsolver(Xi, La, Delta, U);
        EE(i, j) = EEq;
        % Dffout(i, j) = DDq;
    end
end
% mode=floor(mod(EE,10^3));
% E=round(EE-mode+10^3-1,3);
% E=fix((EE*10^(-3))*10^(11))/10^11;
E=EE;
%% 1D plot
% for i=1:length(qqx)
% 
% 
%         [Xi,La] = XiLaMaker(kkx, kky, qqx(i), qqy(j), t, J, mu);
% 
%         Delta = Dff_self(i);
% 
%         [EEq, ~, ~, ~, ~, ~] = FFgapsolver(Xi, La, Delta, U);
%         EE(i, j) = EEq;
%         % Dffout(i, j) = DDq;
% 
% end
figure
plot(qqx, E,'k'), axis square
xlim([qqx(1), qqx(end)])
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
surf(qqx,qqy,E,'LineStyle','none');view(2);colormap(jet);
% xlabel('q_x');ylabel('q_y');
% xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
% yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
xticks([-pi/2 0 pi/2]);xticklabels({'', '', ''});
yticks([-pi/2 0 pi/2]);yticklabels({'', '', ''});
xlim([-pi/8 pi/8]);ylim([-pi/8 pi/8])
axis square
%% 2D q case contour and scatter
figure
contour(qqx,qqy,E);hold on;scatter(qxmin,qymin,'k','filled')
xlabel('q_x');ylabel('q_y');
% xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
% yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
xlim([-pi/2 pi/2]);ylim([-pi/2 pi/2])
xticks([-pi/2 0 pi/2]);xticklabels({'-\pi/2', '0', '\pi/2'});
yticks([-pi/2 0 pi/2]);yticklabels({'-\pi/2', '0', '\pi/2'});
title(sprintf(str, J, mu, U))
axis square
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
%% Superfluid density (second derivative of free energy)

[delEx,delEy]=gradient(E);
[delExx,delExy]=gradient(delEx);
[delEyx,delEyy]=gradient(delEy);
% Superfluid_Density_tensor = zeros(Np, Np, 2, 2);
Superfluid_Density_tensor(:, :, 1, 1) = delExx;  % 텐서의 (1,1) 위치에 delExx
colorlim_Exx = colorlim(delExx);
Superfluid_Density_tensor(:, :, 1, 2) = delExy;  % 텐서의 (1,2) 위치에 delExy
colorlim_Exy = colorlim(delExy);
Superfluid_Density_tensor(:, :, 2, 1) = delEyx;  % 텐서의 (2,1) 위치에 delEyx
colorlim_Eyx = colorlim(delEyx);
Superfluid_Density_tensor(:, :, 2, 2) = delEyy;  % 텐서의 (2,2) 위치에 delEyy
colorlim_Eyy = colorlim(delEyy);

xyend = pi/8;
figure;surf(qqx, qqy, squeeze(Superfluid_Density_tensor(:, :, 1, 1)),'Linestyle','none');view(2);colormap(custom_colormap);clim([-colorlim_Exx, colorlim_Exx]);
axis square
xlim([-xyend xyend]);ylim([-xyend xyend])
xticks([-xyend 0 xyend]);
yticks([-xyend 0 xyend]);
xticklabels({'', '', ''});
yticklabels({'', '', ''});
figure;surf(qqx, qqy, squeeze(Superfluid_Density_tensor(:, :, 1, 2)),'Linestyle','none');view(2);colormap(custom_colormap);clim([-colorlim_Exy, colorlim_Exy]);
axis square
xlim([-xyend xyend]);ylim([-xyend xyend])
xticks([-xyend 0 xyend]);
yticks([-xyend 0 xyend]);
xticklabels({'', '', ''});
yticklabels({'', '', ''});
figure;surf(qqx, qqy, squeeze(Superfluid_Density_tensor(:, :, 2, 1)),'Linestyle','none');view(2);colormap(custom_colormap);clim([-colorlim_Eyx, colorlim_Eyx]);
axis square
xlim([-xyend xyend]);ylim([-xyend xyend])
xticks([-xyend 0 xyend]);
yticks([-xyend 0 xyend]);
xticklabels({'', '', ''});
yticklabels({'', '', ''});
figure;surf(qqx, qqy, squeeze(Superfluid_Density_tensor(:, :, 2, 2)),'Linestyle','none');view(2);colormap(custom_colormap);clim([-colorlim_Eyy, colorlim_Eyy]);
axis square
xlim([-xyend xyend]);ylim([-xyend xyend])
xticks([-xyend 0 xyend]);
yticks([-xyend 0 xyend]);
xticklabels({'', '', ''});
yticklabels({'', '', ''});
%% 1D q case plot and scatter
% figure
% plot(qqx,E)
% xlabel('q_x');ylabel('q_y');
% hold on
% scatter(qqx(minx),ones(1,4)*minE)
%% Degenerated four points overlap
% 
% phiA=linspace(-pi,pi,10);
% phiB=phiA;phiC=phiA;phiD=phiA;
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


function [E_0, Dffout, correl, Ek_u, ind_u, ind_d] = FFgapsolver(Xi, La, Delta, U)
Ek = sqrt(Xi .^2 + Delta^2);
Ek_u = Ek + La;
Ek_d = Ek - La;

ind_u = Ek_u > 0;
ind_d = Ek_d > 0;
Delta_k = U * Delta /2 ./ Ek .* (ind_u - (~ind_d));

Dffout = mean(Delta_k, 'omitnan');

E_0_k = Xi - Ek + Delta^2 ./ Ek /2 .* (ind_u -(~ind_d))...
    + Ek_u .* (~ind_u) + Ek_d .* (~ind_d);
correl = Delta ./ Ek /2 .* (ind_u -(~ind_d));
E_0 = mean(mean(E_0_k, 'omitnan'), 'omitnan');

end

function draw_gapstructure(gap_1, gap_2, N)
figure
hold on
plot([gap_1(N/2+1:N,N/2+1);gap_1(1,N/2+1:N)';flip(diag(gap_1(N/2+1:N,N/2+1:N)));diag(gap_1(N/2+1:-1:2,N/2+1:N));gap_1(1:N/2+1,1);gap_1(N/2+1,N:-1:N/2+1)'],'k','LineWidth',1)
plot([gap_2(N/2+1:N,N/2+1);gap_2(1,N/2+1:N)';flip(diag(gap_2(N/2+1:N,N/2+1:N)));diag(gap_2(N/2+1:-1:2,N/2+1:N));gap_2(1:N/2+1,1);gap_2(N/2+1,N:-1:N/2+1)'],'k','LineWidth',1)
plot([-gap_1(N/2+1:N,N/2+1);-gap_1(1,N/2+1:N)';flip(diag(-gap_1(N/2+1:N,N/2+1:N)));diag(-gap_1(N/2+1:-1:2,N/2+1:N));-gap_1(1:N/2+1,1);-gap_1(N/2+1,N:-1:51)'],'k','LineWidth',1)
plot([-gap_2(N/2+1:N,N/2+1);-gap_2(1,N/2+1:N)';flip(diag(-gap_2(N/2+1:N,N/2+1:N)));diag(-gap_2(N/2+1:-1:2,N/2+1:N));-gap_2(1:N/2+1,1);-gap_2(N/2+1,N:-1:N/2+1)'],'k','LineWidth',1)
plot(zeros(1,3*N+1),'k')
xticks([0 N/2 N 3*N/2 2*N 5*N/2 3*N])
xticklabels({'\Gamma', 'X', 'M','\Gamma','M', 'Y', '\Gamma'})
xlim([0 3*N])
ax1=gca;
ax1.XGrid='on';
ax1.GridAlpha=1;
% ax1.FontName='Times New Roman';
ax1.FontSize=12;
box on
end

function draw_toponodalstructure(mu,ek,Jk,D_self,N,gk,gka)
figure
hold on
surf(reshape(ek-mu,[N,N]),'FaceColor','r','EdgeColor','none')
surf(reshape(ek-mu,[N,N]),'FaceColor','b','EdgeColor','none')
surf(reshape((gk).^2,[N,N]),'FaceColor','k','EdgeColor','none')
surf(reshape(sqrt((ek - mu) .^2 + D_self ^2 * gka)+Jk,[N,N]),'FaceColor','m','EdgeColor','none')
surf(reshape(sqrt((ek - mu) .^2 + D_self ^2 * gka)-Jk,[N,N]),'FaceColor','c','EdgeColor','none')
% surf(reshape(Jk,[N,N]))
zlim([-0.01,0.01])
box on
end


function colorlimit = colorlim(data)
% 데이터의 최댓값과 최솟값 구하기
data_min = min(data(:));
data_max = max(data(:));

% 최솟값과 최댓값의 절대값 중 더 큰 값을 기준으로 대칭 설정
colorlimit = max(abs(data_min), abs(data_max));
end