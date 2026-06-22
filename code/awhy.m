%% plot ground state energy


% figure  % E0s - E0p
% hold on
% xlabel('J')
% ylabel('V')
% title( 'E_{0s}-E_{0p}')
% surf(J_ran, V_ran, E0s_self - E0p_self)
% colorbar
% [Csp, hsp]=contour(J_ran, V_ran, E0p_self - E0s_self, [0, 0], 'k', Linewidth=1.5);
% axis equal
% axis square
% 
% figure  % E0s - E0ff
% hold on
% xlabel('J')
% ylabel('V')
% title( 'E_{0s}-E_{0eff}')
% surf(J_ran, V_ran, E0s_self - E0ff_self)
% colorbar
% [Csff, hsff]=contour(J_ran, V_ran, E0s_self - E0ff_self, [1e-8, 1e-8], 'k', Linewidth=1.5);
% axis equal
% axis square


% custom colormap
customCMap = [
    1, 0.83, 0.83;  % Color near 1: #F0D4D4
    1, 0.90, 0.80;  % Color near 2: #FFE6CC
    1, 1, 1;        % Color near 3: White
    0.90, 0.86, 0.86; % Color near 4: #E6DBDC
    0.90, 0.94, 1;  % Color near 5: #E6EFFF
];



% figure
% ax = gca;
% hold on
% xlabel('J')
% ylabel('V')
% title( 'phase boundary')
% % [Cbps, hbps]=contour(J_ran, V_ran, E0p_self - E0s_self, [0, 0], 'k', Linewidth=1.5);
% [Cbsp, hbsp]=contour(J_ran, V_ran, E0p_self - E0ff_self, [0, 0], 'k', Linewidth=1.5);
% [Cbsff, hbsff]=contour(J_ran, V_ran, E0s_self - E0ff_self, [1e-8, 1e-8], 'k', Linewidth=1.5);
% [Cbfs, hbfs]=contour(J_ran, V_ran, BFSfinder, [0, 0], 'r', Linewidth=1.5);
% [Cbss, hbsd]=contour(J_ran, V_ran, E0s_self - E0d_self, [1e-8, 1e-8], 'k', Linewidth=1.5);
% [Cbpd, hbpd]=contour(J_ran, V_ran, E0p_self - E0d_self, [1e-8, 1e-8], 'k', Linewidth=1.5);
% [Cbdff, hbdff]=contour(J_ran, V_ran, E0ff_self - E0d_self, [1e-8, 1e-8], 'k', Linewidth=1.5);
% [Cbesff, hbesff]=contour(J_ran,V_ran V_ran, E0ff_self - E0es_self, [1e-8, 1e-8], 'k', Linewidth=1.5);
% [Cbpes, hbpes]=contour(J_ran, V_ran, E0p_self - E0es_self, [1e-8, 1e-8], 'k', Linewidth=1.5);
% [Cbses, hbses]=contour(J_ran, V_ran, E0s_self - E0es_self, [1e-8, 1e-8], 'k', Linewidth=1.5);
% 
% box on
% ax.LineWidth = 1;

% figure
% q_space = linspace(qyline(1), pi, length(V_ran));
% ax2 = gca;
% hold on
% plot(J_ran,q_self,'r','LineStyle','none','LineWidth',1,'Marker','v','MarkerSize',12);
% [Cbsff2, hbsff2]=contour(J_ran, q_space, E0s_self - E0ff_self, [1e-8, 1e-8], 'k', Linewidth=1.5);
% xlabel('J')
% ylabel('q')
% yticks([0, pi/2, pi])
% yticklabels({'0', '\pi/2', '\pi'})
% ylim([0,pi])
% box on
% ax2.LineWidth = 1;
% ax2.FontSize = 12;

% q_self sorting process..

for i = 1:length(J_ran)
    if pi/2 < q_self(i) && q_self(i) < pi
        q_self(i) = -1 * (q_self(i)-pi);
    else
    end

end



% phase sorting process..
allMatrices = cat(3, E0s_self, E0ff_self, E0es_self, E0d_self, E0p_self);

[~, minIndex] = min(allMatrices, [], 3);
for vvv = 1: length(V_ran)
    for jjj = 1: length(J_ran)
        aaa = minIndex(vvv,jjj);
        q_check = q_self(jjj);
        if aaa == 2 && q_check == 0
            minIndex(vvv,jjj) = 1;
        else
        end
    end
end
figure
surf(J_ran,V_ran,minIndex,'EdgeColor','none')
xlabel('J')
ylabel('V')
view(2)
xlim([J_ran(1),J_ran(end)])
ylim([V_ran(1),V_ran(end)])
colormap(customCMap);
clim([1,5])
box on

axis normal
box on

% hold on
% [cc1,hh1]=contour(J_ran,V_ran,minIndex,[1,1],'k',Linewidth=1.5);hh1.ContourZLevel=15;
% [cc,hh]=contour(J_ran,V_ran,minIndex,[5,5],'k',Linewidth=1.5);hh.ContourZLevel=15;


% cut J and plot q-J
[mm,ii]=find(minIndex==1 | minIndex == 2);
J_cut = max(J_ran(ii));

%% pairing amplitude

figure
hold on;
plot(J_ran,Ds_self,'o','MarkerSize',8,'LineWidth',2)
plot(J_ran,Dff_self,'+','MarkerSize',8,'LineWidth',2)
plot(J_ran,Dp_self(23,:),'o','MarkerSize',8,'LineWidth',2)
plot(J_ran,Dd_self(23,:),'o','MarkerSize',8,'LineWidth',2)
plot(J_ran,Des_self(23,:),'o','MarkerSize',8,'LineWidth',2)
xlabel('J/2t')
ylabel('Pairing amplitude')

legend({'\Delta^{s}_{\uparrow\downarrow}','\Delta^{FF}_{\uparrow\downarrow}','\Delta^{p}_{\uparrow\uparrow}','\Delta^{d}_{\uparrow\downarrow}','\Delta^{es}_{\uparrow\downarrow}'})
axis square
%% 

figure
yyaxis left
xlim([0.28, J_cut])
yyaxis left
ylim([-0.01,pi/4])
yyaxis left
xlabel('J/2t')
yyaxis left
ylabel('q_{c}')
yyaxis left
yticks([0 pi/4])
yyaxis left
yticklabels({'0', '\pi/4'})
plot(J_ran,q_self,'LineStyle','-')
hold on
plot(J_ran, zeros(1,length(q_self)),'LineStyle','--')

yyaxis right
xlim([0.28, J_cut])
yyaxis right
ylim([-0.15, 0.15])
yyaxis right
ylabel('BFS area')
plot(J_ran,BFS_areas1(9, :), 'LineStyle','--')
hold on
plot(J_ran,BFS_areaff1(9, :), 'LineStyle','-')

%% topological nodal structure

jj = 21;
Jk = J_ran(jj) * sin(kx) .* sin(ky);
kkx=linspace(-pi,pi,N+1);kkx(end)=[];


draw_singlettoponodalstructure(mu, ek, Jk, kkx, N, ones(1,N^2))
title('s-wave nodal structure')

draw_singlettoponodalstructure(mu, ek, Jk, kkx, N, gesk)
title('es-wave nodal structure')

draw_singlettoponodalstructure(mu, ek, Jk, kkx, N, gdk)
title('d-wave nodal structure')

draw_triplettoponodalstructure(mu, ek, Jk, kkx, N, gpk)

jff=21;
[xik1,Jk1,xik2,Jk2,Xi,~] = XiLaMaker(kx, ky, 0, q_self(jff), t, J_ran(jff), mu);
draw_FFtoponodalstructure(Xi,xik1+Jk1,xik2-Jk2,kkx,N)
title('FFLO nodal structure')
%% gap structure
vv = 23;% V_ran(23)=1.00 in V:0.6~1.5 version
jjs = 11;% J_ran(11)=0.2, J_ran(16)=0.3, J_ran(21)=0.41
jjff = 21;
jjp = 40;
jj=11;
Jk = J_ran(jj) * sin(kx) .* sin(ky);
mu=-1;
kx1 = kx;
ky1 = ky + q_self(jjff);
kx2 = -kx;
ky2 = -ky + q_self(jjff);

Jk1 = J_ran(jjff) * sin(kx1) .* sin(ky1);
Jk2 = J_ran(jjff) * sin(kx2) .* sin(ky2);
xik1 = - 2 * t * cos(kx1) - 2 * t * cos(ky1) - mu;
xik2 = - 2 * t * cos(kx2) - 2 * t * cos(ky2) - mu;
Xi = .5 * (xik1 + Jk1 +xik2 - Jk2);
La = .5 * (xik1 + Jk1 -xik2 + Jk2);


gap_s1 = reshape(sqrt((ek - mu) .^2 + Ds_self(jjs) ^2) + J_ran(jjs) * sin(kx) .* sin(ky), [N, N]);
gap_s2 = reshape(sqrt((ek - mu) .^2 + Ds_self(jjs) ^2) - J_ran(jjs) * sin(kx) .* sin(ky), [N, N]);
gap_ff1 = reshape(sqrt(Xi .^2 + Dff_self(jjff)^2) + La, [N, N]);
gap_ff2 = reshape(sqrt(Xi .^2 + Dff_self(jjff)^2) - La, [N, N]);

gap_d1 = reshape(sqrt((ek - mu) .^2 + Dd_self(vv, jj) ^2 * gdka) + Jk, [N, N]);
gap_d2 = reshape(sqrt((ek - mu) .^2 + Dd_self(vv, jj) ^2 * gdka) - Jk, [N, N]);
gap_es1 = reshape(sqrt((ek - mu) .^2 + Des_self(vv, jj) ^2 * geska) + Jk, [N, N]);
gap_es2 = reshape(sqrt((ek - mu) .^2 + Des_self(vv, jj) ^2 * geska) - Jk, [N, N]);

gap_p1 = reshape(sqrt((ek - mu + J_ran(jjp) * sin(kx) .* sin(ky)) .^ 2 + Dp_self(vv, jjp) ^2 * gpka), [N, N]);
gap_p2 = reshape(sqrt((ek - mu - J_ran(jjp) * sin(kx) .* sin(ky)) .^ 2 + Dp_self(vv, jjp) ^2 * gpka), [N, N]);

[s1, s2] = BFS_archiver(gap_s1, gap_s2);
[d1, d2] = BFS_archiver(gap_d1, gap_d2);
[es1, es2] = BFS_archiver(gap_es1, gap_es2);
[p1, p2] = BFS_archiver(gap_p1, gap_p2);

% draw_gapstructure(gap_s1, gap_s2, N)
% title('s-wave gap structure')
% ylim([-4,4])
% draw_gapstructure(gap_ff1, gap_ff2, N)
% title('FF-phase gap structure')
% ylim([-4,4])
% draw_gapstructure(gap_d1, gap_d2, N)
% title('d-wave gap structure')
% ylim([-4,4])
% draw_gapstructure(gap_es1, gap_es2, N)
% title('es-wave gap structure')
% ylim([-4,4])
% draw_gapstructure(gap_p1, gap_p2, N)
% title('p-wave gap structure')
% ylim([-4,4])

kkx = linspace(-pi, pi, N+1);
kky = linspace(-pi, pi, N+1);
kkx(end) = [];
kky(end) = [];

figure;ax=gca;
plot(kkx, gap_s1, 'k')
hold on;ax.LineWidth=1;fontsize(21,'points');
plot(kkx, -gap_s1, 'k')
xlim([-pi,kkx(end)]);ylim([-4,4])
xticks([-pi,-pi/2,0,pi/2,kkx(end)]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});xlabel('k_{x}');ylabel('E^{s}_{k\uparrow}')

figure;ax=gca;
plot(kkx, flip(gap_ff1), 'k')
hold on;ax.LineWidth=1;fontsize(21,'points');
plot(kkx, -flip(gap_ff1), 'k')
xlim([-pi,kkx(end)]);ylim([-4,4])
xticks([-pi,-pi/2,0,pi/2,kkx(end)]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});xlabel('k_{x}');ylabel('E^{FF}_{k\uparrow}')

figure;ax=gca;
plot(kkx, gap_p1, 'k')
hold on;ax.LineWidth=1;fontsize(21,'points');
plot(kkx, -gap_p1, 'k')
xlim([-pi,kkx(end)]);ylim([-4,4])
xticks([-pi,-pi/2,0,pi/2,kkx(end)]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});xlabel('k_{x}');ylabel('E^{p+ip}_{k\uparrow}')

% Define the color shifting range to avoid pure white
shift_start = 0.1;  % Starting point of the colormap to avoid the lightest colors
shift_end = 0.9;    % Ending point of the colormap to avoid the lightest colors

% Generate the 'hot' colormap with a higher resolution for smoother gradients
numColors = 256;  % Increase the number of colors for a smoother gradient
original_hot = hot(numColors);
original_hot = original_hot(round(shift_start * numColors):round(shift_end * numColors), :);

% Create the complementary colormap by inverting each color channel
complementary_hot = 1 - original_hot;

% Combine the complementary colormap with the original colormap
custom_colormap = [complementary_hot; original_hot];

% Display the custom colormap to visually inspect it
figure;
colormap(custom_colormap);
colorbar;


figure;ax=gca;
surf(kkx,kkx,gap_s1,'EdgeColor','none')
xlim([-pi,kkx(end)]);ylim([-pi,kkx(end)]);zlim([-0.5,0.5]);
colormap(custom_colormap);clim([-0.5,0.5]);colorbar
grid off;view(2);fontsize(21,'points');box on;ax.LineWidth=1;
xticks([-pi,0,kkx(end)]);xticklabels({'-\pi','0','\pi'});xlabel('k_{x}')
yticks([-pi,0,kkx(end)]);yticklabels({'-\pi','0','\pi'});ylabel('k_{y}')

figure;ax=gca;
surf(kkx,kkx,gap_s2,'EdgeColor','none')
xlim([-pi,kkx(end)]);ylim([-pi,kkx(end)]);zlim([-0.5,0.5]);
colormap(custom_colormap);clim([-0.5,0.5]);colorbar
grid off;view(2);fontsize(21,'points');box on;ax.LineWidth=1;
xticks([-pi,0,kkx(end)]);xticklabels({'-\pi','0','\pi'});xlabel('k_{x}')
yticks([-pi,0,kkx(end)]);yticklabels({'-\pi','0','\pi'});ylabel('k_{y}')

figure;ax=gca;
surf(kkx,kkx,gap_ff1,'EdgeColor','none')
xlim([-pi,kkx(end)]);ylim([-pi,kkx(end)]);zlim([-0.5,0.5]);
colormap(custom_colormap);clim([-0.5,0.5]);colorbar
grid off;view(2);fontsize(21,'points');box on;ax.LineWidth=1;
xticks([-pi,0,kkx(end)]);xticklabels({'-\pi','0','\pi'});xlabel('k_{x}')
yticks([-pi,0,kkx(end)]);yticklabels({'-\pi','0','\pi'});ylabel('k_{y}')

figure;ax=gca;
surf(kkx,kkx,gap_ff2,'EdgeColor','none')
xlim([-pi,kkx(end)]);ylim([-pi,kkx(end)]);zlim([-0.5,0.5]);
colormap(custom_colormap);clim([-0.5,0.5]);colorbar
grid off;view(2);fontsize(21,'points');box on;ax.LineWidth=1;
xticks([-pi,0,kkx(end)]);xticklabels({'-\pi','0','\pi'});xlabel('k_{x}')
yticks([-pi,0,kkx(end)]);yticklabels({'-\pi','0','\pi'});ylabel('k_{y}')

figure;ax=gca;
surf(kkx,kkx,gap_p1,'EdgeColor','none')
xlim([-pi,kkx(end)]);ylim([-pi,kkx(end)]);zlim([-0.5,0.5]);
colormap(custom_colormap);clim([-0.5,0.5]);colorbar
grid off;view(2);fontsize(21,'points');box on;ax.LineWidth=1;
xticks([-pi,0,kkx(end)]);xticklabels({'-\pi','0','\pi'});xlabel('k_{x}')
yticks([-pi,0,kkx(end)]);yticklabels({'-\pi','0','\pi'});ylabel('k_{y}')

figure;ax=gca;
surf(kkx,kkx,gap_p2,'EdgeColor','none')
xlim([-pi,kkx(end)]);ylim([-pi,kkx(end)]);zlim([-0.5,0.5]);
colormap(custom_colormap);clim([-0.5,0.5]);colorbar
grid off;view(2);fontsize(21,'points');box on;ax.LineWidth=1;
xticks([-pi,0,kkx(end)]);xticklabels({'-\pi','0','\pi'});xlabel('k_{x}')
yticks([-pi,0,kkx(end)]);yticklabels({'-\pi','0','\pi'});ylabel('k_{y}')
%% test
jj=16;
Jk = J_ran(jj) .* sin(kx) .* sin(ky);
Ektest = sqrt((ek - mu) .^2 + Ds_self(jj));
Ektest_u = sqrt((ek - mu) .^2 + Ds_self(jj)) + Jk;
Ektest_d = sqrt((ek - mu) .^2 + Ds_self(jj)) - Jk;
ind_u = Ektest_u > 0;
ind_d = Ektest_d > 0;

c1 = Ds_self(1) ^2 ./ Ektest / 2 .* (ind_u - (~ind_d));
c2 = Ektest_u .* (~ind_u) + Ektest_d .* (~ind_d);
figure
surf(reshape(c1,[N,N]))
hold on
figure
surf(reshape(c2,[N,N]))
view(90,0)

for jj = 1: length(J_ran)
    
    Jk = J_ran(jj) .* sin(kx) .* sin(ky);
    Ektest = sqrt((ek - mu) .^2 + Ds_self(jj) ^2);
    
    nonc(jj) = mean((ek - mu) .^2 - Ektest);
    
    Ektest_u = sqrt((ek - mu) .^2 + Ds_self(jj) ^2) + Jk;
    Ektest_d = sqrt((ek - mu) .^2 + Ds_self(jj) ^2) - Jk;
    ind_u = Ektest_u > 0;
    ind_d = Ektest_d > 0;

    c1 = Ds_self(1) ^2 ./ Ektest / 2 .* (ind_u - (~ind_d));
    c2 = Ektest_u .* (~ind_u) + Ektest_d .* (~ind_d);
    csum(jj) = mean(c1+c2);
    c11(jj) = mean(c1);
    c22(jj) = mean(c2);
end


%% pairing amplitude and ground state energy
% change box linewidth to 1
% change axis fondsize 12

%V(23)=1.00 for V_ran 0.6~1.5
%V(35)=1.00 for V_ran 0~1.5
vv=25;

figure
yyaxis left
fontsize(12,"points")
hold on
plot(J_ran,Ds_self)
plot(J_ran, Dff_self)
% plot(J_ran, Des_self(vv,:), '^-', 'MarkerSize', 7, 'Color','#F0D4D4','LineWidth',1)
% plot(J_ran, Dd_self(vv,:), 'square-', 'MarkerSize', 7, 'Color','#b8adae','LineWidth',1)
plot(J_ran, Dp_self(vv,:))
ylabel('Pairing amplitude')
xlabel('J/2t')
yyaxis right
plot(J_ran,(E0s_self(vv,:)-E0ff_self(vv,:))*N*N)
hold on
plot(J_ran,(E0s_self(vv,:)-E0p_self(vv,:))*N*N)
ylabel('Free energy')
% legend({'\Delta^{\rm{s}}_{\uparrow\downarrow}', '\Delta^{\rm{FF}}_{\uparrow\downarrow}', '\Delta^{\rm{es}}_{\uparrow\downarrow}', '\Delta^{\rm{d}}_{\uparrow\downarrow}', '\Delta^{\rm{p}}_{\uparrow\uparrow}'})
legend({'\Delta^{\rm{s}}_{\uparrow\downarrow}', '\Delta^{\rm{FF}}_{\uparrow\downarrow}', '\Delta^{\rm{p}}_{\uparrow\uparrow}'})
%\it{F^{ s}_{0} - F^{ FF}_{0}}
%\it{F^{ s}_{0} - F^{ p+ip}_{0}}
box on
axis square
%% Correlation function
V_ran = linspace(0.6, 1.5, ParamN);

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



% J_ran(11) = 0.20 
% J_ran(17) = 0.33 
% J_ran(21) = 0.41
% % % % % % % % % % % % J_ran(30) = 0.60
% J_ran(40) = 0.80
vv=23;
jj=40;
NN=1000;
kkx=linspace(-pi,pi,NN+1);kkx(end)=[];
kx = linspace(-pi, pi, NN+1);
ky = linspace(-pi, pi, NN+1);
kx(end) = [];
ky(end) = [];
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

[xik1,Jk1,xik2,Jk2,Xi,La] = XiLaMaker(kx, ky, 0, q_self(jj), t, J_ran(jj), mu);

[~, ~, correls_self, ~] = singletgapsolver(ek, mu, Ds_self(jj), J_ran(jj) * sin(kx) .* sin(ky), V_ran(vv), 1, 1);
[~, ~, correles_self, ~] = singletgapsolver(ek, mu, Des_self(vv, jj), J_ran(jj) * sin(kx) .* sin(ky), V_ran(vv), gesk, geska);
[~, ~, correld_self, ~] = singletgapsolver(ek, mu, Dd_self(vv, jj), J_ran(jj) * sin(kx) .* sin(ky), V_ran(vv), gdk, gdka);
[~, ~, ~, correlff_self] = FFgapsolver(Xi, La, Dff_self(jj), U);

correls_self=reshape(correls_self,[NN,NN]);
correles_self=reshape(correles_self,[NN,NN]);
correld_self=reshape(correld_self,[NN,NN]);
correlff_self=reshape(correlff_self,[NN,NN]);

figure % s-wave
surf(kkx,kkx,squeeze(correls_self),'EdgeColor','none');view(2)
% title('\it{s}-wave pair correlation')
xticks([-pi,-pi/2,0,pi/2,kkx(end)])
yticks([-pi,-pi/2,0,pi/2,kkx(end)])
xlim([-pi,kkx(end)])
ylim([-pi,kkx(end)])

% xticklabels({'', '', '', '', ''})
% yticklabels({'', '', '', '', ''})

xticklabels({'-\pi', '', '0', '', '\pi'})
yticklabels({'-\pi', '', '0', '', '\pi'})
fontsize(16,'points')

colormap(custom_colormap);
clim([-0.5,0.5])
axis square

figure % ff-phase
surf(kkx,kkx,squeeze(correlff_self),'EdgeColor','none');view(2)
% title('FFLO pair correlation')
xticks([-pi,-pi/2,0,pi/2,kkx(end)])
yticks([-pi,-pi/2,0,pi/2,kkx(end)])
xlim([-pi,kkx(end)])
ylim([-pi,kkx(end)])
% xticklabels({'', '', '', '', ''})
% yticklabels({'', '', '', '', ''})
xticklabels({'-\pi', '', '0', '', '\pi'})
yticklabels({'-\pi', '', '0', '', '\pi'})
fontsize(16,'points')
colormap(custom_colormap);
clim([-0.5,0.5])
axis square

figure; % es-wave
surf(kkx,kkx,correles_self,'EdgeColor','none');view(2)
% title('\it{es}-wave pair correlation')
xticks([-pi,-pi/2,0,pi/2,kkx(end)])
yticks([-pi,-pi/2,0,pi/2,kkx(end)])
xlim([-pi,kkx(end)])
ylim([-pi,kkx(end)])
% xticklabels({'', '', '', '', ''})
% yticklabels({'', '', '', '', ''})
xticklabels({'-\pi', '', '0', '', '\pi'})
yticklabels({'-\pi', '', '0', '', '\pi'})
fontsize(16,'points')
colormap(custom_colormap);
clim([-0.5,0.5])
axis square

figure; % d-wave
surf(kkx,kkx,correld_self,'EdgeColor','none');view(2)
% title('\it{d}-wave pair correlation')
xticks([-pi,-pi/2,0,pi/2,kkx(end)])
yticks([-pi,-pi/2,0,pi/2,kkx(end)])
xlim([-pi,kkx(end)])
ylim([-pi,kkx(end)])
% xticklabels({'', '', '', '', ''})
% yticklabels({'', '', '', '', ''})
xticklabels({'-\pi', '', '0', '', '\pi'})
yticklabels({'-\pi', '', '0', '', '\pi'})
fontsize(16,'points')
colormap(custom_colormap);
clim([-0.5,0.5])
axis square

%% p-wave correlation real, imag
V_ran = linspace(0.6, 1.5, ParamN);


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



% J_ran(11) = 0.20 
% J_ran(17) = 0.33 
% J_ran(21) = 0.41
% J_ran(30) = 0.60
% J_ran(40) = 0.80
% vv=23;
% jj=17;
NN=1000;
kkx=linspace(-pi,pi,NN+1);kkx(end)=[];
kx = linspace(-pi, pi, NN+1);
ky = linspace(-pi, pi, NN+1);
kx(end) = [];
ky(end) = [];
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

[~, ~, ~, correlp_self, ~] = tripletgapsolver(ek, mu, Dp_self(vv, jj), J_ran(jj) .* sin(kx) .* sin(ky), V_ran(vv), gpk, gpka);

figure % p-wave
surf(kkx,kkx,reshape(squeeze(abs(correlp_self)),[NN,NN]),'EdgeColor','none');view(2)
% title('\it{p}-wave pair correlation')
xticks([-pi,-pi/2,0,pi/2,kkx(end)])
yticks([-pi,-pi/2,0,pi/2,kkx(end)])
xlim([-pi,kkx(end)])
ylim([-pi,kkx(end)])
% xticklabels({'', '', '', '', ''})
% yticklabels({'', '', '', '', ''})
xticklabels({'-\pi', '', '0', '', '\pi'})
yticklabels({'-\pi', '', '0', '', '\pi'})
fontsize(16,'points')
colormap(custom_colormap);
clim([-0.5,0.5])
axis square

figure % p-wave
surf(kkx,kkx,reshape(squeeze(real(correlp_self)),[NN, NN]),'EdgeColor','none');view(2)
% title('\it{p}-wave pair correlation')
xticks([-pi,-pi/2,0,pi/2,kkx(end)])
yticks([-pi,-pi/2,0,pi/2,kkx(end)])
xlim([-pi,kkx(end)])
ylim([-pi,kkx(end)])
% xticklabels({'', '', '', '', ''})
% yticklabels({'', '', '', '', ''})
xticklabels({'-\pi', '', '0', '', '\pi'})
yticklabels({'-\pi', '', '0', '', '\pi'})
fontsize(16,'points')
colormap(custom_colormap);
clim([-0.5,0.5])
axis square

figure % p-wave
surf(kkx,kkx,reshape(squeeze(imag(correlp_self)),[NN, NN]),'EdgeColor','none');view(2)
% title('\it{p}-wave pair correlation')
xticks([-pi,-pi/2,0,pi/2,kkx(end)])
yticks([-pi,-pi/2,0,pi/2,kkx(end)])
xlim([-pi,kkx(end)])
ylim([-pi,kkx(end)])
% xticklabels({'', '', '', '', ''})
% yticklabels({'', '', '', '', ''})
xticklabels({'-\pi', '', '0', '', '\pi'})
yticklabels({'-\pi', '', '0', '', '\pi'})
fontsize(16,'points')
colormap(custom_colormap);
clim([-0.5,0.5])
axis square

%% More plots

for vv = 1: length(V_ran)


    for jj = 1: length(J_ran)


        Jk = J_ran(jj) * sin(kx) .* sin(ky);
        kx1 = kx;
        kx2 = -kx;
        ky1 = ky + q_self(jj);
        ky2 = -ky + q_self(jj);
        Jk1 = J_ran(jj) * sin(kx1) .* sin(ky1);
        Jk2 = J_ran(jj) * sin(kx2) .* sin(ky2);
        xik1 = - 2 * t * cos(kx1) - 2 * t * cos(ky1) - mu;
        xik2 = - 2 * t * cos(kx2) - 2 * t * cos(ky2) - mu;
        Xi = .5 * (xik1 + Jk1 +xik2 - Jk2);
        La = .5 * (xik1 + Jk1 -xik2 + Jk2);

        gap_s1 = reshape(sqrt((ek - mu) .^2 + Ds_self(jj) ^2) + Jk, [N, N]);
        gap_s2 = reshape(sqrt((ek - mu) .^2 + Ds_self(jj) ^2) - Jk, [N, N]);
        gap_ff1 = reshape(sqrt(Xi .^2 + Dff_self(jj)^2) + La, [N, N]);
        gap_ff2 = reshape(sqrt(Xi .^2 + Dff_self(jj)^2) - La, [N, N]);

        gap_d1 = reshape(sqrt((ek - mu) .^2 + Dd_self(vv, jj) ^2 * gdka) + Jk, [N, N]);
        gap_d2 = reshape(sqrt((ek - mu) .^2 + Dd_self(vv, jj) ^2 * gdka) - Jk, [N, N]);
        gap_es1 = reshape(sqrt((ek - mu) .^2 + Des_self(vv, jj) ^2 * geska) + Jk, [N, N]);
        gap_es2 = reshape(sqrt((ek - mu) .^2 + Des_self(vv, jj) ^2 * geska) - Jk, [N, N]);

        gap_p1 = reshape(sqrt((ek - mu + Jk) .^ 2 + Dp_self(vv, jj) ^2 * gpka), [N, N]);
        gap_p2 = reshape(sqrt((ek - mu - Jk) .^ 2 + Dp_self(vv, jj) ^2 * gpka), [N, N]);

        [s1, s2] = BFS_archiver(gap_s1, gap_s2);
        BFS_areas1(vv, jj) = sum(sum(s1,"omitnan"),"omitnan") / (N ^2);
        BFS_areas2(vv, jj) = sum(sum(s2,"omitnan"),"omitnan") / (N ^2);
        [ff1, ff2] = BFS_archiver(gap_ff1, gap_ff2);
        BFS_areaff1(vv, jj) = sum(sum(ff1,"omitnan"),"omitnan") / (N ^2);
        BFS_areaff2(vv, jj) = sum(sum(ff2,"omitnan"),"omitnan") / (N ^2);
        
        [d1, d2] = BFS_archiver(gap_d1, gap_d2);
        BFS_aread1(vv, jj) = sum(sum(d1,"omitnan"),"omitnan") / (N ^2);
        BFS_aread2(vv, jj) = sum(sum(d2,"omitnan"),"omitnan") / (N ^2);
        [es1, es2] = BFS_archiver(gap_es1, gap_es2);
        BFS_areaes1(vv, jj) = sum(sum(es1,"omitnan"),"omitnan") / (N ^2);
        BFS_areaes2(vv, jj) = sum(sum(es2,"omitnan"),"omitnan") / (N ^2);
        [p1, p2] = BFS_archiver(gap_p1, gap_p2);
        BFS_areap1(vv, jj) = sum(sum(p1,"omitnan"),"omitnan") / (N ^2);
        BFS_areap2(vv, jj) = sum(sum(p2,"omitnan"),"omitnan") / (N ^2);

    end
end

%% function
function [xik1,Jk1,xik2,Jk2,Xi,La] = XiLaMaker(kx, ky, qx, qy, t, J, mu)

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

function [E0out, Duout, Ddout, correl, minEk_d] = tripletgapsolver(ek, mu, Delta, Jk, V, gk, gka)
Ek_u = sqrt((ek - mu + Jk) .^ 2 + Delta^2 * gka);    % p-wave
Ek_d = sqrt((ek - mu - Jk) .^ 2 + Delta^2 * gka);
Delta_u_k = V * Delta * gka / 2 ./ Ek_u;
Delta_d_k = V * Delta * gka / 2 ./ Ek_d;
E_0_k = (ek - mu) - Ek_u / 2 - Ek_d / 2 + Delta ^2 * gka ./ Ek_u /4 + Delta ^ 2 * gka ./ Ek_d / 4;
correl = Delta  * gk ./ Ek_u /2;
Duout = mean(Delta_u_k);
Ddout = mean(Delta_d_k);
E0out = mean(E_0_k);
minEk_d = min(Ek_d);
end

function [E0out, Dout,correl,minEk_d] = singletgapsolver(ek, mu, Delta, Jk, U, gk, gka)
Ek = sqrt((ek-mu).^2 + Delta^2*gka);
Ek_u = Ek + Jk;
Ek_d = Ek - Jk;
ind_u = Ek_u > 0;
ind_d = Ek_d > 0;
Delta_k = U * Delta * gka / 2 ./ Ek .* (ind_d - (~ind_u));
E_0_s_k = (ek - mu) - Ek + Delta^2 * gka ./ Ek / 2 .* (ind_u - (~ind_d))...
    + Ek_u .* (~ind_u) + Ek_d .* (~ind_d);
correl = Delta * gk ./ Ek / 2 .* (ind_u - (~ind_d));
Dout = mean(Delta_k);
E0out = mean(E_0_s_k);
minEk_d = min(Ek_d);
end
function [E_0, Dffout, Ek, correl] = FFgapsolver(Xi, La, Delta, U)
Ek = sqrt(Xi .^2 + Delta^2);
Ek_u = Ek + La;
Ek_d = Ek - La;

ind_u = Ek_u > 0;
ind_d = Ek_d > 0;
Delta_k = U * Delta /2 ./ Ek .* (ind_u - (~ind_d));

Dffout = mean(Delta_k);

E_0_k = Xi - Ek + Delta^2 ./ Ek /2 .* (ind_u -(~ind_d))...
    + Ek_u .* (~ind_u) + Ek_d .* (~ind_d);

correl = Delta ./ Ek / 2 .* (ind_u - (~ind_d));

E_0 = mean(E_0_k);


% E1ff=mean(Xi - Ek - Delta^2 ./ Ek /2 .* (ind_u -(~ind_d)));
% E2ff=mean(Ek_u .* (~ind_u) + Ek_d .* (~ind_d));
end

function draw_gapstructure(gap_1, gap_2, N)
figure
hold on
plot([gap_1(N/2+1:N,N/2+1);gap_1(1,N/2+1:N)';flip(diag(gap_1(N/2+1:N,N/2+1:N)));diag(gap_1(N/2+1:-1:2,N/2+1:N));gap_1(1:N/2+1,1);gap_1(N/2+1,N:-1:N/2+1)'],'r','LineWidth',1)
plot([gap_2(N/2+1:N,N/2+1);gap_2(1,N/2+1:N)';flip(diag(gap_2(N/2+1:N,N/2+1:N)));diag(gap_2(N/2+1:-1:2,N/2+1:N));gap_2(1:N/2+1,1);gap_2(N/2+1,N:-1:N/2+1)'],'b','LineWidth',1)
plot([-gap_1(N/2+1:N,N/2+1);-gap_1(1,N/2+1:N)';flip(diag(-gap_1(N/2+1:N,N/2+1:N)));diag(-gap_1(N/2+1:-1:2,N/2+1:N));-gap_1(1:N/2+1,1);-gap_1(N/2+1,N:-1:51)'],'r','LineWidth',1)
plot([-gap_2(N/2+1:N,N/2+1);-gap_2(1,N/2+1:N)';flip(diag(-gap_2(N/2+1:N,N/2+1:N)));diag(-gap_2(N/2+1:-1:2,N/2+1:N));-gap_2(1:N/2+1,1);-gap_2(N/2+1,N:-1:N/2+1)'],'b','LineWidth',1)
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

function draw_singlettoponodalstructure(mu,ek,Jk,kkx,N,gk)
figure
hold on
surf(kkx, kkx, reshape(gk,[N,N]), EdgeColor="none");
[~, h1]=contour(kkx, kkx, reshape(gk,[N,N]), [0,0], 'k', Linewidth=1.5);
[~, h2]=contour(kkx, kkx, reshape(ek-mu,[N,N]), [0,0], 'c', Linewidth=1.5);
[~, h3]=contour(kkx, kkx, reshape(ek-mu+Jk,[N,N]), [0,0], 'r', Linewidth=1.5);
[~, h4]=contour(kkx, kkx, reshape(ek-mu-Jk,[N,N]), [0,0], 'b', Linewidth=1.5);

h1.ContourZLevel=10;
h2.ContourZLevel=10;
h3.ContourZLevel=10;
h4.ContourZLevel=10;
box on
axis square
clim([-2,2]);
xticks([-pi,-pi/2,0,pi/2,kkx(end)])
yticks([-pi,-pi/2,0,pi/2,kkx(end)])
xlim([-pi,kkx(end)])
ylim([-pi,kkx(end)])
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
fontsize(14,"points")
end

function draw_FFtoponodalstructure(xi1k,xi1Jk1,xi2Jk2,kkx,N)
figure
hold on
surf(kkx, kkx, reshape(ones(1,N^2),[N,N]), EdgeColor="none");
[~, h0]=contour(kkx, kkx, reshape(xi1k,[N,N]), [0,0], 'c', Linewidth=1.5);
[~, h1]=contour(kkx, kkx, reshape(xi1Jk1,[N,N]), [0,0], 'r', Linewidth=1.5);
[~, h2]=contour(kkx, kkx, reshape(xi2Jk2,[N,N]), [0,0], 'b', Linewidth=1.5);

h0.ContourZLevel=10;
h1.ContourZLevel=10;
h2.ContourZLevel=10;
box on
axis square
clim([-2,2]);
xticks([-pi,-pi/2,0,pi/2,kkx(end)])
yticks([-pi,-pi/2,0,pi/2,kkx(end)])
xlim([-pi,kkx(end)])
ylim([-pi,kkx(end)])
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
fontsize(14,"points")
end

function draw_triplettoponodalstructure(mu,ek,Jk,kkx,N,gk)
figure
hold on
surf(kkx, kkx, reshape(real(gk),[N,N]), EdgeColor="none");
[~, h3]=contour(kkx, kkx, reshape(real(gk),[N,N]), [0,0], 'k', Linewidth=1.5);
[~, h4]=contour(kkx, kkx, reshape(ek-mu+Jk,[N,N]), [0,0], 'r', Linewidth=1.5);
h3.ContourZLevel=10;
h4.ContourZLevel=10;
box on
axis square
clim([-2,2]);
title('p_{x}-wave nodal structure')
xticks([-pi,-pi/2,0,pi/2,kkx(end)])
yticks([-pi,-pi/2,0,pi/2,kkx(end)])
xlim([-pi,kkx(end)])
ylim([-pi,kkx(end)])
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
fontsize(14,"points")

figure
hold on
surf(kkx,kkx, reshape(imag(gk),[N,N]), EdgeColor="none");
[~, h5]=contour(kkx, kkx, reshape(imag(gk),[N,N]), [0,0], 'k', Linewidth=1.5);
[~, h6]=contour(kkx, kkx, reshape(ek-mu+Jk,[N,N]), [0,0], 'r', Linewidth=1.5);
h5.ContourZLevel=10;
h6.ContourZLevel=10;
box on
axis square
clim([-2,2]);
title('p_{y}-wave nodal structure')
xticks([-pi,-pi/2,0,pi/2,kkx(end)])
yticks([-pi,-pi/2,0,pi/2,kkx(end)])
xlim([-pi,kkx(end)])
ylim([-pi,kkx(end)])
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
fontsize(14,"points")
end

function [ind1, ind2] = BFS_archiver(gap_1, gap_2)
ind1 = (gap_1 < 0) * 1;
ind1(ind1==0) = NaN;
ind2 = (gap_2 < 0) * 1;
ind2(ind2==0) = NaN;
end

function FermiDirac = nF(beta, E, mu)
FermiDirac = 1 / (1 + exp(beta * (E - mu)));
end