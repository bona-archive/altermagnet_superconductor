% clear all

N = 300;
ParamN = 50;
t = 0.5;
mu = -1.2;
mu_ran = ones(1,1)*(mu);
U_ran = ones(1,1) * 1.5;
U = U_ran(1);
J_ran = linspace(0, 1, ParamN);
% J_ran = ones(1,1) * 0.4;
V_ran = linspace(0.6, 1.5, ParamN);
V_ran = linspace(0, 1.5, ParamN);

Dend = 1;
Dinput = linspace(0, Dend, 500);


kx = linspace(-pi, pi, N+1);
ky = linspace(-pi, pi, N+1);
kx(end) = [];
ky(end) = [];
[kx, ky] = meshgrid(kx, ky);
kx = reshape(kx, [numel(kx), 1]);
ky = reshape(ky, [numel(ky), 1]);

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

ek = - 2* t * cos(kx) - 2 * t * cos(ky);
%% s-wave
disp('start s-wave')
figure
hold on
plot(Dinput,Dinput)
for mmu = 1: length(mu_ran)
    mu = mu_ran(mmu);
    for j = 1: length(J_ran)

        J = J_ran(j);
        Jk = J .* sin(kx) .* sin(ky);

        for i = 1: length(Dinput)

            Delta = Dinput(i);

            [Ea, Da,~,mina] = singletgapsolver(ek, mu, Delta, Jk, U);
           

            Dsout(i) = Da;
            E0sout(i) = Ea;
            minEk_d(i) = mina;
            Dsout_arc(j, i) = Da;
            E0sout_arc(j, i) = Ea;

        end
        plot(Dinput,Dsout)

        [~ , I] = min(abs(Dsout ./ (Dinput + epsilon) - 1));
        if (Dsout(2) - 0) / (Dinput(2) - Dinput(1)) < 1
            I = 1;
            Ds = 0;
            E0s = E0sout(I);
            BFS = minEk_d(I);
        else
            Ds = Dinput(I);
            E0s = E0sout(I);
            BFS = minEk_d(I);
        end

        I_self(j) = I;
        Ds_self(j) = Ds;
        E0s_self(j) = E0s;
        BFSfinder(j) = BFS;
        % [~, ~, correla,~] = singletgapsolver(ek, mu, Ds, Jk, U);
        % correls_self(j,:,:) = reshape(correla,[N,N]);
    end
E0s_self = repmat(E0s_self, length(V_ran), 1);
BFSfinder = repmat(BFSfinder, length(V_ran), 1);
end
%% check?
[Ea, Da,~,mina] = singletgapsolver(ek, mu, Delta, Jk, U);
%% quasiparticle band structure and correlation and correlation suppressed area
ii=1;

figure
[Ek_u_test, Ek_d_test] = signeltEkJkMaker(kx, ky, ek, mu, Ds_self(ii), J_ran(ii),1);

surf(reshape(Ek_u_test,[N,N]))
hold on
surf(reshape(Ek_d_test,[N,N]))
zlim([-0.1,0]);view(2)
figure
surf(squeeze(correls_self(ii,:,:)));view(2)

for ii = 1: length(J_ran)
    [Ek_u_test, Ek_d_test] = signeltEkJkMaker(kx, ky, ek, mu, Ds_self(ii), J_ran(ii),1);

    Area_s(ii) = sum((Ek_u_test < 0) * 1) + sum((Ek_d_test < 0) * 1);
end
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

function [E0out, Dout,correl,minEk_d] = singletgapsolver_finiteTemp(ek, mu, Delta, Jk, U, T)
Ek = sqrt((ek-mu).^2 + Delta^2);
Ek_u = Ek + Jk;
Ek_d = Ek - Jk;
ind_u = Ek_u > 0;
ind_d = Ek_d > 0;
Delta_k = U * Delta / 4 ./ Ek .* (tanh(Ek_u/ 2 / T) + tanh(Ek_d / 2/ T));
E_0_s_k = (ek - mu) - Ek + Delta^2 ./ Ek / 4 .* (tanh(Ek_u/ 2 / T) + tanh(Ek_d / 2/ T))...
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

function [ind1, ind2] = BFS_archiver(gap_1, gap_2)
ind1 = (gap_1 < 0.0001 & gap_1 < 0.001) * 1;
ind1(ind1==0) = NaN;
ind2 = (gap_2 < 0.0001 & gap_2 < 0.001) * 1;
ind2(ind2==0) = NaN;
end