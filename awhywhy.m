clear all

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

gs = 1;
gsa = 1;
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

            [Ea, Da,~,mina] = singletgapsolver(ek, mu, Delta, Jk, U, gs, gsa);
           

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

%% FF phase
disp('start FF phase')

mu = mu_ran(1);
% for j = 11: 11
for j = 1: length(J_ran)
    J = J_ran(j);
    Jk = J .* sin(kx) .* sin(ky);
    
    % for qq = 4: 4
    for qq = 1: length(qyline)
        [Xi,La] = XiLaMaker(kx, ky, qx0line(qq), qyline(qq), t, J, mu);

        for i = 1: length(Dinput)
        

            Delta = Dinput(i);

            [EEq, DDq, ~, ~, ~, ~] = FFgapsolver(Xi, La, Delta, U);
            E0ffout(i) = EEq;
            Dffout(i) = DDq;

        end


        [~ , I] = min(abs(Dffout ./ (Dinput + epsilon) - 1));
        if (Dffout(2) - 0) / (Dinput(2) - Dinput(1)) < 1
            I = 1;
            Dff(qq) = 0;
            E0ff(qq) = E0ffout(I);
            E0ff_arc(j,qq) = E0ffout(I);
        else
            Dff(qq) = Dinput(I);
            E0ff(qq) = E0ffout(I);
            E0ff_arc(j,qq) = E0ffout(I);
        end

    end

    [aa, ~] = min(E0ff);
    II = find(E0ff == aa);
    II = II(1);
    Dff_self(j) = Dff(II);
    E0ff_self(j) = aa;
    q_self(j) = qyline(II);
    [Xi_self,La_self] = XiLaMaker(kx, ky, 0, qyline(II), t, J, mu);
    % [~, ~, ca, ~, ~, ~] = FFgapsolver(Xi_self, La_self, Dff(II), U);
    % correlff_self(j, :,:) = reshape(ca,[N,N]);
end
E0ff_self = repmat(E0ff_self, length(V_ran), 1);

%% p-, es-, d- wave

disp('start p-, es-, d-wave')
for v = 1: length(V_ran)

    
    V = V_ran(v);


    for j = 1: length(J_ran)


        J = J_ran(j);
        Jk = J .* sin(kx) .* sin(ky);


        for i = 1: length(Dinput)

            Delta = Dinput(i);

            [Ea, Dua, Dda, ~, ~] = tripletgapsolver(ek, mu, Delta, Jk, V, gpk, gpka); % p wave

            Dpout_u(i) = Dua;
            Dpout_d(i) = Dda;
            E0pout(i) = Ea;
            Dpout_u_arc(v, j, i) = Dua;
            Dpout_d_arc(v, j, i) = Dda;
            E0pout_arc(v, j, i) = Ea;


            [Eb, Db, ~, ~] = singletgapsolver(ek, mu, Delta, Jk, V, gesk, geska); % es wave
            
            Desout(i) = Db;
            E0esout(i) = Eb;
            Desout_arc(v, j, i) = Db;
            E0esout_arc(v, j, i) = Eb;

            [Ec, Dc, ~, ~] = singletgapsolver(ek, mu, Delta, Jk, V, gdk, gdka); % d wave

            Ddout(i) = Dc;
            E0dout(i) = Ec;
            Ddout_arc(v, j, i) = Dc;
            E0dout_arc(v, j, i) = Ec;

        end


        [~ , I] = min(abs(Dpout_u ./ (Dinput + epsilon) - 1));
        if (Dpout_u(2) - 0) / (Dinput(2) - Dinput(1)) < 1
            I = 1;
            Dp = 0;
            E0p = E0pout(I);
        else
            Dp = Dinput(I);
            E0p = E0pout(I);
        end

        [~ , I] = min(abs(Desout ./ (Dinput + epsilon) - 1));
        if (Desout(2) - 0) / (Dinput(2) - Dinput(1)) < 1
            I = 1;
            Des = 0;
            E0es = E0esout(I);
        else
            Des = Dinput(I);
            E0es = E0esout(I);
        end

        [~ , I] = min(abs(Ddout ./ (Dinput + epsilon) - 1));
        if (Ddout(2) - 0) / (Dinput(2) - Dinput(1)) < 1
            I = 1;
            Dd = 0;
            E0d = E0dout(I);
        else
            Dd = Dinput(I);
            E0d = E0dout(I);
        end

        Dp_self(v,j) = Dp;
        E0p_self(v,j) = E0p;
        % [~, ~, ~, correlp, ~] = tripletgapsolver(ek, mu, Dp, Jk, V, gpk, gpka);
        % correlp_self(v,j,:,:) = reshape(correlp,[N,N]);

        Des_self(v,j) = Des;
        E0es_self(v,j) = E0es;
        % [~, ~, correles, ~] = singletgapsolver(ek, mu, Des, Jk, V, gesk, geska);
        % correles_self(v,j,:,:) = reshape(correles,[N,N]);
        
        Dd_self(v,j) = Dd;
        E0d_self(v,j) = E0d;
        % [~, ~, correld, ~] = singletgapsolver(ek, mu, Dd, Jk, V, gdk, gdka);
        % correld_self(v,j,:,:) = reshape(correld,[N,N]);

    end
end

%% plot ground state energy

% custom colormap
customCMap = [
    1, 0.83, 0.83;  % Color near 1: #F0D4D4
    1, 0.90, 0.80;  % Color near 2: #FFE6CC
    1, 1, 1;        % Color near 3: White
    0.90, 0.86, 0.86; % Color near 4: #E6DBDC
    0.90, 0.94, 1;  % Color near 5: #E6EFFF
];

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


%% gap structure
vv = 23;% V_ran(23)=1.00 in V:0.6~1.5 version
jjs = 11;% J_ran(11)=0.2, J_ran(16)=0.3, J_ran(21)=0.41
jjff = 21;
jjp = 40;
jj=11;
Jk = J_ran(jj) * sin(kx) .* sin(ky);
mu=-1;
kx1 = kx;
kx2 = -kx;

ky1 = ky + q_self(jjff);
ky2 = -ky + q_self(jjff);
Jk1 = J_ran(jjff) * sin(kx1) .* sin(ky1);
Jk2 = J_ran(jjff) * sin(kx2) .* sin(ky2);
xik1 = - 2 * t * cos(kx1) - 2 * t * cos(ky1) - mu;
xik2 = - 2 * t * cos(kx2) - 2 * t * cos(ky2) - mu;
Xi = .5 * (xik1 + Jk1 +xik2 - Jk2);
La = .5 * (xik1 + Jk1 -xik2 + Jk2);


gap_s1 = reshape(sqrt((ek - mu) .^2 + Ds_self(jj) ^2) + J_ran(jjs) * sin(kx) .* sin(ky), [N, N]);
gap_s2 = reshape(sqrt((ek - mu) .^2 + Ds_self(jj) ^2) - J_ran(jjs) * sin(kx) .* sin(ky), [N, N]);
gap_ff1 = reshape(sqrt(Xi .^2 + Dff_self(jjff)^2) + La, [N, N]);
gap_ff2 = reshape(sqrt(Xi .^2 + Dff_self(jjff)^2) - La, [N, N]);

gap_d1 = reshape(sqrt((ek - mu) .^2 + Dd_self(vv, jj) ^2 * gdka) + Jk, [N, N]);
gap_d2 = reshape(sqrt((ek - mu) .^2 + Dd_self(vv, jj) ^2 * gdka) - Jk, [N, N]);
gap_es1 = reshape(sqrt((ek - mu) .^2 + Des_self(vv, jj) ^2 * geska) + Jk, [N, N]);
gap_es2 = reshape(sqrt((ek - mu) .^2 + Des_self(vv, jj) ^2 * geska) - Jk, [N, N]);

gap_p1 = reshape(sqrt((ek - mu + J_ran(jjp) * sin(kx) .* sin(ky)) .^ 2 + Dp_self(vv, jjp) ^2 * gpka), [N, N]);
gap_p2 = reshape(sqrt((ek - mu - J_ran(jjp) * sin(kx) .* sin(ky)) .^ 2 + Dp_self(vv, jjp) ^2 * gpka), [N, N]);



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
%% function

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
Dout = mean(Delta_k, 'omitnan');
E0out = mean(E_0_s_k, 'omitnan');
minEk_d = min(Ek_d);
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
E_0 = mean(E_0_k, 'omitnan');

end


function [E0out, Duout, Ddout, correl, minEk_d] = tripletgapsolver(ek, mu, Delta, Jk, V, gk, gka)
Ek_u = sqrt((ek - mu + Jk) .^ 2 + Delta^2 * gka);    % p-wave
Ek_d = sqrt((ek - mu - Jk) .^ 2 + Delta^2 * gka);
Delta_u_k = V * Delta * gka / 2 ./ Ek_u;
Delta_d_k = V * Delta * gka / 2 ./ Ek_d;
E_0_k = (ek - mu) - Ek_u / 2 - Ek_d / 2 + Delta ^2 * gka ./ Ek_u /4 + Delta ^ 2 * gka ./ Ek_d / 4;
correl = Delta  * gk ./ Ek_u /2;
Duout = mean(Delta_u_k, 'omitnan');
Ddout = mean(Delta_d_k, 'omitnan');
E0out = mean(E_0_k, 'omitnan');
minEk_d = min(Ek_d);
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