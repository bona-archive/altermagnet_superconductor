%% parameters
% N = 200;
% t = 0.5;
% JJ=2.5;
% U_ran=linspace(5,5,1);
% V_ran=linspace(0,5,100);
% J_ran = linspace(0, 1, 100);
% x = J_ran;
% xl = 'J';
% 
% Dend = 3;
% Dinput=linspace(0, Dend, Dend * 200);
% 
% clc
% % close all
% 
% mu = -1;
% tl = sprintf('t=%.3f, mu=%.3f', t, mu);
% 
% kx = linspace(-pi, pi, N+1);
% kx(end) = [];
% ky = linspace(-pi, pi, N+1);
% ky(end) = [];
% 
% fk=findBFS(kx,1,0.5,-1,0.792);
% 
% [kx, ky] = meshgrid(kx, ky);
% kx = reshape(kx, [numel(kx), 1]);
% ky = reshape(ky, [numel(ky), 1]);
% gk = sin(kx) + 1i * sin(ky);
% gka = abs(gk).^2;
% epsilon=1e-5;
% ek = e_k(kx, ky, 't', t);



%% s-wave (1d phase diagram)
% clc
% close all
clear all

N = 200;
t = 0.5;
mu = -1;
J=2.5;

U_ran=linspace(1,1,1);
% U=5;
J_ran = linspace(0, 1, 11);

Dend = 3;
Dinput=linspace(0, 3, Dend * 200);

% tl2 = sprintf('t=%.3f, mu=%.3f, U=%.3f,', t, mu,U);

kx = linspace(-pi, pi, N+1);
kx(end) = [];
ky = linspace(-pi, pi, N+1);
ky(end) = [];
[kx, ky] = meshgrid(kx, ky);
kx = reshape(kx, [numel(kx), 1]);
ky = reshape(ky, [numel(ky), 1]);

gk = sin(kx) + 1i * sin(ky);
gka = abs(gk).^2;
epsilon=1e-5;
ek = e_k(kx, ky, 't', t);
Jk = J_k(kx, ky, 'J', J);

Dsout = zeros(1,length(Dinput));
Dsfinal = zeros(1,length(U_ran));

figure
plot(Dinput, Dinput)
hold on
xlabel('\Delta_{s input}')
ylabel('\Delta_{s output}')
% title(tl2)

for u = 1: length(U_ran)
    for i = 1: length(Dinput)
        % J = J_ran(u);
        U = U_ran(u);
        % Jk = J_k(kx, ky, 'J', J);
        Delta = Dinput(i);
        %%% analytic
        Ek = sqrt((ek - mu).^2 + Delta^2);
        Ek_u = sqrt((ek - mu).^2 + Delta^2) + Jk;
        Ek_d = sqrt((ek - mu).^2 + Delta^2) - Jk;
        ind_u = Ek_u > 0;
        ind_d = Ek_d > 0;
        Delta_k = U * Delta / 2 ./ Ek .* (ind_u - (~ind_d));
        mean(Delta_k);
        Dsout(i) = mean(Delta_k);
        %%%% eig ED
        % pause(2)
    end
    plot(Dinput, Dsout)
    [~ , I] = min(abs(Dsout ./ (Dinput+epsilon) - 1));
    if (Dsout(2) - 0) / (Dinput(2) - Dinput(1)) < 1
        I = 1;
    else
        I = I;
    end
    Dsfinal(u) = Dsout(I); %%% self consistent Delta
end

E_0_s = zeros(1, length(Dsfinal));
for i = 1: length(Dsfinal)
    Delta = Dsfinal(i);
    % J = J_ran(i);
    % Jk = J_k(kx, ky, 'J', J);
    Ek = sqrt((ek - mu).^2 + Delta^2);
    Ek_u = sqrt((ek - mu).^2 + Delta^2) + Jk;
    Ek_d = sqrt((ek - mu).^2 + Delta^2) - Jk;
    ind_u = Ek_u > 0;
    ind_d = Ek_d > 0;
    E_0_s_k = (ek - mu) - Ek + Delta^2 ./ Ek / 2 .* (ind_u - (~ind_d))...
        + Ek_u .* (~ind_u) + Ek_d .* (~ind_d);
    E_0_s(i) = mean(E_0_s_k);
    %%% archive
    EEk_u(i,:)=Ek_u;
    EEk_d(i,:)=Ek_d;
end
figure
plot(U_ran, E_0_s, 'o')
% xlabel('J')
ylabel('E_{0s}')
title(tl2)

%% p-wave (1d phase diagram)
% clc
% close all

N = 200;
t = 0.5;
mu = -1;
% U_ran=linspace(5,5,1);
V_ran=linspace(0,3,10);
% V=5;
J_ran = linspace(0, 1, 100);
J=0.5;
Dend = 3;
Dinput=linspace(0, Dend, Dend * 200);

tl2 = sprintf('t=%.3f, mu=%.3f, J=%.3f,', t, mu,J);

kx = linspace(-pi, pi, N+1);
kx(end) = [];
ky = linspace(-pi, pi, N+1);
ky(end) = [];
[kx, ky] = meshgrid(kx, ky);
kx = reshape(kx, [numel(kx), 1]);
ky = reshape(ky, [numel(ky), 1]);

gk = sin(kx) + 1i * sin(ky);
gka = abs(gk).^2;
epsilon=1e-5;
ek = e_k(kx, ky, 't', t);
Jk = J_k(kx, ky, 'J', J);

Dpout_u = zeros(1, length(Dinput));
Dpout_d = zeros(1, length(Dinput));
Dpfinal_u = zeros(1, length(V_ran));
Dpfinal_d = zeros(1, length(V_ran));
E_0_p = zeros(1, length(V_ran));

figure
hold on
plot(Dinput, Dinput)
xlabel('\Delta_{p input}')
ylabel('\Delta_{p out}')
title(tl2)
for u = 1: length(V_ran)
    for i = 1: length(Dinput)
        % J = J_ran(u);
        % Jk = J_k(kx, ky, 'J', J);
        V = V_ran(u);
        Delta_u = Dinput(i);
        Delta_d = Dinput(i);
        %%% analytic
        Ek_u = sqrt((ek - mu + Jk).^2 + Delta_u^2 * gka);
        Ek_d = sqrt((ek - mu - Jk).^2 + Delta_d^2 * gka);
        Delta_u_k = V * Delta_u * gka / 2 ./ Ek_u;
        Delta_d_k = V * Delta_d * gka / 2 ./ Ek_d;
        Dpout_u(i) = mean(Delta_u_k);
        Dpout_d(i) = mean(Delta_d_k);
        %%%% eig ED
    end
    
    [~, I] = min(abs(Dpout_u ./ (Dinput + epsilon) - 1));
    if (Dpout_u(2) - 0) / (Dinput(2) - Dinput(1)) < 1
        I = 1;
    else
        I = I;
    end
    Dpfinal_u(u) = Dpout_u(I);

    % [~,I] = min(abs(Dpout_d./(Dinput+epsilon)-1));
    % if (Dpout_d(2)-0) / (Dinput(2)-Dinput(1)) < 1
    %     I = 1;
    % else
    %     I = I;
    % end
    Dpfinal_d(u) = Dpout_d(I);

    plot(Dinput, Dpout_u);
end

for i=1:length(Dpfinal_u)
    Delta_u=Dpfinal_u(i);
    Delta_d=Dpfinal_d(i);
    % J = J_ran(i);
    % Jk = J_k(kx, ky, 'J', J);
    Ek_u = sqrt((ek - mu + Jk).^2 + Delta_u^2 * gka);
    Ek_d = sqrt((ek - mu - Jk).^2 + Delta_d^2 * gka);
    E_0_p_k = (ek - mu) - Ek_u / 2 - Ek_d / 2 ...
        + Delta_u^2 * gka ./ Ek_u / 4 + Delta_d^2 * gka ./ Ek_d / 4;
    E_0_p(i) = mean(E_0_p_k);

    Epk_u(i,:)=Ek_u;
    Epk_d(i,:)=Ek_d;
end
figure
plot(V_ran, E_0_p)
xlabel('J')
ylabel('E_{0s}')
title(tl2)


figure
plot(V_ran,E_0_p)
hold on
plot(U_ran,E_0_s)

%% High symmetry point plot(s-wave)
Delta=Dsfinal(1);  % select Delta in (~)
Ek_u = sqrt((ek - mu).^2 + Delta^2) + Jk;
Ek_d = sqrt((ek - mu).^2 + Delta^2) - Jk;

EE1 = reshape(Ek_u, [N,N]);
EE2 = reshape(Ek_d, [N,N]);

figure
plot([EE1(101:200,101);EE1(1,101:200)';flip(diag(EE1(101:200,101:200)));diag(EE1(101:-1:2,101:200));EE1(1:101,1);EE1(101,200:-1:101)'],'k','LineWidth',1)
hold on
plot([EE2(101:200,101);EE2(1,101:200)';flip(diag(EE2(101:200,101:200)));diag(EE2(101:-1:2,101:200));EE2(1:101,1);EE2(101,200:-1:101)'],'k','LineWidth',1)
% plot([EE3(101:200,101);EE3(1,101:200)';flip(diag(EE3(101:200,101:200)));diag(EE3(101:-1:2,101:200));EE3(1:101,1);EE3(101,200:-1:101)'],'k','LineWidth',1)
% plot([EEE4(101:200,101);EEE4(1,101:200)';flip(diag(EEE4(101:200,101:200)));diag(EEE4(101:-1:2,101:200));EEE4(1:101,1);EEE4(101,200:-1:101)'],'k','LineWidth',1)
plot(zeros(1,601),'k')
xticks([0 N/2 N 3*N/2 2*N 5*N/2 3*N])
xticklabels({'\Gamma', 'X', 'M','\Gamma','M', 'Y', '\Gamma'})
xlim([0 3*N])
ylabel('E_{\uparrow/\downarrow}(k)')
ax1=gca;
ax1.XGrid='on';
ax1.GridAlpha=1;
ax1.FontName='Times New Roman';
ax1.FontSize=12;
% title(sprintf(str, t, J, mu, V, Delta))

%% High symmetry point plot(p-wave)
Delta=Dinput(~);  % select Delta in (~)
Ek_u = sqrt((ek - mu + Jk).^2 + Delta_u^2 * gka);
Ek_d = sqrt((ek - mu - Jk).^2 + Delta_d^2 * gka);

EE1 = reshape(Ek_u, [N,N]);
EE2 = reshape(Ek_d, [N,N]);

figure
plot([EE1(101:200,101);EE1(1,101:200)';flip(diag(EE1(101:200,101:200)));diag(EE1(101:-1:2,101:200));EE1(1:101,1);EE1(101,200:-1:101)'],'k','LineWidth',1)
hold on
plot([EE2(101:200,101);EE2(1,101:200)';flip(diag(EE2(101:200,101:200)));diag(EE2(101:-1:2,101:200));EE2(1:101,1);EE2(101,200:-1:101)'],'k','LineWidth',1)
% plot([EE3(101:200,101);EE3(1,101:200)';flip(diag(EE3(101:200,101:200)));diag(EE3(101:-1:2,101:200));EE3(1:101,1);EE3(101,200:-1:101)'],'k','LineWidth',1)
% plot([EEE4(101:200,101);EEE4(1,101:200)';flip(diag(EEE4(101:200,101:200)));diag(EEE4(101:-1:2,101:200));EEE4(1:101,1);EEE4(101,200:-1:101)'],'k','LineWidth',1)
plot(zeros(1,601),'k')
xticks([0 N/2 N 3*N/2 2*N 5*N/2 3*N])
xticklabels({'\Gamma', 'X', 'M','\Gamma','M', 'Y', '\Gamma'})
xlim([0 3*N])
ylabel('E_{\uparrow/\downarrow}(k)')
ax1=gca;
ax1.XGrid='on';
ax1.GridAlpha=1;
ax1.FontName='Times New Roman';
ax1.FontSize=12;
% title(sprintf(str, t, J, mu, V, Delta))
%% 2d phase diagram
clc
% close all

N = 100;
t = 0.5;
mu = -1;

% U_ran=linspace(5,5,1);
U=1;
ParamN=10;
V_ran=linspace(0,3,ParamN);
J_ran = linspace(0, 3, ParamN+1);
x = J_ran;
xl = 'J';
[JJ, VV] = meshgrid(J_ran, V_ran);
VV = reshape(VV,[numel(VV), 1]);
JJ = reshape(JJ, [numel(JJ), 1]);

Dend = 3;
Dinput=linspace(0, Dend, Dend * 200)';

tl = sprintf('t=%.3f, mu=%.3f', t, mu);

kx = linspace(-pi, pi, N+1);
kx(end) = [];
ky = linspace(-pi, pi, N+1);
ky(end) = [];
[kx, ky] = meshgrid(kx, ky);
kx = reshape(kx, [numel(kx), 1]);
ky = reshape(ky, [numel(ky), 1]);

gk = sin(kx) + 1i * sin(ky);
gka = abs(gk).^2;
epsilon=1e-5;
ek = e_k(kx, ky, 't', t);

Dsout = zeros(length(Dinput),1);
Dpout_u = zeros(length(Dinput), 1);
Dpout_d = zeros(length(Dinput), 1);
Dsfinal = zeros(length(J_ran), length(V_ran));
Dpfinal_u = zeros(length(J_ran), length(V_ran));
Dpfinal_d = zeros(length(J_ran), length(V_ran));
E_0_s=zeros(numel(Dsfinal), 1);
E_0_p=zeros(numel(Dsfinal), 1);

for j = 1: length(J_ran)
    for v = 1: length(V_ran)
        J = J_ran(j)
        Jk = J_k(kx, ky, 'J', J);
        V = V_ran(v)
        for d = 1: length(Dinput)
            Delta = Dinput(d);
            Delta_u = Dinput(d);
            Delta_d = Dinput(d);
            Ek = sqrt((ek - mu).^2 + Delta^2);
            Ek_u = sqrt((ek - mu).^2 + Delta^2) + Jk;
            Ek_d = sqrt((ek - mu).^2 + Delta^2) - Jk;
            ind_u = Ek_u > 0;
            ind_d = Ek_d > 0;
            Delta_k = U * Delta / 2 ./ Ek .* (ind_u - (~ind_d));
            Dsout(d) = mean(Delta_k);

            Ek_u = sqrt((ek - mu + Jk).^2 + Delta_u^2 * gka);
            Ek_d = sqrt((ek - mu - Jk).^2 + Delta_d^2 * gka);
            Delta_u_k = V * Delta_u * gka / 2 ./ Ek_u;
            Delta_d_k = V * Delta_d * gka / 2 ./ Ek_d;
            Dpout_u(d) = mean(Delta_u_k);
            Dpout_d(d) = mean(Delta_d_k);
            [~ , Is] = min(abs(Dsout ./ (Dinput+epsilon) - 1));
            if (Dsout(2) - 0) / (Dinput(2) - Dinput(1)) < 1
                Dsfinal(v,j) = 0;
            else
                Dsfinal(v,j) = Dsout(Is);
            end
            %%%% Bogoliubov Fermi surface
            % f = @(k) J^2 * sin(k)^4 - (2*t*cos(k) - mu)^2 - Dsout(Is)^2;
            %
            % k_range = linspace(-pi, pi, N+1);k_range(end)=[];
            % roots = [];
            %
            % for i = 1:length(k_range)-1
            %     try
            %         k_root = fzero(f, [k_range(i), k_range(i+1)]);
            %
            %         if k_root >= k_range(i) && k_root <= k_range(i+1)
            %             roots(end+1) = k_root;
            %         end
            %     catch
            %         %
            %     end
            % end
            % roots = unique(roots);
            [~, Ipu] = min(abs(Dpout_u ./ (Dinput + epsilon) - 1));
            if (Dpout_u(2) - 0) / (Dinput(2) - Dinput(1)) < 1
                Dpfinal_u(v,j) = 0;
            else
                Dpfinal_u(v,j) = Dpout_u(Ipu);
            end
            % Dpfinal_u(j,v) = Dpout_u(Ipu);

            [~, Ipd] = min(abs(Dpout_d ./ (Dinput + epsilon) - 1));
            if (Dpout_d(2) - 0) / (Dinput(2) - Dinput(1)) < 1
                Dpfinal_d(v,j) = 0;
            else
                Dpfinal_d(v,j) = Dpout_d(Ipd);
            end
            % Dpfinal_d(j,v) = Dpout_d(Ipd);
        end
    end
end

Dsfinal = reshape(Dsfinal, [numel(Dsfinal), 1]);
Dpfinal_u = reshape(Dpfinal_u, [numel(Dpfinal_u), 1]);
Dpfinal_d = reshape(Dpfinal_d, [numel(Dpfinal_d), 1]);

for i = 1: length(Dsfinal)
    Delta = Dsfinal(i);
    J = JJ(i);
    Jk = J_k(kx, ky, 'J', J);    
    Ek = sqrt((ek - mu).^2 + Delta^2);
    Ek_u = sqrt((ek - mu).^2 + Delta^2) + Jk;
    Ek_d = sqrt((ek - mu).^2 + Delta^2) - Jk;
    ind_u = Ek_u > 0;
    ind_d = Ek_d > 0;
    E_0_s_k = (ek - mu) - Ek + Delta^2 ./ Ek / 2 .* (ind_u - (~ind_d))...
        + Ek_u .* (~ind_u) + Ek_d .* (~ind_d);
    E_0_s(i) = mean(E_0_s_k);
end
E_0_s = reshape(E_0_s, [length(V_ran), length(J_ran)]);

for i = 1: length(Dpfinal_u)
    Delta_u = Dpfinal_u(i);
    Delta_d = Dpfinal_d(i);
    J = JJ(i);
    Jk = J_k(kx, ky, 'J', J);
    Ek_u = sqrt((ek - mu + Jk).^2 + Delta_u^2 * gka);
    Ek_d = sqrt((ek - mu - Jk).^2 + Delta_d^2 * gka);
    E_0_p_k = (ek - mu) - Ek_u / 2 - Ek_d / 2 ...
        + Delta_u^2 * gka ./ Ek_u / 4 + Delta_d^2 * gka ./ Ek_d / 4;
    E_0_p(i) = mean(E_0_p_k);
end
E_0_p = reshape(E_0_p, [length(V_ran), length(J_ran)]);
JJ = reshape(JJ, [length(V_ran), length(J_ran)]);


figure
hold on
xlabel('J')
ylabel('V')
title(tl)
surf(JJ, VV, E_0_s)
surf(JJ, VV, E_0_p)

figure
hold on
xlabel('J')
ylabel('V')
title(tl, ', E_{0s}-E_{0p}')
surf(JJ, VV, E_0_s - E_0_p)

allMatrices = cat(3, E_0_s, E_0_p);
[~, minIndex] = min(allMatrices, [], 3);
figure
surf(JJ,VV,minIndex)
xlabel('J')
ylabel('V')
title(tl)
view(2)

%% NaN?
Ek_u = sqrt((ek - mu + Jk).^2 + Dinput(2)^2 * gka);
Delta_u_k = V * Dinput(2) * gka / 2 ./ Ek_u;
a=Delta_u_k==0;
find(a==false)
%% functions

function ek = e_k(kx, ky, varargin)
t = 0.5;
for i = 1: 2: length(varargin)
    if strcmpi(varargin{i}, 't')
        t = varargin{i+1};
    else
        error('Unrecognized input argument');
    end
end
ek = - 2 * t * cos(kx) - 2 * t * cos(ky);
end

function Jk = J_k(kx, ky, varargin)
J = 1;
for i = 1: 2: length(varargin)
    if strcmpi(varargin{i}, 'J')
        J = varargin{i+1};
    else
        error('Unrecognized input argument');
    end
end
Jk = J * sin(kx) .* sin(ky);
end

function fk=findBFS(kx,J,t,mu,Delta)
    fk=J^2 * sin(kx) .^4 - abs(4 * t * cos(kx) - mu) .^2 - abs(Delta)^2;
end