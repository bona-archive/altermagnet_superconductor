%% Superfluid density check for mu/2t=--1 (BFS in s-wave, FF state might induce negative superfluid density)
% clear all
% % load('mu=-1 real.mat')
% load('s_mu=-1_N1000.mat')
% load('ff_mu-1_N1000.mat')
% load('E0p_self.mat')
%% 실수는 성공의 어머니이다....
q_self_corrected = zeros(size(q_self)); % 수정된 q_self 배열
% Dff_self_corrected = zeros(size(Dff_self));

% q_self의 각 값에 대해 qyline에서 일치하는 인덱스를 찾음
for jj = 1:length(q_self)
    % q_self(jj) 값이 qyline에서 어느 인덱스에 해당하는지 찾음
    II = find(qyline == q_self(jj)); 
    % III = find(Dff == Dff_self(jj));
    % 만약 여러 개의 일치하는 값이 있을 경우 첫 번째 값을 사용
    if ~isempty(II)
        II = II(1);  % 첫 번째 일치하는 인덱스 선택
        % III = III(1);
        % 올바른 인덱스를 기반으로 10*(II-1)+1을 계산하여 qyline에서 값 가져오기
        q_self_corrected(jj) = qyline(10*(II-1) + 1);
        % Dff_self_corrected(jj) = Dff(10*(III-1) + 1);
    else
        warning('qyline에서 일치하는 값을 찾지 못했습니다: q_self(%d)', jj);
    end
end
q_self = q_self_corrected;
%% Parameter

% clear all

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

betaend = 50;
betainput = linspace(1, betaend, betaend);
betainput = betainput * 0.001;

Tempinput = 1./betainput;

Dend = 1;
Dinput = linspace(0, Dend, 500);


kx = linspace(-pi, pi, N+1);
ky = linspace(-pi, pi, N+1);
kx(end) = [];
ky(end) = [];
[kx, ky] = meshgrid(kx, ky);
% kx = reshape(kx, [numel(kx), 1]);
% ky = reshape(ky, [numel(ky), 1]);

qx = linspace(-pi, pi, N+1);qx(end)=[];
qy = qx;
[qx, qy] = meshgrid(qx, qy);
qx = reshape(qx, [numel(qx), 1]);
qy = reshape(qy, [numel(qy), 1]);

%% pre-calculation for q_self and Dff_self needs here
for i = 1:length(J_ran)
    if pi/2 < q_self(i) && q_self(i) < pi
        q_self(i) = -1 * (q_self(i)-pi);
    else
    end

end
%% Calculation

% J_cLeft = 0.3061(16th in J_ran), J_cRight = 0.4082(21st in J_ran)

beta = 100; % adequately big number for finite temperature?
figure
for ii = 1:21
    
    normal_up = -2 * t * (cos(kx) + cos(ky)) - mu + J_ran(ii) .* sin(kx) .* sin(ky);
    normal_down = -2 * t * (cos(kx) + cos(ky)) - mu - J_ran(ii) .* sin(kx) .* sin(ky);
    
    n_up = normal_density(normal_up, N);
    n_down = normal_density(normal_down, N);
    n_up_arc(ii) = n_up;
    n_down_arc(ii) = n_down;
    

    [~, ~, ~, ~, Xi_test, La_test] = XiLaMaker(kx, ky, 0, q_self(ii), t, J_ran(ii), mu);
    EFF_kupq = sqrt(Xi_test .^2 + Dff_self(ii) ^2) + La_test;
    surf(EFF_kupq,'EdgeColor','none')
    view(2)
    pause(0.01)
    superfluid_density_up(ii) = n_up + (mu/N/N) * sum(sum(beta .* nF(beta, EFF_kupq, 0) .* (1 - nF(beta, EFF_kupq, 0))));
    EFF_kdownq = sqrt(Xi_test .^2 + Dff_self(ii) ^2) - La_test;
    superfluid_density_down(ii) = n_down + (mu/N/N) * sum(sum(beta .* nF(beta, EFF_kdownq, 0) .* (1 - nF(beta, EFF_kdownq, 0))));
    condensation_density_down(ii) = (mu/N/N) * sum(sum(beta .* nF(beta, EFF_kdownq, 0) .* (1 - nF(beta, EFF_kdownq, 0))));

end

%s-wave J_ran(5)
%ff-phase J_ran(15)

% for ii = 1:length(betainput)
%     beta = betainput(ii);
%     temp = Tempinput(ii);
%     normal_up_betas = -2 * t * (cos(kx) + cos(ky)) - mu + J_ran(5) .* sin(kx) .* sin(ky);
%     normal_up_betaff = -2 * t * (cos(kx) + cos(ky)) - mu + J_ran(15) .* sin(kx) .* sin(ky);
% 
%     n_up_betas = normal_density(normal_up_betas, N);
%     n_up_betaff = normal_density(normal_up_betaff, N);
% 
%     [~, ~, ~, ~, Xi_test, La_test] = XiLaMaker(kx, ky, 0, q_self(15), t, J_ran(15), mu);
%     EFF_kupq = sqrt(Xi_test .^2 + Dff_self(15) ^2) + La_test;
%     superfluid_density_up_betaff(ii) = n_up_betas + (mu/N/N) * sum(sum(beta .* nF(beta, EFF_kupq, 0) .* (1 - nF(beta, EFF_kupq, 0))));
% 
%     [~, ~, ~, ~, Xi_test, La_test] = XiLaMaker(kx, ky, 0, q_self(5), t, J_ran(5), mu);
%     EFF_kupq = sqrt(Xi_test .^2 + Dff_self(5) ^2) + La_test;
%     superfluid_density_up_betas(ii) = n_up_betaff + (mu/N/N) * sum(sum(beta .* nF(beta, EFF_kupq, 0) .* (1 - nF(beta, EFF_kupq, 0))));
% end
%% figure

figure
plot(J_ran(1:21),superfluid_density_down)
hold on
plot(J_ran(1:21), zeros(1,21))
xlabel('J/2t')
ylabel('n_{s\downarrow}')

figure
plot(J_ran(1:21),n_down_arc)
xlabel('J/2t')
ylabel('n_\downarrow(J)')

figure
plot(J_ran(1:21),condensation_density_down)
xlabel('J/2t')
ylabel('n_{condensation, \downarrow}(J)')

figure
plot(J_ran(1:21),superfluid_density_up)
hold on
plot(J_ran(1:21), zeros(1,21))
xlabel('J/2t')
ylabel('n_{s\uparrow}')

figure
plot(betainput,superfluid_density_up_betas)
hold on
plot(betainput, zeros(1,length(betainput)))
xlabel('\beta=1/k_{B}T')
ylabel('n_{s\uparrow}, s-wave')

figure
plot(betainput,superfluid_density_up_betaff)
hold on
plot(betainput, zeros(1,length(betainput)))
xlabel('\beta=1/k_{B}T')
ylabel('n_{s\uparrow}, FF phase')

figure
plot(Tempinput,superfluid_density_up_betas)
hold on
plot(Tempinput, zeros(1,length(Tempinput)))
xlabel('Temperature')
ylabel('n_{s\uparrow}, s-wave')

figure
plot(Tempinput,superfluid_density_up_betaff)
hold on
plot(Tempinput, zeros(1,length(Tempinput)))
xlabel('Temperature')
ylabel('n_{s\uparrow}, FF phase')
%% function

function FermiDirac = nF(beta, E, mu)
FermiDirac = 1 ./ (1 + exp(beta .* (E - mu)));
end

function n_density = normal_density(normal, N)
    normal_neg = normal < 0;
    normal_occupied = normal_neg * 1;
    % figure
    % surf(reshape(normal_occupied,[N,N]))
    % view(2)
    n_density = sum(sum(normal_occupied))/N/N;
end

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

E_0 = mean(E_0_k);


% E1ff=mean(Xi - Ek - Delta^2 ./ Ek /2 .* (ind_u -(~ind_d)));
% E2ff=mean(Ek_u .* (~ind_u) + Ek_d .* (~ind_d));
end