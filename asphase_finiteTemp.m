clear all
load('C:\Users\hsb83\Desktop\ff_1000_V0to1.5.mat')
load('C:\Users\hsb83\Desktop\s_1000_V0to1.5.mat')
%% q_self corrected..
q_self_corrected = zeros(size(q_self)); % 수정된 q_self 배열
% Dff_self_corrected = zeros(size(Dff_self));

% q_self의 각 값에 대해 qyline에서 일치하는 인덱스를 찾음
for jj = 1:length(q_self)
    jj
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
%%
for i = 1:length(J_ran)
    if pi/2 < q_self(i) && q_self(i) < pi
        q_self(i) = -1 * (q_self(i)-pi);
    else
    end

end
figure
plot(q_self)

%%
Ere = E0ff_arc(15:20,:);
Ere = Ere(:,[251:500,1:250]);
kkkx = linspace(-pi/4,pi/4,125+1);kkkx(end)=[];
figure;surf(kkkx,J_ran(15:20),Ere(:,189:313),'Edgecolor','none')
for i = 15: 20
hold on
scatter(q_self(i),J_ran(i),'k','filled')
end
for i = 15: 20
hold on
scatter(-q_self(i),J_ran(i),'k','filled')
end
axis square

% figure
% plot(E0ff_arc(20,:))

%%
ParamN = 50;
t = 0.5;
mu = -1;
V = 1.5;
J_ran = linspace(0, 1, ParamN);
% T_ran = linspace(10^(-3), 0.12, 200); % s-wave = 0.116 FF = 0.00124
T_ran = linspace(0, 3, 5); % s-wave = 0.116 FF = 0.00124
% 구간: 밀도가 높은 부분
DN = 10;
DN_dense = 90;
% D1end = 0.02; % s-wave는 0.2 FF = 0.027
% D2end = 0.03;
D1end = 4;
D2end = 5;
D1 = linspace(0, D1end, DN+DN_dense);  % 예: 50개의 샘플

% 0.02 ~ 0.3 구간: 밀도가 낮은 부분
D2 = linspace(D1end, D2end, DN+1); % 예: 50개의 샘플

% 두 구간 병합
Dinput = [D1, D2(2:end)]; % 중복되는 값 제거

Dinput = linspace(0,D2end, 100);

kx = linspace(-pi, pi, N+1);
ky = linspace(-pi, pi, N+1);
kx(end) = [];
ky(end) = [];
qx_scan = kx;

[kx, ky] = meshgrid(kx, ky);
kx = reshape(kx, [numel(kx), 1]);
ky = reshape(ky, [numel(ky), 1]);

epsilon = 1e-6;

gpk = sin(kx) + 1i * sin(ky);
gpka = abs(gpk) .^ 2;
gesk = cos(kx) + cos(ky);
geska = abs(gesk) .^ 2;
gdk = cos(kx) - cos(ky);
gdka = abs(gdk) .^ 2;
ek = - 2* t * cos(kx) - 2 * t * cos(ky);
Jk = J * sin(kx) .* sin(ky);

Dout1 = zeros(1, length(Dinput));
E0out1 = zeros(1, length(Dinput));
Dout2 = zeros(1, length(Dinput));
E0out2 = zeros(1, length(Dinput));
Dout3 = zeros(1, length(Dinput));
E0out3 = zeros(1, length(Dinput));
Fq = zeros(1, 3);
rho = zeros(1, length(T_ran));

%% s-wave
disp('start s-wave')

dqx = qx_scan(2)-qx_scan(1);
jjj=1;
% jjj = jjff;

figure
hold on
plot(Dinput,Dinput)
for j = 1: length(T_ran)
    j
    T = T_ran(j);
    J = J_ran(jjj);
    % Delta = Dff_self(jjj);
    qx = q_self(jjj);


    for i = 1: length(Dinput)
        Delta = Dinput(i);

        J = J_ran(jjj);

        [Xi, La] = XiLaMaker(kx, ky, qx-dqx, 0, t, J, mu);
        [Ea1, Da1, ~, ~] = FFgapsolver_finiteTemp(Xi, La, Delta, V, T);

        Dout1(i) = Da1;
        E0out1(i) = Ea1;

    end

    [~ , I] = min(abs(Dout1 ./ (Dinput + epsilon) - 1));
    if (Dout1(2) - 0) / (Dinput(2) - Dinput(1)) < 1
        I = 1;
        Dsq1 = 0;
        E0sq1 = E0out1(I);
    else
        Dsq1 = Dinput(I);
        E0sq1 = E0out1(I);
    end
    Fq(1) = E0sq1;

    % [Xi, La] = XiLaMaker(kx, ky, qx-dqx, 0, t, J, mu);
    % [Ea1, ~, ~] = FFgapsolver_finiteTemp(Xi, La, Delta, V, T);
    % Fq(1) = Ea1;
    for i = 1: length(Dinput)
        Delta = Dinput(i);

        J = J_ran(jjj);

        [Xi, La] = XiLaMaker(kx, ky, qx, 0, t, J, mu);
        [Ea2, Da2, ~ ,~] = FFgapsolver_finiteTemp(Xi, La, Delta, V, T);

        Dout2(i) = Da2;
        E0out2(i) = Ea2;

    end
    plot(Dinput,Dout1)

    [~ , I] = min(abs(Dout2 ./ (Dinput + epsilon) - 1));
    if (Dout2(2) - 0) / (Dinput(2) - Dinput(1)) < 1
        I = 1;
        Dsq2 = 0;
        E0sq2 = E0out2(I);
    else
        Dsq2 = Dinput(I);
        E0sq2 = E0out2(I);
    end
    Fq(2) = E0sq2;

    % [Xi, La] = XiLaMaker(kx, ky, qx, 0, t, J, mu);
    % [Ea2, Da2, ~] = FFgapsolver_finiteTemp(Xi, La, Delta, V, T);
    % Fq(2) = Ea2;
    for i = 1: length(Dinput)
        Delta = Dinput(i);

        J = J_ran(jjj);

        [Xi, La] = XiLaMaker(kx, ky, qx+dqx, 0, t, J, mu);
        [Ea3, Da3, ~, ~] = FFgapsolver_finiteTemp(Xi, La, Delta, V, T);

        Dout3(i) = Da3;
        E0out3(i) = Ea3;

    end

    [~ , I] = min(abs(Dout2 ./ (Dinput + epsilon) - 1));
    if (Dout3(2) - 0) / (Dinput(2) - Dinput(1)) < 1
        I = 1;
        Dsq3 = 0;
        E0sq3 = E0out3(I);
    else
        Dsq3 = Dinput(I);
        E0sq3 = E0out3(I);
    end
    Fq(3) = E0sq3;
    % [Xi, La] = XiLaMaker(kx, ky, qx+dqx, 0, t, J, mu);
    % [Ea3, ~, ~] = FFgapsolver_finiteTemp(Xi, La, Delta, V, T);
    % Fq(3) = Ea3;

    rho(j) = (Fq(3) - 2 * Fq(2) + Fq(1)) / dqx / dqx;
    DeltaT(j) = Dsq2;
    
end
%% 

for i = 1: length(DeltaT)
    i
    Delta = DeltaT(i);
    TT = T_ran(i);
    en1 = sqrt(ek.^2 + Delta.^2) + Jk;
    en2 = sqrt(ek.^2 + Delta.^2) - Jk;
    en = [en1 en2 -en1 -en2];
    K_0T(i) = calculate_PenetrationDepth_FD(TT, en, Delta);
end

%% 

lam = 1 ./ sqrt(abs(rho));

Tcindex = find(DeltaT == 0);
Tcindex = Tcindex(1);
Tc = T_ran(Tcindex);

figure
plot(T_ran/Tc, DeltaT);xlabel('T/T_{c} [K]');ylabel('\Delta(T)')

figure
plot(T_ran(1:Tcindex-1)/Tc, rho(1:Tcindex-1)/rho(1));xlabel('T/T_{c} [K]');ylabel('K(0,T)/K(0, 0)')


figure
plot(T_ran(1:Tcindex-1)/Tc, lam(1:Tcindex-1) / lam(1));xlabel('T/T_{c} [K]');ylabel('\lambda(T)/\lambda(0)')
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

function [E_0, Dffout,Ek, F] = FFgapsolver_finiteTemp(Xi, La, Delta, U, T)
Ek = sqrt(Xi .^2 + Delta^2);
Ek_u = Ek + La;
Ek_d = Ek - La;

Ek_merged = [-Ek_u;-Ek_d];
% ind_u = Ek_u > 0;
% ind_d = Ek_d > 0;
Delta_k = U * Delta /4 ./ Ek .* (tanh(Ek_u / 2 / T) + tanh(Ek_d / 2 / T));

Dffout = mean(Delta_k);

if T == 0
    F = Ek_u.*1./(1+exp(Ek_u/T))+Ek_d.*1./(1+exp(Ek_d/T))+(-Ek_u).*1./(1+exp(-Ek_u/T))+(-Ek_d).*1./(1+exp(-Ek_d/T));
    F = mean(F, 'omitnan');
else
    c = max(-Ek_merged/T);
    logZ = c + log(sum(exp(-Ek_merged/T-c)));
    F = -T * logZ;
end
% E_0_k = Xi - Ek + Delta^2 ./ Ek /4 .* (tanh(Ek_u / 2 / T) + tanh(Ek_d / 2 /T))...
%     + Ek_u.*1./(1+exp(Ek_u/T))+Ek_d.*1./(1+exp(Ek_d/T))+(-Ek_u).*1./(1+exp(-Ek_u/T))+(-Ek_d).*1./(1+exp(-Ek_d/T));
E_0_k = Ek_u.*1./(1+exp(Ek_u/T))+Ek_d.*1./(1+exp(Ek_d/T))+(-Ek_u).*1./(1+exp(-Ek_u/T))+(-Ek_d).*1./(1+exp(-Ek_d/T));
E_0 = mean(E_0_k, 'omitnan');

end


function K_0T = calculate_PenetrationDepth_FD(t, en, Delta)
    K_0T = zeros(1, length(t));
    % 양수 에너지 (혹은 반대로 점유에너지) 만 골라내기
    en = en(en > 0);
    for i=1:length(t)
        % lambda(i) = 1 + 2 * sum(exp(-1 / t(i) * en) ./ ((exp(-1 / t(i) * en) + 1) .^ 2) .* en ./ sqrt(en.^2 - Delta^2));
        K_0T(i) = 1 + 2 * sum(1/t(i) .* exp(1 / t(i) * en) ./ ((exp(1 / t(i) * en) + 1) .^ 2 ));
        % eq0(i)=sum(en.*(1./(exp(en./t(i))+1))-en.*(1./(exp(-en./t(i))+1)));
    end
    % plot(t/t(end), K_0T/K_0T(1))
    % % plot(t,1./sqrt(eq0))
    % hold on
    % xlabel('T/T_{c} [K]')
    % ylabel('K(0, T)/K(0, 0) [a.u.]')
    % fontsize(16,"points")
    
end