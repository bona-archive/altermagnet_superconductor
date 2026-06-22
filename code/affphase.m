% clear all

N = 300;
ParamN = 50;
t = 0.5;
mu = -1.8;
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

E0ffout = zeros(1, length(Dinput));
Dffout = zeros(1, length(Dinput));
Dff = zeros(1, length(qyline));
E0ff = zeros(1, length(qyline));

%% FF phase
disp('start FF phase')

mu = mu_ran(1);
% for j = 16: 26
for j = 1: length(J_ran)
    J = J_ran(j);
    Jk = J .* sin(kx) .* sin(ky);
    
    % for qq2 = 1: length(qyline)
    for qq = 1: length(qyline)

        % [Xi,La] = XiLaMaker(kx, ky, qyline(qq2), qyline(qq), t, J, mu);
        [Xi,La] = XiLaMaker(kx, ky, qx0line(qq2), qyline(qq), t, J, mu);

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

%% 실수는 성공의 어머니이다....
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
%% quasiparticle band structure and correlation
ii=17;
[Xi_test,La_test] = XiLaMaker(kx, ky, 0, q_self(ii), t, J_ran(ii), mu);
[~, ~, correlff_test, ~, iu, id] = FFgapsolver(Xi_test, La_test, Dff_self(ii), U);

figure
surf(reshape(sqrt(Xi_test.^2+Dff_self(ii)^2)+La_test,[N,N]))
hold on
surf(reshape(sqrt(Xi_test.^2+Dff_self(ii)^2)-La_test,[N,N]))
zlim([-1,0]);view(2)
figure
surf(reshape(correlff_test,[N,N]));view(2)

for ii = 1: length(J_ran)
    [Xi_test,La_test] = XiLaMaker(kx, ky, 0, q_self(ii), t, J_ran(ii), mu);
    Area_ff(ii) = sum((sqrt(Xi_test .^2 + Dff_self(ii) ^2) + La_test < 0) * 1) + sum((sqrt(Xi_test .^2 + Dff_self(ii) ^2) - La_test < 0) * 1);
end



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
E_0 = mean(E_0_k, 'omitnan');

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