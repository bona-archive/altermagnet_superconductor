clear all

N = 30;
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
Dinput = linspace(0, 0.3, 100);


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

% 데이터 조각 수 정의
num_chunks = 10;
chunk_size = N * N / num_chunks;

% 병렬 연산을 위한 배열 초기화
partial_means = zeros(1, num_chunks);
%% FF phase
disp('start FF phase')

mu = mu_ran(1);
for j = 1: 20  %26
    %check point
    j
% for j = 1: length(J_ran)
    J = J_ran(j);
    Jk = J .* sin(kx) .* sin(ky);
     % 임시 배열 초기화
    temp_E0ffout = zeros(1, length(Dinput));
    temp_Dffout = zeros(1, length(Dinput));
    temp_Nout = zeros(1, length(Dinput));
    temp_Dff = zeros(length(qyline)/3, length(qyline)/3);
    temp_E0ff = zeros(length(qyline)/3, length(qyline)/3);
    temp_N = zeros(length(qyline)/3, length(qyline)/3);
    for qq2 = 1: length(qyline)/3
        qq2
        for qq = 1: length(qyline)/3

            [Xi,La] = XiLaMaker(kx, ky, qyline(qq2), qyline(qq), t, J, mu);
 
            for i = 1: length(Dinput)


                Delta = Dinput(i);

                [EEq, DDq, ~, ~, ind_u, ind_d] = FFgapsolver(Xi, La, Delta, U, num_chunks);
                temp_E0ffout(i) = EEq;
                temp_Dffout(i) = DDq;
                temp_Nout(i) = sum(ind_u*1) + sum(ind_d*1);
            end


            [~ , I] = min(abs(Dffout ./ (Dinput + epsilon) - 1));
            if (Dffout(2) - 0) / (Dinput(2) - Dinput(1)) < 1
                I = 1;
                temp_Dff(qq2, qq) = 0;
                temp_E0ff(qq2, qq) = temp_E0ffout(I);
                % E0ff_arc(j, qq2, qq) = E0ffout(I);
                temp_N(qq2, qq) = temp_Nout(I);
            else
                temp_Dff(qq2, qq) = Dinput(I);
                temp_E0ff(qq2, qq) = temp_E0ffout(I);
                % E0ff_arc(j, qq2, qq) = E0ffout(I);
                temp_N(qq2, qq) = temp_Nout(I);

            end
        end
    end
    % parfor 루프 종료 후, 결과를 메인 배열로 복사
    E0ff_arc(j, :, :) = temp_E0ff;
    N_arc(j,:,:) = temp_N;
    % [aa, ~] = min(E0ff);
    % II = find(E0ff == aa);
    % II = II(1);
    % Dff_self(j) = Dff(II);
    % E0ff_self(j) = aa;
    % q_self(j) = qyline(II);
    % [Xi_self,La_self] = XiLaMaker(kx, ky, 0, qyline(II), t, J, mu);
    % [~, ~, ca, ~, ~, ~] = FFgapsolver(Xi_self, La_self, Dff(II), U);
    % correlff_self(j, :,:) = reshape(ca,[N,N]);
end
% E0ff_self = repmat(E0ff_self, length(V_ran), 1);

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


function [overall_E_0, overall_Dffout, correl, Ek_u, ind_u, ind_d] = FFgapsolver(Xi, La, Delta, U, num_chunks)
% GPU 연산 수행
Xi_gpu = gpuArray(Xi);  % Xi 배열을 GPU로 전송
Delta_gpu = gpuArray(Delta);  % Delta를 GPU 배열로 전환


% Ek = sqrt(Xi .^2 + Delta^2)
Ek_gpu = hypot(Xi_gpu, Delta_gpu);

% Ek 결과를 CPU로 다시 가져오기
Ek = gather(Ek_gpu);  % gather를 사용하여 GPU에서 CPU로 전송

Ek_u = Ek + La;
Ek_d = Ek - La;

ind_u = Ek_u > 0;
ind_d = Ek_d > 0;
Delta_k = U * Delta /2 ./ Ek .* (ind_u - (~ind_d));
E_0_k = Xi - Ek + Delta^2 ./ Ek /2 .* (ind_u -(~ind_d))...
    + Ek_u .* (~ind_u) + Ek_d .* (~ind_d);
correl = Delta ./ Ek /2 .* (ind_u -(~ind_d));

chunk_size = length(E_0_k) / num_chunks;
for chunk = 1: num_chunks
    start_idx = round((chunk-1) * chunk_size) + 1;
    end_idx = round(chunk * chunk_size);
    partial_Dffout(chunk) = mean(Delta_k(start_idx:end_idx),'omitnan');
    partial_E_0(chunk) = mean(E_0_k(start_idx:end_idx),'omitnan');
end
overall_Dffout = mean(partial_Dffout,'omitnan');
overall_E_0 = mean(partial_E_0,'omitnan');
% E_0 = mean(E_0_k, 'omitnan');
% Dffout = mean(Delta_k, 'omitnan');
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