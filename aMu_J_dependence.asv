% MATLAB 코드: J 변화에 따른 mu 조정
% clear all

% 파라미터 설정
t = 0.5; % 홉핑 파라미터
Nx = 500; % k 공간에서의 격자 포인트 수
Ny = 500;
kx_vals = linspace(-pi, pi, Nx); % kx의 값들
ky_vals = linspace(-pi, pi, Ny); % ky의 값들
[kx_vals, ky_vals]=meshgrid(kx_vals,ky_vals);
J_vals = linspace(0, 1, 50); % J 값을 0에서 0.4까지 변화시킴
tolerance = 0.0001; % 수렴 허용 오차

% 초기 화학 퍼텐셜
mu_initial = -1; % 초기 추정치
initial=13;
[Xi_initial, La_initial] = XiLaMaker(reshape(kx_vals, [1, Ny*Nx]), reshape(ky_vals, [1, Nx*Ny]), J_vals(initial), q_self(initial), t, J, mu);
density_target = calculate_density(Xi_initial, La_initial, Dff_self(initial), Nx, Ny);%초기 정상 상태 밀도

% J가 변화할 때 mu를 조정
mu_vals = zeros(size(J_vals)); % 각 J에 대한 mu 저장
mu_vals(1) = mu_initial;
%% 

for idx = (initial+1):length(J_vals)
    J = J_vals(idx);
    mu = mu_vals(idx - 1); % 이전 mu 값을 사용하여 새로운 mu 계산 시작
    [Xi,La] = XiLaMaker(reshape(kx_vals, [1, Ny*Nx]), reshape(ky_vals, [1, Nx*Ny]), 0, q_self(idx), t, J, mu);
    density_current = calculate_density(Xi, La, Dff_self(idx), Nx, Ny);
    
    % 목표 밀도에 도달할 때까지 화학 퍼텐셜을 수정
    while abs(density_current - density_target) > tolerance
        if density_current > density_target
            mu = mu - 0.001; % 밀도가 높으면 mu를 줄임
        else
            mu = mu + 0.001; % 밀도가 낮으면 mu를 증가
        end
        [Xi,La] = XiLaMaker(reshape(kx_vals, [1, Ny*Nx]), reshape(ky_vals, [1, Nx*Ny]), 0, q_self(idx), t, J, mu);        
        density_current = calculate_density(Xi, La, Dff_self(idx), Nx, Ny);
    end
    
    mu_vals(idx) = mu; % 조정된 mu 값 저장
end

% 결과 그래프
figure;
plot(J_vals, mu_vals, '-o');
xlabel('J');
ylabel('Chemical potential (\mu)');
title('Chemical Potential Adjustment to Maintain Constant Density');
grid on;
%% function

% 화학 퍼텐셜 조정 함수
function density = calculate_density(Xi, La, Delta, Nx, Ny)
    density = 0; % 정상 상태 밀도 초기화

    Ek = hypot(Xi, Delta);
    % Ek_u = Ek + La;
    Ek_d = Ek - La;

    energy = Ek_d;
    s = sign(energy);
    density = sum(s(:) == -1);

    density = density / (Nx * Ny); % 평균 밀도
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