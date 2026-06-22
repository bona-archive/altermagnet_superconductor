% MATLAB 코드: J 변화에 따른 mu 조정
clear all

% 파라미터 설정
t = 0.5; % 홉핑 파라미터
Nx = 500; % k 공간에서의 격자 포인트 수
Ny = 500;
kx_vals = linspace(-pi, pi, Nx + 1); % kx의 값들
ky_vals = linspace(-pi, pi, Ny + 1); % ky의 값들
kx_vals(end) = [];
ky_vals(end) = [];
qx_vals = linspace(-pi, pi, Nx + 1); % kx의 값들
qy_vals = linspace(-pi, pi, Ny + 1); % ky의 값들
qx_vals(end) = [];
qy_vals(end) = [];

[kx_vals, ky_vals]=meshgrid(kx_vals, ky_vals);

J = 0.4;

tolerance = 0.0001; % 수렴 허용 오차

% 초기 화학 퍼텐셜
mu_initial = -1; % 초기 추정치
target_density = calculate_density(mu_initial, J, t, kx_vals, ky_vals, qx_vals(Nx/2 + 1), qy_vals(Nx/2 + 1));%초기 정상 상태 밀도

% J가 변화할 때 mu를 조정
mu_vals = zeros(Nx, Ny); % 각 J에 대한 mu 저장
mu_vals(Nx/2 + 1, Ny/2 + 1) = mu_initial;
%% 

for idx = ((Nx/2 + 1) + 1): Nx
    for idy = ((Ny/2 + 1) + 1): Ny
    qx = qx_vals(idx);
    qy = qy_vals(idy);
    mu = mu_vals(idx - 1, idy - 1); % 이전 mu 값을 사용하여 새로운 mu 계산 시작
    current_density = calculate_density(mu, J, t, kx_vals, ky_vals, qx, qy);
    
    % 목표 밀도에 도달할 때까지 화학 퍼텐셜을 수정
    while abs(current_density - target_density) > tolerance
        if current_density > target_density
            mu = mu - 0.00001; % 밀도가 높으면 mu를 줄임
        else
            mu = mu + 0.00001; % 밀도가 낮으면 mu를 증가
        end
        current_density = calculate_density(mu, J, t, kx_vals, ky_vals, qx, qy);
    end
    
    mu_vals(idx) = mu; % 조정된 mu 값 저장
    end
end

% 결과 그래프
figure;
plot(J_vals, mu_vals, '-o');
xlabel('J');
ylabel('Chemical potential (\mu)');
title('Chemical Potential Adjustment to Maintain Constant Density');
grid on;
%% 

% 화학 퍼텐셜 조정 함수
function density = calculate_density(mu, J, t, kx_vals, ky_vals, qx, qy)
    density = 0; % 정상 상태 밀도 초기화
    
    
kx1 = kx_vals + qx;
kx2 = -kx_vals + qx;
ky1 = ky_vals + qy;
ky2 = -ky_vals + qy;

Jk1 = J * sin(kx1) .* sin(ky1);
Jk2 = J * sin(kx2) .* sin(ky2);

xik1 = - 2 * t * cos(kx1) - 2 * t * cos(ky1) - mu;
xik2 = - 2 * t * cos(kx2) - 2 * t * cos(ky2) - mu;

Xi = .5 * (xik1 + Jk1 +xik2 - Jk2);                            
La = .5 * (xik1 + Jk1 -xik2 + Jk2);


    energy = xi_k + J_k;
    s = sign(energy);
    density = sum(s(:) == -1);

    density = density / (Nx * Ny); % 평균 밀도
end