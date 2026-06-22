%% parameters

clear all
N = 1000;
kx = linspace(-pi,pi, N+1);
kx(end)=[];
ky=kx;
[kx,ky]=meshgrid(kx,ky);

t = 0.5;
mu = -1;
J_ran = linspace(0,1,50);
jjs = 11; % J/2t = 0.2
jjff = 21; % J/2t = 0.4
jjp = 40; % J/2t = 0.8

jj = jjff; % s, p+ip, ff 결정

J = J_ran(jj);

q_self = [0	 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.0565486677646163	0.201061929829747	0.263893782901543	0.295309709437441	0.333008821280518	0.358141562509236	0.383274303737955	0.402123859659494	0.420973415581032	0.439822971502571	0.458672527424110	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];


if jj == jjs
    % Delta = Ds_self(jj);
    Delta = 0.2305;
    qqx = 0;
    T_ran = linspace(0, 0.02, 100);
elseif jj == jjff
    % Delta = Dff_self(jj);
    Delta = 0.0261;
    qqx = q_self(jjff);
    T_ran = linspace(0, 0.01, 100);
    % T_ran_start = linspace(0, 0.005, 300);
    % T_ran_end = linspace(0.005, 0.01, 11);T_ran_end(1) = [];
    % T_ran = [T_ran_start T_ran_end];
elseif jj == jjp
    % Delta = Dp_self(23, jj);
    Delta = 0.0120;
    qqx = 0;
    T_ran = linspace(0, 0.004, 100);
end
Delta
qqx
%% penetration depth calculation

A = 0.001;
Jkq = J * (kx+qqx) .* (ky);
ekq = -2*t*((kx+qqx).^2+(ky).^2)-mu;

Jkqm = J * (-kx + qqx) .* (-ky);
ekqm = -2*t*((-kx + qqx).^2 + (ky).^2) - mu;

% s-wave spectrum
Eks1 = sqrt(ekq.^2 + Delta^2) + Jkq;
Eks2 = sqrt(ekq.^2 + Delta^2) - Jkq;

% p+ip-wave spectrum
gpka = (kx).^2 + (ky).^2;
gpka = kx.^2;
Ekp1 = sqrt((ekq + Jkq).^2 + Delta^2 * gpka);
Ekp2 = sqrt((ekq - Jkq).^2 + Delta^2 * gpka);

% gapless FF spectrum
Xikq = (1/2) * (ekq + Jkq + ekqm - Jkqm);
Lakq = (1/2) * (ekq + Jkq - ekqm + Jkqm);
Ekq1 = sqrt(Xikq.^2 + Delta^2) + Lakq;
Ekq2 = sqrt(Xikq.^2 + Delta^2) - Lakq;

check_num = 100;
for j = 1: length(T_ran)
% j = 300
    T = T_ran(j);
    j
    % Delta = sqrt(1 - (T/T_ran(end)));
    % Delta = tanh(T_ran(end) * sqrt((T_ran(end)./T_ran(j)-1)));
    % DeltaT(j) = Delta;
    J0(j)=0;
    for i = 1: length(kx)
        for l = 1: length(ky)
        kkx = kx(i);
        kky = ky(l);
        
        if jj == jjff
            J0(j) = J0(j) + (kky) * (FD(-abs(Ekq1(i, l)) + A*(kky), T) - FD(-abs(Ekq2(i, l)) - A * (kky), T));
            if j == check_num
            check(i,l) = (N^3 * kky) *(FD(Ekq1(i,l)+A*(kky), T) - FD(Ekq1(i,l)-A*(kky),T)) + (N^3 * kky) *(FD(Ekq2(i,l)+A*(kky), T) - FD(Ekq2(i,l)-A*(kky),T));
            end

        elseif jj == jjs
            J0(j) = J0(j) + (kky) *(FD(Eks1(i,l)+A*(kky), T) - FD(Eks1(i,l)-A*(kky),T)) + (kky) *(FD(Eks2(i,l)+A*(kky), T) - FD(Eks2(i,l)-A*(kky),T));
            if j == check_num 
            check(i,l) = (N^3 * kky) *(FD(Eks1(i,l)+A*(kky), T) - FD(Eks1(i,l)-A*(kky),T)) + (N^3 * kky) *(FD(Eks2(i,l)+A*(kky), T) - FD(Eks2(i,l)-A*(kky),T));
            end

        elseif jj == jjp
            J0(j) = J0(j) + (kky) *(FD(Ekp1(i,l)+A*(kky), T) - FD(Ekp2(i,l)-A*(kky), T)) + (kky) *(FD(Ekp1(i,l)+A*(kky), T) - FD(Ekp2(i,l)-A*(kky), T));
            if j == check_num
            check(i,l) = (N^3 * kky) *(FD(Ekp1(i,l)+A*(kky), T) - FD(Ekp2(i,l)-A*(kky), T)) + (N^3 * kky) *(FD(Ekp1(i,l)+A*(kky), T) - FD(Ekp2(i,l)-A*(kky), T));
            end
        end
        
        end
    end
    if j == check_num
        figure
        surf(check*10^(-8),'LineStyle','none')
    else
    end
end
K = abs(abs(J0(end))+J0);
lam = 1./sqrt(abs(abs(J0(end))+J0));
dT = T_ran(2)-T_ran(1);
dl = diff(lam)/dT;
dl(dl==Inf)=[];
sizedl = length(dl);
if mod(sizedl,2) == 0
    sizedl = sizedl/2;
elseif mod(sizedl,2) ==1
    sizedl = (sizedl+1)/2;
end

figure;plot(T_ran, J0);title('K_1(T)')
figure;plot(T_ran, K/K(1));title('K(T)')
figure;plot(T_ran, lam/lam(1));title('\lambda')
figure;plot(T_ran(1:end-1)/T_ran(end), diff(lam)/dT);title('d\lambda/dT');ylim([-0.1,5]);

% figure
% hold on
% surf(reshape(-Ekq2,[N,N]),'Linestyle','none')
% surf(zeros(N,N),'Linestyle','none')
% surf(reshape(-Ekq1,[N,N]),'Linestyle','none')
% title('gapless FF spectrum')
% camlight

x = T_ran(1:sizedl)/T_ran(end);
y=dl(1:sizedl);

% figure
% aa = Ekq1 < 0;
% surf(aa*1.*Ekq1,'linestyle','none')
% aa = (Ekq2<0);
% hold on
% surf(aa*1.*Ekq2,'linestyle','none')

% %% fitting
% 
% % fitting1: polynomial
% 
% p = polyfit(x, y, 3);
% poly_fit_y = polyval(p, x);
% 
% a = p(4); % 상수항
% b = p(3); % x의 계수
% c = p(2); % x^2의 계수
% d = p(1); % x^3의 계수
% 
% disp(['a = ', num2str(a)]);
% disp(['b = ', num2str(b)]);
% disp(['c = ', num2str(c)]);
% disp(['d = ', num2str(d)]);
% 
% % R^2
% SS_res = sum((y - poly_fit_y).^2);
% SS_tot = sum((y - mean(y)).^2);
% R_squared = 1 - (SS_res / SS_tot);
% 
% disp(['R^2 for Polynomial Fit: ', num2str(R_squared)]);
% 
% % fitting2: exponential
% 
% model = fittype('a*exp(b*x)');
% fit_result = fit(x', y', model, 'StartPoint', [0, 0.1]); % 초기값 설정
% exp_fit_y = fit_result(x)';
% 
% % check a, b
% disp('Fitting Results:');
% disp(fit_result);
% 
% a = fit_result.a;
% b = fit_result.b;
% 
% disp(['a = ', num2str(a)]);
% disp(['b = ', num2str(b)]);
% 
% % R^2
% SS_res = sum((y - exp_fit_y).^2); % Residual Sum of Squares
% SS_tot = sum((y - mean(y)).^2); % Total Sum of Squares
% R_squared = 1 - (SS_res / SS_tot); % R^2 계산
% 
% % 결과 출력
% disp(['R^2 for Exponential Fit: ', num2str(R_squared)]);
% 
% figure
% hold on
% plot(x, y, 'o-', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'Original Data (Red)');
% plot(x, poly_fit_y, '.', 'Color', 'g', 'MarkerSize', 12, 'DisplayName', 'Polynomial Fit');
% plot(x, exp_fit_y, '.', 'Color', 'm', 'MarkerSize', 12, 'DisplayName', 'Exponential Fit');
% legend('Location', 'best');
% xlabel('T/T_c');
% ylabel('d\lambda/dT');
% title('Original Data and Fits');
%% function


function Fd = FD(en, T)
 Fd = 1/(exp(en/T) + 1);
end

function F_pr = FD_prime(en, T)
    F_pr = -(1/T) * (exp(en/T)) * (exp(en/T)+1)^2; 
end