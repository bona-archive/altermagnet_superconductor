%%
clear all
load('C:\Users\hsb83\Desktop\ff_1000_V0to1.5.mat')
load('C:\Users\hsb83\Desktop\s_1000_V0to1.5.mat')
load('C:\Users\hsb83\Desktop\p_1000_V0to1.5\E0_self, D_self, pesd1000_V0to1.5\Dp_self.mat')

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
clear all
N=1000;
kx = linspace(-pi,pi, N+1);
kx(end)=[];
t = 0.5;
mu = -1;
J_ran = linspace(0,1,50);
jjs = 11; % J/2t = 0.2
jjff = 21; % J/2t = 0.4
jjp = 40; % J/2t = 0.8

jj = jjff; % 여기서 s, p+ip, ff 결정

J = J_ran(jj);
if jj == jjs
    T_ran = linspace(0, 0.02, 100);
elseif jj == jjff
    T_ran = linspace(0, 0.01, 50);
    % T_ran_start = linspace(0, 0.005, 200);
    % T_ran_end = linspace(0.005, 0.01, 10);T_ran_end(1) = [];
    % T_ran = [T_ran_start T_ran_end];
elseif jj == jjp
    T_ran = linspace(0, 0.004, 100);
end

q_self = [0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.0565486677646163	0.201061929829747	0.263893782901543	0.295309709437441	0.333008821280518	0.358141562509236	0.383274303737955	0.402123859659494	0.420973415581032	0.439822971502571	0.458672527424110	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];

% 
% for i = 1: length(T_ran)
%     t = T_ran(i);
%     Delta = sqrt(1 - (t/T_ran(end))^2);
%     Delta_asd(i) = sqrt(1 - (t/T_ran(end))^2);
% 
%     Ek = sqrt(ek.^2 + Delta.^2);
%     en = Ek(:);
%     K_0T(i) = 1 + 2 * sum(1/t .* exp(1 / t * en) ./ ((exp(1 / t * en) + 1) .^ 2 ));
% 
% end




if jj == jjs
    spinsymm = 1; % 1 = singlet,-1=triplet
    % Delta = Ds_self(jj);
    Delta = 0.2305;
    qqx = 0;
elseif jj == jjff
    spinsymm = 1;
    % Delta = Dff_self(jj);
    Delta = 0.0261;
    qqx = q_self(jjff);
elseif jj == jjp
    spinsymm = -1;
    % Delta = Dp_self(23, jj);
    Delta = 0.0120;
    qqx = 0;
end
spinsymm
Delta
qqx
    %% 
    

A = 0.001;
for j = 1: length(T_ran)
    T = T_ran(j);
    j
    % Delta = sqrt(1 - (T/T_ran(end)));
    % Delta = tanh(T_ran(end) * sqrt((T_ran(end)./T_ran(j)-1)));
    % DeltaT(j) = Delta;
    J0(j)=0;
    for i = 1: length(kx)
        for l = 1: length(kx)
        kkx = kx(i);
        kky = kx(l);
        % ek = -2*t*cos(kkx)-mu;
        % ek = -2*t*cos(kkx) - 2*t*cos(kky)-mu;
        ekq = -2*t*(kkx + qqx)^2 - 2*t*(kky)^2 - mu;
        emkq = -2*t*(-kkx + qqx)^2 - 2*t*(-kky)^2 - mu;
        Jkq = J * (kkx + qqx) * (kky);
        Jmkq = J * (-kkx + qqx) * (-kky);
        if jj == jjs
        Deltak = Delta * (1);
        elseif jj == jjff
        Deltak = Delta * (1);
        elseif jj == jjp
        Deltak = Delta * (kkx + 1i * kky);
        end
        Hbdg = [ekq-spinsymm*Jkq Deltak;conj(Deltak) -emkq-Jmkq];
        [u,d] = eig(Hbdg);
        J0(j) = J0(j) + (kky) *(FD(-abs(d(1,1))+A*(kky), T)-FD(-abs(d(2,2))-A*(kky),T));
        en1(i, l)=d(1,1);
        en2(i, l)=d(2,2);
        ekqq(i,l) = ekq+Jkq;
        ekqmq(i,l) = emkq-Jmkq;
        end
    end
    if j == 1
        figure;surf(en1,'EdgeColor','none');hold on;surf(zeros(N,N));surf(en2,'EdgeColor','none');zlim([-0.1,0.1]);view(2)
    end
end
K = abs(abs(J0(end))+J0);
lam = 1./sqrt(abs(abs(J0(end))+J0));
dT = T_ran(2)-T_ran(1);

% figure;plot(DeltaT);title('\Delta(T)')
figure;plot(J0);title('K_1(T)')
figure;plot(K/K(1));title('K(T)')
% figure;hold on;surf(en1);surf(en2);title('E(\Delta=0)')
figure;plot(lam);title('\lambda')
figure;plot(T_ran(1:end-1)/T_ran(end),diff(lam)/dT);title('d\lambda/dT');ylim([-0.1,5]);
% figure
% plot(log(lam));title('log(\lambda)');hold on
%% 
plot(T_ran/T(end),lam/lam(1))
hold on
xlabel('T/T_c')
ylabel('\lambda(T)/\lambda(0)')
%% 
plot(exp(T_ran/T(end)),(lam/lam(1)))
hold on
xlabel('exp(T/T_c)')
ylabel('\lambda(T)/\lambda(0)')
%% 

% figure
plot(T_ran/T(end),K/K(1))
hold on
xlabel('T/T_c')
ylabel('K(0,T)/K(0,0)')
%% 

A = 0.001;
ky=kx;
[kx,ky]=meshgrid(kx,ky);
%%
Jkq = J * (kx+qqx) .* (ky);
ekq = -2*t*((kx+qqx).^2+(ky).^2)-mu;
Jkqm = J * (-kx+qqx) .* (-ky);
ekqm = -2*t*((-kx+qqx).^2+(ky).^2)-mu;
Xikq=(1/2)*(ekq+Jkq+ekqm-Jkqm);
Lakq=(1/2)*(ekq+Jkq-ekqm+Jkqm);
Ekq = sqrt(Xikq.^2+Delta^2);
Ekq1 = Ekq + Lakq;
Ekq2 = Ekq - Lakq;

for j = 1: length(T_ran)
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
        
        J0(j) = J0(j) + (kky) *(FD(-abs(Ekq1(i,l))+A*(kky), T)-FD(-abs(Ekq2(i,l))-A*(kky),T));
        
        end
    end
end
K = abs(abs(J0(end))+J0);
lam = 1./sqrt(abs(abs(J0(end))+J0));
dT = T_ran(2)-T_ran(1);

% figure;plot(DeltaT);title('\Delta(T)')
figure;plot(J0);title('K_1(T)')
figure;plot(K/K(1));title('K(T)')
% figure;hold on;surf(en1);surf(en2);title('E(\Delta=0)')
figure;plot(lam);title('\lambda')
figure;plot(T_ran(1:end-1)/T_ran(end),diff(lam)/dT);title('d\lambda/dT');ylim([-0.1,5]);
% figure
% plot(log(lam));title('log(\lambda)');hold on

%% 


function Fd = FD(en, T)
 Fd = 1/(exp(en/T) + 1);
end

function F_pr = FD_prime(en, T)
    F_pr = -(1/T) * (exp(en/T)) * (exp(en/T)+1)^2; 
end