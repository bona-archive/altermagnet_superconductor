%% parameters
% clc
% close all
clear all

N = 50;
t = 0.5;
mu = -0.5;
J = 1;
ParamN = 30;


% U_ran = linspace(0, 2, ParamN+1);
U_ran = ones(1,1) * 1.5;


V_ran = linspace(0, 2, ParamN);


J_ran = linspace(0, 1, ParamN);
% J_ran = ones(1,1) * J;


% mu_ran = linspace(mu, -mu, ParamN);
% mu_ran = ones(1,1) * mu;


Dend = 0.5;
Dinput=linspace(0, Dend, Dend * 1000);


kx = linspace(-pi, pi, N+1);
kx(end) = [];
ky = linspace(-pi, pi, N+1);
ky(end) = [];
[kx, ky] = meshgrid(kx, ky);
kx = reshape(kx, [numel(kx), 1]);
ky = reshape(ky, [numel(ky), 1]);


gesk = cos(kx) + cos(ky);
geska = abs(gesk).^2;


gdk = cos(kx) - cos(ky);
gdka = abs(gdk).^2;


epsilon=1e-5;
ek = e_k(kx, ky, 't', t);
%% s-wave

disp('start s-wave')

figure
hold on
plot(Dinput,Dinput,'k')
for u = 1: length(U_ran)


    U = U_ran(u);


    for j = 1: length(J_ran)


        J = J_ran(j);
        Jk = J_k(kx, ky, 'J', J);
        % mu = mu_ran(1);

        % for j = 1: length(mu_ran)
        % 
        % 
        %     mu = mu_ran(j);
        % 
        % 
        %     if length(mu_ran) == 1
        %         clear j
        %         j = jj;
        %     else
        %     end


            for i = 1: length(Dinput)


                Delta = Dinput(i);
                % Delta_es = Dinput(i);
                % Delta_d = Dinput(i);


                %%% analytic


                Ek = sqrt((ek - mu).^2 + Delta^2);
                Ek_u = sqrt((ek - mu).^2 + Delta^2) + Jk; %% s-wave
                Ek_d = sqrt((ek - mu).^2 + Delta^2) - Jk;
                ind_u = Ek_u > 0;
                ind_d = Ek_d > 0;
                Delta_k = U * Delta / 2 ./ Ek .* (ind_u - (~ind_d));
                Dsout(i) = mean(Delta_k);
                E_0_s_k = (ek - mu) - Ek - Delta^2 ./ Ek / 2 .* (ind_u - (~ind_d))...
                    + Ek_u .* (~ind_u) + Ek_d .* (~ind_d);
                condsout(i)=mean(E_0_s_k);
                % Eesk = sqrt((ek - mu).^2 + Delta^2 * geska); %% es-wave
                % Eesk_u = sqrt((ek - mu).^2 + Delta^2 * geska) + Jk;
                % Eesk_d = sqrt((ek - mu).^2 + Delta^2 * geska) - Jk;
                % indes_u = Eesk_u > 0;
                % indes_d = Eesk_d > 0;
                % Deltaes_k = U * Delta * geska / 2 ./ Eesk .* (indes_u - (~indes_d));
                % mean(Deltaes_k);
                % Desout(i) = mean(Deltaes_k);
                
                % Edk = sqrt((ek - mu).^2 + Delta^2 * gdka); %% d-wave
                % Edk_u = sqrt((ek - mu).^2 + Delta^2 * gdka) + Jk;
                % Edk_d = sqrt((ek - mu).^2 + Delta^2 * gdka) - Jk;
                % indd_u = Edk_u > 0;
                % indd_d = Edk_d > 0;
                % Deltad_k = U * Delta * gdka / 2 ./ Edk .* (indd_u - (~indd_d));
                % mean(Deltad_k);
                % Ddout(i) = mean(Deltad_k);
                
                Esout(i)=-sum(Ek_u(Ek_u>0))-sum(Ek_d(Ek_d>0));

            end
            plot(Dinput,Dsout);


            [~ , I] = min(abs(Dsout ./ (Dinput+epsilon) - 1));
            if (Dsout(2) - 0) / (Dinput(2) - Dinput(1)) < 1
                I = 1;
                Ds=0;
            else
                Ds=Dsout(I);
            end
            
            I_arc(u,j)=I;
            Dss(u,j) = Ds;
            Ek = sqrt((ek - mu).^2 + Ds^2);
            Ek_u = sqrt((ek - mu).^2 + Ds^2) + Jk;
            Ek_d = sqrt((ek - mu).^2 + Ds^2) - Jk;
            ind_u = Ek_u > 0;
            ind_d = Ek_d > 0;
            E_0_s_k = (ek - mu) - Ek - Ds^2 ./ Ek / 2 .* (ind_u - (~ind_d))...
                + Ek_u .* (~ind_u) + Ek_d .* (~ind_d);   
            
            conds(u,j) = mean(E_0_s_k);
            E_0_ss(u,j) = -sum(Ek_u(Ek_u>0)) - sum(Ek_d(Ek_d>0)); %% particle hole symmetry for Bogoliubon is used here
            BFSfinder(u,j) = min(Ek_d);


            % [~ , I] = min(abs(Desout ./ (Dinput+epsilon) - 1)); % es wave
            % if (Desout(2) - 0) / (Dinput(2) - Dinput(1)) < 1
            %     I = 1;
            %     Des=0;
            % else
            %     Des=Desout(I);
            % end
            % Eesk = sqrt((ek - mu).^2 + Des^2 * geska);
            % Eesk_u = sqrt((ek - mu).^2 + Des^2 * geska) + Jk;
            % Eesk_d = sqrt((ek - mu).^2 + Des^2 * geska) - Jk;
            % indes_u = Eesk_u > 0;
            % indes_d = Eesk_d > 0;
            % E_0_es_k = (ek - mu) - Eesk + Des^2 * geska ./ Eesk / 2 .* (indes_u - (~indes_d))...
            %     + Eesk_u .* (~indes_u) + Eesk_d .* (~indes_d);
            % E_0_eses(u,j) = mean(E_0_es_k);
            %
            % [~ , I] = min(abs(Ddout ./ (Dinput+epsilon) - 1)); % d-wave
            % if (Ddout(2) - 0) / (Dinput(2) - Dinput(1)) < 1
            %     I = 1;
            %     Dd=0;
            % else
            %     Dd=Ddout(I);
            % end
            % Edk = sqrt((ek - mu).^2 + Dd^2 * gdka);
            % Edk_u = sqrt((ek - mu).^2 + Dd^2 * gdka) + Jk;
            % Edk_d = sqrt((ek - mu).^2 + Dd^2 * gdka) - Jk;
            % indd_u = Edk_u > 0;
            % indd_d = Edk_d > 0;
            % E_0_d_k = (ek - mu) - Edk + Ds^2 * gdka ./ Edk / 2 .* (indd_u - (~indd_d))...
            %     + Edk_u .* (~indd_u) + Edk_d .* (~indd_d);
            % E_0_dd(u,j) = mean(E_0_d_k);


        % end
    end
end


% E_0_s=zeros(length(E_0_ss),length(E_0_ss));


for i=1:length(V_ran)
% E_0_s(i,:)=E_0_ss;
BFSfinders(i,:)=BFSfinder;
condEs(i,:)=conds;
end


% figure;surf(E_0_s);xlabel('J')


%% FFLO ver1 (select q, scan Delta)


% 
% qx_ran = linspace(-pi, pi, N+1);qx_ran(end)=[];
% qy_ran = qx_ran;
% 
% qx = linspace(-pi, pi, N+1);qx(end)=[];
% qy = qx;
% [qx, qy] = meshgrid(qx, qy);
% qx = reshape(qx, [numel(qx), 1]);
% qy = reshape(qy, [numel(qy), 1]);
% qx0line = qx(N * (N / 2) + 1 : N * (N / 2 + 1));
% qyline = qy(N * (N / 2) + 1 : N * (N / 2 + 1));
% 
% for u = 1: length(U_ran)
% 
% 
%     U = U_ran(u);
% 
%     qq = N / 2 ; %select q
%     kx1 = kx + qx0line(qq);
%     kx2 = kx - qx0line(qq);
%     ky1 = ky + qyline(qq);
%     ky2 = ky - qyline(qq);
% 
%     for j = 1: length(J_ran)
% 
% 
%         J = J_ran(j);
%         Jk1 = J * sin(kx1) .* sin(ky1);                     % J(k+q)
%         Jk2 = J * sin(kx2) .* sin(ky2);                     % J(k-q)
% 
% 
%         xik1 = - 2 * t * cos(kx1) - 2 * t * cos(ky1) - mu;  % xi(k+q)
%         xik2 = - 2 * t * cos(kx2) - 2 * t * cos(ky2) - mu;  % xi(k-q)
% 
% 
%         Xi = .5 * (xik1 + Jk1 + xik2 - Jk2);                % Xi(k;q)
%         La = .5 * (xik1 + Jk1 - xik2 + Jk2);                % Lambda(k;q)
% 
% 
%             for d = 1: length(Dinput) %% scan Delta
%                 Delta = Dinput(d);
%                 [~, D_ffout] = FFgapsolver(Xi, La, Delta, U);
%                 Dffout(d) = D_ffout;
%             end
% 
% 
%             [~ , I] = min(abs(Dffout ./ (Dinput+epsilon) - 1));
% 
% 
%             if (Dffout(2) - 0) / (Dinput(2) - Dinput(1)) < 1
% 
% 
%                 I = 1;
%                 Df=0;
% 
% 
%             else
% 
% 
%                 Df=Dffout(I);
% 
% 
%             end
% 
% 
%             [Eff, Dffout] = FFgapsolver(Xi, La, Df, U);
%             self_E0ff(u, j)=Eff;
% 
% 
% 
%     end
% end
% 
% for i=1:length(V_ran)
% E_0_ff(i,:)=self_E0ff;
% end

% figure
% hold on
% surf(J_ran,U_ran,self_E0ff-E_0_ss);
% [C,H]=contour(J_ran,U_ran,self_E0ff-E_0_ss,[0,0],'k',Linewidth=1.5);
% xlabel('J')
% ylabel('U')
%% FFLO ver2 (select input Delta, scan q)
disp('start FF')

% N = 30;
% t = 0.5;
% mu = -0.5;
% ParamN = 20;
% 
% U_ran=linspace(0, 3.5, ParamN);
% U_ran=ones(1,1) * 2.5;
% J_ran = linspace(0, 1, ParamN);
% J_ran = ones(1,1) * 1;


qx = linspace(-pi, pi, N+1);qx(end)=[];
qy = qx;
[qx, qy] = meshgrid(qx, qy);
qx = reshape(qx, [numel(qx), 1]);
qy = reshape(qy, [numel(qy), 1]);
qx0line = qx(N * (N / 2) + 1 : N * (N / 2 + 1));
qyline = qy(N * (N / 2) + 1 : N * (N / 2 + 1));
kx0line = kx(N * (N / 2) + 1 : N * (N / 2 + 1));
kyline = ky(N * (N / 2) + 1 : N * (N / 2 + 1));
% Dend = 10;
% Dinput=linspace(0, Dend, Dend * 50);
% 
% kx = linspace(-pi, pi, N+1);
% kx(end) = [];
% ky = linspace(-pi, pi, N+1);
% ky(end) = [];
% [kx, ky] = meshgrid(kx, ky);
% kx = reshape(kx, [numel(kx), 1]);
% ky = reshape(ky, [numel(ky), 1]);


gsk = ones(size(kx));


gesk = cos(kx) + cos(ky);
geska = abs(gesk).^2;


gdk = cos(kx) - cos(ky);
gdka = abs(gdk).^2;


epsilon=1e-5;
ek = e_k(kx, ky, 't', t);


figure
hold on
plot(Dinput,Dinput,'k')

for u = 1: length(U_ran)


    U = U_ran(u);


    for j = 1: length(J_ran) 


        J = J_ran(j);
        mu = mu_ran(1);

        % for j = 1: length(mu_ran)
        % 
        % 
        %     mu = mu_ran(j)
        % 
        % 
        %     if length(mu_ran) == 1
        %         clear j
        %         j = jj;
        %     else
        %     end
        % 
        %     j

        for d = 1: length(Dinput) %% select input Delta


            Delta = Dinput(d);


            for qq = 1: length(qx0line) %% scan q


                kx1 = kx + qx0line(qq);
                kx2 = kx - qx0line(qq);
                ky1 = ky + qyline(qq);
                ky2 = ky - qyline(qq);


                Jk1 = J * sin(kx1) .* sin(ky1);                     % J(k+q)
                Jk2 = J * sin(kx2) .* sin(ky2);                     % J(k-q)


                xik1 = - 2 * t * cos(kx1) - 2 * t * cos(ky1) - mu;  % xi(k+q)
                xik2 = - 2 * t * cos(kx2) - 2 * t * cos(ky2) - mu;  % xi(k-q)


                Xi = .5 * (xik1 + Jk1 + xik2 - Jk2);                % Xi(k;q)
                La = .5 * (xik1 + Jk1 - xik2 + Jk2);                % Lambda(k;q)




                [E_0ff, ~, cond0fff] = FFgapsolver(Xi, La, Delta, U);
                E0ff(qq) = E_0ff;
                cond0ff(qq) = cond0fff;


            end

            
            % q_index_minE = find(E0ff == min(E0ff));  %% "find q" with minimum occupied Energy
            % qx_minE = qx0line(q_index_minE(1));
            % qy_minE = qyline(q_index_minE(1));
            % kx1 = kx + qx_minE;
            % kx2 = kx - qx_minE;
            % ky1 = ky + qy_minE;
            % ky2 = ky - qy_minE;
            % Jk1 = J * sin(kx1) .* sin(ky1);                     % J(k+q)
            % Jk2 = J * sin(kx2) .* sin(ky2);                     % J(k-q)
            % xik1 = - 2 * t * cos(kx1) - 2 * t * cos(ky1) - mu;  % xi(k+q)
            % xik2 = - 2 * t * cos(kx2) - 2 * t * cos(ky2) - mu;  % xi(k-q)
            % Xi = .5 * (xik1 + Jk1 + xik2 - Jk2);                % Xi(k;q)
            % La = .5 * (xik1 + Jk1 - xik2 + Jk2);                % Lambda(k;q)
            % [E_0ff, D_ffout, condff] = FFgapsolver(Xi, La, Delta, U);
            % Effout(d) = E_0ff;
            % Dffout(d) = D_ffout;
            % qyminE(d) = qy_minE;
            % condffout(d) = condff;


            q_index_minE = find(cond0ff == min(cond0ff));  %% "find q" with minimum condensation Energy
            qx_minEc = qx0line(q_index_minE(1));
            qy_minEc = qyline(q_index_minE(1));
            kx1 = kx + qx_minEc;
            kx2 = - kx + qx_minEc;
            ky1 = ky + qy_minEc;
            ky2 = - ky + qy_minEc;
            Jk1 = J * sin(kx1) .* sin(ky1);                     % J(k+q)
            Jk2 = J * sin(kx2) .* sin(ky2);                     % J(k-q)
            xik1 = - 2 * t * cos(kx1) - 2 * t * cos(ky1) - mu;  % xi(k+q)
            xik2 = - 2 * t * cos(kx2) - 2 * t * cos(ky2) - mu;  % xi(k-q)
            Xi = .5 * (xik1 + Jk1 + xik2 - Jk2);                % Xi(k;q)
            La = .5 * (xik1 + Jk1 - xik2 + Jk2);                % Lambda(k;q)
            [E_0ffc, D_ffoutc, condffc] = FFgapsolver(Xi, La, Delta, U);
            Effoutc(d) = E_0ffc;
            Dffoutc(d) = D_ffoutc;
            qyminEc(d) = qy_minEc;
            condffoutc(d) = condffc;


        end


        % [~ , I] = min(abs(Dffout ./ (Dinput+epsilon) - 1)); % solve gapeq for q from the minimum occupied energy
        % if (Dffout(2) - 0) / (Dinput(2) - Dinput(1)) < 1
        %     I = 1;
        %     Dff=0;
        %     E_0_fff(u,j) = Effout(1);
        %     condfff(u,j) = condffout(1);
        %     qymin(u,j) = qymin(I);
        %     Dffff(u,j) = 0;
        % else
        %     Dff=Dffout(I);
        %     E_0_fff(u,j) = Effout(I);
        %     condfff(u,j) = condffout(I);
        %     qymin(u,j) = qyminE(I);
        %     Dffff(u,j) = Dffout(I);
        % end

        [~ , I] = min(abs(Dffoutc ./ (Dinput+epsilon) - 1)); % solve gapeq for q from the minimum condensation energy
        if (Dffoutc(2) - 0) / (Dinput(2) - Dinput(1)) < 1
            I = 1;
            Dffc=0;
            E_0_fffc(u,j) = Effoutc(1);
            condfffc(u,j) = condffoutc(1);
            qyminc(u,j) = qyminc(I);
            Dffffc(u,j)=0;
        else
            Dffc=Dffoutc(I);
            E_0_fffc(u,j) = Effoutc(I);
            condfffc(u,j) = condffoutc(I);
            qyminc(u,j) = qyminEc(I);
            Dffffc(u,j)=Dffoutc(I);
        end
        kx1 = kx + qyminEc(I);
        kx2 = kx - qyminEc(I);
        ky1 = ky + qyminEc(I);
        ky2 = ky - qyminEc(I);
        Jk1 = J * sin(kx1) .* sin(ky1);                     % J(k+q)
        Jk2 = J * sin(kx2) .* sin(ky2);                     % J(k-q)
        xik1 = - 2 * t * cos(kx1) - 2 * t * cos(ky1) - mu;  % xi(k+q)
        xik2 = - 2 * t * cos(kx2) - 2 * t * cos(ky2) - mu;  % xi(k-q)
        Xi = .5 * (xik1 + Jk1 + xik2 - Jk2);                % Xi(k;q)
        La = .5 * (xik1 + Jk1 - xik2 + Jk2);                % Lambda(k;q)
        Ek = sqrt(Xi.^2 + (Delta * 1).^2);             % E_eta(k;q)
        EEk_uff(j,:) = Ek + La;                                 % E_u^eta(k;q)
        EEk_dff(j,:)= Ek - La;                                 % E_d^eta(k;q)
        
        [E_0ff, D_ffout, condf] = FFgapsolver(Xi, La, Dffc, U);
        
        condff(u,j) = condf;
        BFSfinderfff(u,j) = min(sqrt(Xi.^2 + (Dffc * 1).^2)-La);

        % plot(Dinput,Dffoutc)
    end
end


for i=1:length(V_ran)
% E_0_ff(i,:)=E_0_fff;
% BFSfinderff(i,:)=BFSfinderfff;
% condEff(i,:)=condfff;
% qyminn(i,:)=qymin;

E_0_ffc(i,:)=E_0_fffc;
condEff(i,:)=condff;
qyminnc(i,:)=qyminc;
end


% figure;surf(E_0_ff);xlabel('J')
%% p, es, d-wave without mixing
% clc
% close all


disp('start p-wave')


gpk = sin(kx) + 1i * sin(ky);
gpka = abs(gpk).^2;
gesk = cos(kx) + cos(ky);
geska = abs(gesk).^2;
gdk = cos(kx) - cos(ky);
gdka = abs(gdk).^2;


epsilon=1e-5;
ek = e_k(kx, ky, 't', t);


for u = 1: length(V_ran)


    V = V_ran(u);


    for j = 1: length(J_ran)


        jj = j;
        J = J_ran(j);
        Jk = J_k(kx, ky, 'J', J);
        % mu = mu_ran(1);
        
        % for j = 1: length(mu_ran)
        % 
        % 
        %     mu = mu_ran(j);
        % 
        % 
        %     if length(J_ran) == 1
        %         clear j
        %         j = jj;
        %     else
        %     end


        for i = 1: length(Dinput)


            Delta_u = Dinput(i);
            Delta_d = Dinput(i);
            Delta = Dinput(i);


            %%% analytic


            Ek_u = sqrt((ek - mu + Jk).^2 + Delta_u^2 * gpka);
            Ek_d = sqrt((ek - mu - Jk).^2 + Delta_d^2 * gpka);
            Delta_u_k = V * Delta_u * gpka / 2 ./ Ek_u;
            Delta_d_k = V * Delta_d * gpka / 2 ./ Ek_d;
            Dpout_u(i) = mean(Delta_u_k);
            Dpout_d(i) = mean(Delta_d_k);


            Eesk = sqrt((ek - mu).^2 + Delta^2 * geska); %% es-wave
            Eesk_u = sqrt((ek - mu).^2 + Delta^2 * geska) + Jk;
            Eesk_d = sqrt((ek - mu).^2 + Delta^2 * geska) - Jk;
            indes_u = Eesk_u > 0;
            indes_d = Eesk_d > 0;
            Deltaes_k = V * Delta * geska / 2 ./ Eesk .* (indes_u - (~indes_d));
            Desout(i) = mean(Deltaes_k);


            Edk = sqrt((ek - mu).^2 + Delta^2 * gdka); %% d-wave
            Edk_u = sqrt((ek - mu).^2 + Delta^2 * gdka) + Jk;
            Edk_d = sqrt((ek - mu).^2 + Delta^2 * gdka) - Jk;
            indd_u = Edk_u > 0;
            indd_d = Edk_d > 0;
            Deltad_k = V * Delta * gdka / 2 ./ Edk .* (indd_u - (~indd_d));
            Ddout(i) = mean(Deltad_k);      


        end


        [~, I] = min(abs(Dpout_u ./ (Dinput + epsilon) - 1)); % p+ip-wave
        if (Dpout_u(2) - 0) / (Dinput(2) - Dinput(1)) < 1
            I = 1;
            Dpu=0;
        else
            Dpu = Dpout_u(I);
        end


        [~,I] = min(abs(Dpout_d./(Dinput+epsilon)-1));
        if (Dpout_d(2)-0) / (Dinput(2)-Dinput(1)) < 1
            I = 1;
            Dpd=0;
        else
            Dpd = Dpout_d(I);
        end
        

        Ek_u = sqrt((ek - mu + Jk).^2 + Dpu^2 * gpka);
        Ek_d = sqrt((ek - mu - Jk).^2 + Dpd^2 * gpka);
        E_0_p_k = (ek - mu) - Ek_u / 2 - Ek_d / 2 ...
            + Dpu^2 * gpka ./ Ek_u / 4 + Dpd^2 * gpka ./ Ek_d / 4;
        % E_0_p_k = (ek - mu) - Ek_u / 2 - Ek_d / 2;        
        condEp(u,j) = mean(E_0_p_k);
        E_0_p(u,j) = -sum(Ek_u(Ek_u>0))-sum(Ek_d(Ek_d>0)); %% particle hole symmetry for Bogoliubon is used here








        [~ , I] = min(abs(Desout ./ (Dinput+epsilon) - 1)); % es wave
        if (Desout(2) - 0) / (Dinput(2) - Dinput(1)) < 1
            I = 1;
            Des=0;
        else
            Des=Desout(I);
        end

        
        Eesk = sqrt((ek - mu).^2 + Des^2 * geska);
        Eesk_u = sqrt((ek - mu).^2 + Des^2 * geska) + Jk;
        Eesk_d = sqrt((ek - mu).^2 + Des^2 * geska) - Jk;
        indes_u = Eesk_u > 0;
        indes_d = Eesk_d > 0;
        E_0_es_k = (ek - mu) - Eesk + Des^2 * geska ./ Eesk / 2 .* (indes_u - (~indes_d))...
            + Eesk_u .* (~indes_u) + Eesk_d .* (~indes_d);
        condEes(u,j) = mean(E_0_es_k);
        E_0_es(u,j) = -sum(Eesk_u(Eesk_u>0))-sum(Eesk_d(Eesk_d>0)); %% particle hole symmetry for Bogoliubon is used here







        
        [~ , I] = min(abs(Ddout ./ (Dinput+epsilon) - 1)); % d-wave
        if (Ddout(2) - 0) / (Dinput(2) - Dinput(1)) < 1
            I = 1;
            Dd=0;
        else
            Dd=Ddout(I);
        end


        Edk = sqrt((ek - mu).^2 + Dd^2 * gdka);
        Edk_u = sqrt((ek - mu).^2 + Dd^2 * gdka) + Jk;
        Edk_d = sqrt((ek - mu).^2 + Dd^2 * gdka) - Jk;
        indd_u = Edk_u > 0;
        indd_d = Edk_d > 0;
        E_0_d_k = (ek - mu) - Edk + Dd^2 * gdka ./ Edk / 2 .* (indd_u - (~indd_d))...
            + Edk_u .* (~indd_u) + Edk_d .* (~indd_d);
        condEd(u,j) = mean(E_0_d_k);   
        E_0_d(u,j) = -sum(Edk_u(Edk_u>0))-sum(Edk_d(Edk_d>0)); %% particle hole symmetry for Bogoliubon is used here
        % end
    end
end


% figure;surf(E_0_p);title('p-wave');xlabel('J');ylabel('V');
% figure;surf(E_0_es);title('es-wave');xlabel('J');ylabel('V');
% figure;surf(E_0_d);title('d-wave');xlabel('J');ylabel('V');

%% occupied energy

% figure
% hold on
% xlabel('J')
% ylabel('V')
% title( 'Occupied E_{s}-E_{p}')
% surf(J_ran, V_ran, E_0_s / N^2 - E_0_p / N^2)
% colorbar
% [Cop, hop]=contour(J_ran, V_ran, E_0_s / N^2 - E_0_p / N^2, [0, 0], 'k', Linewidth=1.5);
% axis equal
% axis square
% 
% figure
% hold on
% xlabel('J')
% ylabel('V')
% title( 'Occupied E_{s}-E_{eff}')
% surf(J_ran, V_ran, E_0_s / N^2 - E_0_ff / N^2)
% colorbar
% [Coes, hoes]=contour(J_ran, V_ran, E_0_s / N^2 - E_0_ff / N^2, [0, 0], 'k', Linewidth=1.5);
% axis equal
% axis square

% figure
% hold on
% xlabel('J')
% ylabel('V')
% title( 'E_{0s}-E_{0d}')
% surf(mu_ran, V_ran, E_0_s - E_0_d)
% colorbar
% [Cd, hd]=contour(mu_ran, V_ran, E_0_s - E_0_d, [0, 0], 'k', Linewidth=1.5);
% axis equal
% axis square
% 
% figure
% hold on
% xlabel('J')
% ylabel('V')
% title( 'E_{0p}-E_{0es}')
% surf(mu_ran, V_ran, E_0_p - E_0_es)
% colorbar
% [Cpes, hpes]=contour(mu_ran, V_ran, E_0_p - E_0_es, [0, 0], 'k', Linewidth=1.5);
% axis equal
% axis square
% 
% figure
% hold on
% xlabel('J')
% ylabel('V')
% title( 'E_{0p}-E_{0d}')
% surf(mu_ran, V_ran, E_0_p - E_0_d)
% colorbar
% [Cpd, hpd]=contour(mu_ran, V_ran, E_0_p - E_0_d, [0, 0], 'k', Linewidth=1.5);
% axis equal
% axis square


% allMatrices = cat(3, E_0_s, E_0_p, E_0_ff);
% [~, minIndex] = min(allMatrices, [], 3);
% figure
% surf(minIndex)
% xlabel('J')
% ylabel('V')
% view(2)
%% BFS finder

% figure
% hold on
% xlabel('J')
% ylabel('V')
% title( 'BFS in s-wave')
% surf(J_ran, V_ran, BFSfinders / N^2)
% colorbar
% [CBFSs, hBFSs]=contour(J_ran, V_ran, BFSfinders / N^2, [0, 0], 'k', Linewidth=1.5);
% axis equal
% axis square
% 
% figure
% hold on
% xlabel('J')
% ylabel('V')
% title( 'BFS in ff-wave')
% surf(J_ran, V_ran, BFSfinderff / N^2)
% colorbar
% [CBFSff, hBFSff]=contour(J_ran, V_ran, BFSfinderff / N^2, [0, 0], 'k', Linewidth=1.5);
% axis equal
% axis square

%% condensation energy

figure
hold on
xlabel('J')
ylabel('V')
title( 'E_{0s}-E_{0p}')
surf(J_ran, V_ran, condEs - condEp)
colorbar
[Csp, hsp]=contour(J_ran, V_ran, condEs - condEp, [0, 0], 'k', Linewidth=1.5);
axis equal
axis square

figure
hold on
xlabel('J')
ylabel('V')
title( 'E_{0s}-E_{0eff}')
surf(J_ran, V_ran, condEs - condEff)
colorbar
[Csff, hsff]=contour(J_ran, V_ran, condEs-condEff, [1e-10, 1e-10], 'k', Linewidth=1.5);
axis equal
axis square

figure
hold on
xlabel('J')
ylabel('V')
title( 'E_{0s}-E_{0es}')
surf(J_ran, V_ran, condEs - condEes)
colorbar
[Cses, hses]=contour(J_ran, V_ran, condEs - condEes, [0, 0], 'k', Linewidth=1.5);
axis equal
axis square
% figure
% hold on
% xlabel('J')
% ylabel('V')
% title( 'E_{0s}-E_{0d}')
% surf(mu_ran, V_ran, E_0_s - E_0_d)
% colorbar
% [Cd, hd]=contour(mu_ran, V_ran, E_0_s - E_0_d, [0, 0], 'k', Linewidth=1.5);
% axis equal
% axis square
% 
% figure
% hold on
% xlabel('J')
% ylabel('V')
% title( 'E_{0p}-E_{0es}')
% surf(mu_ran, V_ran, E_0_p - E_0_es)
% colorbar
% [Cpes, hpes]=contour(mu_ran, V_ran, E_0_p - E_0_es, [0, 0], 'k', Linewidth=1.5);
% axis equal
% axis square
% 
% figure
% hold on
% xlabel('J')
% ylabel('V')
% title( 'E_{0p}-E_{0d}')
% surf(mu_ran, V_ran, E_0_p - E_0_d)
% colorbar
% [Cpd, hpd]=contour(mu_ran, V_ran, E_0_p - E_0_d, [0, 0], 'k', Linewidth=1.5);
% axis equal
% axis square


% allMatrices = cat(3, E_0_s, E_0_p, E_0_ff);
% [~, minIndex] = min(allMatrices, [], 3);
% figure
% surf(minIndex)
% xlabel('J')
% ylabel('V')
% view(2)

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



function [E_0, Dffout, condff] = FFgapsolver(Xi, La, Delta, U)
Ek = sqrt(Xi.^2 + (Delta * 1).^2);             % E_eta(k;q)
Ek_u = Ek + La;                                 % E_u^eta(k;q)
Ek_d = Ek - La;                                 % E_d^eta(k;q)


ind_u = Ek_u > 0;
ind_d = Ek_d > 0;
Delta_k = U * (Delta * 1) / 2 ./ Ek .* (ind_u - (~ind_d));


Dffout = mean(Delta_k, 'omitnan');               % Delta_eta(q)
% compute ground state energy
E_0_k = Xi - Ek - (Delta * 1).^2 ./ Ek / 2 .* (ind_u - (~ind_d))...
    + Ek_u .* (~ind_u) + Ek_d .* (~ind_d);



condff = mean(E_0_k);                          % E_0^eta(q)
E_0 = - sum(Ek_u(Ek_u>0)) - sum(Ek_d(Ek_d>0)); %% particle hole symmetry for Bogoliubon is used here
end