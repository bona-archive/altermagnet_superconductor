%% 
clear all

deldel=linspace(0.6,0.7,10);
figure
for qqq=1:length(deldel)

    Del=deldel(qqq);

%% 
t=1;
J=1.5;
mu=-1;
m=1/40;
Vs=-2;
% Del=0.64;
% Vs=-0.1;
% Del=0.1;
str='t=%.2f, J=%.2f, \\mu=%.2f, V_s=%.2f, \\Delta=%.2f';
beta=5000;
NOGx=51; % Number Of Grid x
NOGy=51; % Number Of Grid y
Np=NOGx-1;
NF=1/((NOGx-1)*(NOGy-1)); % Normalization Factor
kkx=linspace(-pi,pi,NOGx);kkx(end)=[];
kky=linspace(-pi,pi,NOGy);kky(end)=[];

% qend=pi;
qqx=kkx(end/2+1:end);
qqy=kky(end/2+1:end);
% qqx=linspace(-qend,qend,NOGx);qqx(end)=[];
% qqy=linspace(-qend,qend,NOGy);qqy(end)=[];
for i=1:length(qqx)
    % for j=1:length(qqy)
        [correl,en1,en2,en3,en4,delta]=Hbdg_solver('s',NF,t,mu,J,beta,Vs,kkx,kky,Del,true,qqx(i),0);
        EE(i)=sum(en1(en1<0))+sum(en2(en2<0))+sum(en3(en3<0))+sum(en4(en4<0));
    % end
end
% mode=floor(mod(EE,10^3));
% E=round(EE-mode+10^3-1,3);
E=fix((EE*10^(-3))*10^(11))/10^11;
%% 

hold on
plot(E-min(E))
title(sprintf(str, t, J, mu,Vs,Del))
[xx, yy]=find(E==min(E));
qqx(yy(length(yy)/2))
%% 

minE=min(min(E));
[minx, miny]=find(E==minE);
qxmin=qqx(minx);
qymin=qqy(miny);
valid_indices = find(qxmin >= -pi/2 & qxmin <= pi/2 & qymin >= -pi/2 & qymin <= pi/2);
end
% scatter를 이용하여 [-pi/2, pi/2] 범위 내에 속하는 점만 찍기
% scatter(qqx(valid_indices), qqy(valid_indices), 'o');
%% 2D q surf
% figure
% surf(qqx,qqy,E)
% xlabel('q_x');ylabel('q_y');
% % xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
% % yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
% xticks([-pi/2 0 pi/2]);xticklabels({'-\pi/2', '0', '\pi/2'});
% yticks([-pi/2 0 pi/2]);yticklabels({'-\pi/2', '0', '\pi/2'});
% %% 2D q case contour and scatter
% figure
% contour(qqx,qqy,EE);hold on;scatter(qxmin,qymin,'r','filled')
% xlabel('q_x');ylabel('q_y');
% % xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
% % yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
% xlim([-pi/2 pi/2]);ylim([-pi/2 pi/2])
% xticks([-pi/2 0 pi/2]);xticklabels({'-\pi/2', '0', '\pi/2'});
% yticks([-pi/2 0 pi/2]);yticklabels({'-\pi/2', '0', '\pi/2'});
% title(sprintf(str, t, J, mu,Vs,Del))
%% 1D q case plot and scatter
% figure
% plot(qqx,E)
% xlabel('q_x');ylabel('q_y');
% hold on
% scatter(qqx(minx),ones(1,4)*minE)
%% Degenerated four points overlap
% 
% phiA=linspace(-pi,pi,10);
% phiB=phiA;phiC=phiA;phiD=phiA;
