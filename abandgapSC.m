clear all
%% 
t=1;
J=1;
mu=-1;
% Delta=0.1755;
% del1=0.1755;
% del2=0.1755;
% % V=1;
% Delta=0.792;
% del1=0.0658;
% del2=0.0658;
% V=1.5;
Delta=0;
del1=0;
del2=0;
V=0;
str='t=%.2f, J=%.2f, \\mu=%.2f, U=%.2f, \\Delta=%.2f';

% del=0;
% del1=0;
% del2=0;
% V=0;

N=200;
kx=linspace(-pi,pi,N+1);kx(end)=[];
ky=linspace(-pi,pi,N+1);ky(end)=[];

for i=1:length(kx)
    for j=1:length(ky)
        % qx=kx(121);
        qx=0;
        qy=0;
        kkx=kx(i);
        kky=ky(j);
        A=-t*(cos(kkx+qx)+cos(kky+qy))-mu;
        B=J*sin(kkx+qx)*sin(kky+qy);
        AA=-t*(cos(-kkx+qx)+cos(-kky+qy))-mu;
        BB=J*sin(-kkx+qx)*sin(-kky+qy);

        Z=-t*(cos(kkx)+cos(kky))-mu;
        W=J*sin(kkx)*sin(kky);
        
        C=Delta;
        gpp=sin(kkx)+1i*sin(kky);
        gpm=sin(kkx)-1i*sin(kky);
        Dc11=V*del1*gpp;
        Dc22=V*del2*gpp;
        Dh11=V*del1*gpp+0*gpm;
        Dh22=0*gpp+V*del2*gpm;

        Hs=[A+B 0 0 -C;0 A-B C 0;0 C -AA-BB 0;-C 0 0 -AA+BB];
        Hpc=[A+B 0 Dc11 0;0 A-B 0 Dc22;conj(Dc11) 0 -AA-BB 0;0 conj(Dc22) 0 -AA+BB];
        Hph=[A+B 0 Dh11 0;0 A-B 0 Dh22;conj(Dh11) 0 -AA-BB 0;0 conj(Dh22) 0 -AA+BB];
        [u,d]=eig(Hs);
        EEE1(i,j)=d(1,1);EEE2(i,j)=d(2,2);EEE3(i,j)=d(3,3);EEE4(i,j)=d(4,4);
        [upc,dpc]=eig(Hpc);
        EEEc1(i,j)=dpc(1,1);EEEc2(i,j)=dpc(2,2);EEEc3(i,j)=dpc(3,3);EEEc4(i,j)=dpc(4,4);
        [uph,dph]=eig(Hph);
        EEEh1(i,j)=dph(1,1);EEEh2(i,j)=dph(2,2);EEEh3(i,j)=dph(3,3);EEEh4(i,j)=dph(4,4);
        EE1(i,j)=1/2*((-AA-BB+A-B)+sqrt((-AA-BB-A+B)^2+4*C^2));
        EE2(i,j)=1/2*((-AA-BB+A-B)-sqrt((-AA-BB-A+B)^2+4*C^2));
        EE3(i,j)=1/2*((-AA+BB+A+B)+sqrt((-AA+BB-A-B)^2+4*C^2));
        EE4(i,j)=1/2*((-AA+BB+A+B)-sqrt((-AA+BB-A-B)^2+4*C^2));
        % 
        Ep1(i,j)=sqrt((Z+W)^2+Dc11*conj(Dc11));
        Ep2(i,j)=-sqrt((Z+W)^2+Dc11*conj(Dc11));
        Ep3(i,j)=sqrt((Z-W)^2+Dc11*conj(Dc11));
        Ep4(i,j)=-sqrt((Z-W)^2+Dc11*conj(Dc11));
        % Ep1(i,j)=sqrt((A+B)^2+Dc11*conj(Dc11));
        % Ep2(i,j)=-sqrt((A+B)^2+Dc11*conj(Dc11));
        % Ep3(i,j)=sqrt((A-B)^2+Dc11*conj(Dc11));
        % Ep4(i,j)=-sqrt((A-B)^2+Dc11*conj(Dc11));        
    end
end

% figure
% surf(EE1)
% hold on
% surf(EE2)
% surf(EE3)
% surf(EE4)
% 
% figure
% plot(EE1,'k')
% hold on
% plot(EE2,'k')
% plot(EE3,'k')
% plot(EE4,'k')
% 
% figure
% surf(Ep1)
% hold on
% surf(Ep2)
% surf(Ep3)
% surf(Ep4)
% 
% figure
% plot(Ep1,'k')
% hold on
% plot(Ep2,'k')
% plot(Ep3,'k')
% plot(Ep4,'k')

figure
plot([EEE1(101:200,101);EEE1(1,101:200)';flip(diag(EEE1(101:200,101:200)));diag(EEE1(101:-1:2,101:200));EEE1(1:101,1);EEE1(101,200:-1:101)'],'k','LineWidth',1)
hold on
plot([EEE2(101:200,101);EEE2(1,101:200)';flip(diag(EEE2(101:200,101:200)));diag(EEE2(101:-1:2,101:200));EEE2(1:101,1);EEE2(101,200:-1:101)'],'k','LineWidth',1)
plot([EEE3(101:200,101);EEE3(1,101:200)';flip(diag(EEE3(101:200,101:200)));diag(EEE3(101:-1:2,101:200));EEE3(1:101,1);EEE3(101,200:-1:101)'],'k','LineWidth',1)
plot([EEE4(101:200,101);EEE4(1,101:200)';flip(diag(EEE4(101:200,101:200)));diag(EEE4(101:-1:2,101:200));EEE4(1:101,1);EEE4(101,200:-1:101)'],'k','LineWidth',1)
plot(zeros(1,601),'k')
xticks([0 N/2 N 3*N/2 2*N 5*N/2 3*N])
xticklabels({'\Gamma', 'X', 'M','\Gamma','M', 'Y', '\Gamma'})
xlim([0 3*N])
ax1=gca;
ax1.XGrid='on'
ax1.GridAlpha=1;
ax1.FontName='Times New Roman';
ax1.FontSize=12;
% title(sprintf(str, t, J, mu, V, del))
title(sprintf(str, t, J, mu, V, Delta))
% 
% figure
% plot([EE1(N/2+1:N,N/2+1);EE1(1,N/2+1:N)';flip(diag(EE1(N/2+1:N,N/2+1:N)));diag(EE1(N/2+1:-1:2,N/2+1:N));EE1(1:N/2+1,1);EE1(N/2+1,N:-1:N/2+1)'],'k','LineWidth',1)

% figure
% ax1=gca;
% plot([EEEc1(101:200,101);EEEc1(1,101:200)';flip(diag(EEEc1(101:200,101:200)));diag(EEEc1(101:-1:2,101:200));EEEc1(1:101,1);EEEc1(101,200:-1:101)'],'k','LineWidth',0.8)
% hold on
% plot([EEEc2(101:200,101);EEEc2(1,101:200)';flip(diag(EEEc2(101:200,101:200)));diag(EEEc2(101:-1:2,101:200));EEEc2(1:101,1);EEEc2(101,200:-1:101)'],'k','LineWidth',0.8)
% plot([EEEc3(101:200,101);EEEc3(1,101:200)';flip(diag(EEEc3(101:200,101:200)));diag(EEEc3(101:-1:2,101:200));EEEc3(1:101,1);EEEc3(101,200:-1:101)'],'k','LineWidth',0.8)
% plot([EEEc4(101:200,101);EEEc4(1,101:200)';flip(diag(EEEc4(101:200,101:200)));diag(EEEc4(101:-1:2,101:200));EEEc4(1:101,1);EEEc4(101,200:-1:101)'],'k','LineWidth',0.8)
% xticks([0 N/2 N 3*N/2 2*N 5*N/2 3*N])
% xticklabels({'\Gamma', 'X', 'M','\Gamma','M', 'Y', '\Gamma'})
% xlim([0 3*N])
% ax1.XGrid='on'
% ax1.GridAlpha=1;
% ax1.FontName='Times New Roman';
% ax1.FontSize=12;

% figure
% ax1=gca;
% plot([Ep1(101:200,101);Ep1(1,101:200)';flip(diag(Ep1(101:200,101:200)));diag(Ep1(101:-1:2,101:200));EEEc1(1:101,1);Ep1(101,200:-1:101)'],'k','LineWidth',0.8)
% hold on
% plot([Ep2(101:200,101);Ep2(1,101:200)';flip(diag(Ep2(101:200,101:200)));diag(Ep2(101:-1:2,101:200));Ep2(1:101,1);Ep2(101,200:-1:101)'],'k','LineWidth',0.8)
% plot([Ep3(101:200,101);Ep3(1,101:200)';flip(diag(Ep3(101:200,101:200)));diag(Ep3(101:-1:2,101:200));Ep3(1:101,1);Ep3(101,200:-1:101)'],'k','LineWidth',0.8)
% plot([Ep4(101:200,101);Ep4(1,101:200)';flip(diag(Ep4(101:200,101:200)));diag(Ep4(101:-1:2,101:200));Ep4(1:101,1);Ep4(101,200:-1:101)'],'k','LineWidth',0.8)
% xticks([0 N/2 N 3*N/2 2*N 5*N/2 3*N])
% xticklabels({'\Gamma', 'X', 'M','\Gamma','M', 'Y', '\Gamma'})
% xlim([0 3*N])
% ax1.XGrid='on'
% ax1.GridAlpha=1;
% ax1.FontName='Times New Roman';
% ax1.FontSize=12;

% figure
% ax1=gca;
% plot([EEEh1(101:200,101);EEEh1(1,101:200)';flip(diag(EEEh1(101:200,101:200)));diag(EEEh1(101:-1:2,101:200));EEEh1(1:101,1);EEEh1(101,200:-1:101)'],'k','LineWidth',0.8)
% hold on
% plot([EEEh2(101:200,101);EEEh2(1,101:200)';flip(diag(EEEh2(101:200,101:200)));diag(EEEh2(101:-1:2,101:200));EEEh2(1:101,1);EEEh2(101,200:-1:101)'],'k','LineWidth',0.8)
% plot([EEEh3(101:200,101);EEEh3(1,101:200)';flip(diag(EEEh3(101:200,101:200)));diag(EEEh3(101:-1:2,101:200));EEEh3(1:101,1);EEEh3(101,200:-1:101)'],'k','LineWidth',0.8)
% plot([EEEh4(101:200,101);EEEh4(1,101:200)';flip(diag(EEEh4(101:200,101:200)));diag(EEEh4(101:-1:2,101:200));EEEh4(1:101,1);EEEh4(101,200:-1:101)'],'k','LineWidth',0.8)
% xticks([0 N/2 N 3*N/2 2*N 5*N/2 3*N])
% xticklabels({'\Gamma', 'X', 'M','\Gamma','M', 'Y', '\Gamma'})
% xlim([0 3*N])
% ax1.XGrid='on'
% ax1.GridAlpha=1;
% ax1.FontName='Times New Roman';
% ax1.FontSize=12;
% 
% figure
% ax1=gca;
% plot([EE1(101:200,101);EE1(1,101:200)';flip(diag(EE1(101:200,101:200)));diag(EE1(101:-1:2,101:200));EE1(1:101,1);EE1(101,200:-1:101)'],'k','LineWidth',0.8)
% hold on
% plot([EE2(101:200,101);EE2(1,101:200)';flip(diag(EE2(101:200,101:200)));diag(EE2(101:-1:2,101:200));EE2(1:101,1);EE2(101,200:-1:101)'],'k','LineWidth',0.8)
% plot([EE3(101:200,101);EE3(1,101:200)';flip(diag(EE3(101:200,101:200)));diag(EE3(101:-1:2,101:200));EE3(1:101,1);EE3(101,200:-1:101)'],'k','LineWidth',0.8)
% plot([EE4(101:200,101);EE4(1,101:200)';flip(diag(EE4(101:200,101:200)));diag(EE4(101:-1:2,101:200));EE4(1:101,1);EE4(101,200:-1:101)'],'k','LineWidth',0.8)
% xticks([0 N/2 N 3*N/2 2*N 5*N/2 3*N])
% xticklabels({'\Gamma', 'X', 'M','\Gamma','M', 'Y', '\Gamma'})
% xlim([0 3*N])
% ax1.XGrid='on'
% ax1.GridAlpha=1;
% ax1.FontName='Times New Roman';
% ax1.FontSize=12;
%% 
highsymm=1:2*N+4;
figure
ax1=gca;
EE1spg=[EE1(N/2:N,N/2) ; EE1(N,N/2:N)' ; EE1(N:-1:N/2,N) ; EE1(N/2,N:-1:N/2)'];
EE2spg=[EE2(N/2:N,N/2) ; EE2(N,N/2:N)' ; EE2(N:-1:N/2,N) ; EE2(N/2,N:-1:N/2)'];
EE3spg=[EE3(N/2:N,N/2) ; EE3(N,N/2:N)' ; EE3(N:-1:N/2,N) ; EE3(N/2,N:-1:N/2)'];
EE4spg=[EE4(N/2:N,N/2) ; EE4(N,N/2:N)' ; EE4(N:-1:N/2,N) ; EE4(N/2,N:-1:N/2)'];
plot(EE1spg,'r','LineWidth',0.8)
hold on
plot(EE2spg,'r','LineWidth',0.8)
plot(EE3spg,'b','LineWidth',0.8)
plot(EE4spg,'b','LineWidth',0.8)
% overlap_indicesEE1 = (EE1spg == EE3spg);
% plot(highsymm(overlap_indicesEE1),EE1spg(overlap_indicesEE1), 'k', 'LineWidth', 0.8);
% overlap_indicesEE2 = (EE2spg == EE4spg);
% plot(highsymm(overlap_indicesEE2),EE2spg(overlap_indicesEE2), 'k', 'LineWidth', 0.8);

xticks([0 N/2 N 3*N/2 2*N])
xticklabels({'\Gamma', 'X', 'M', 'Y', '\Gamma'})
xlim([0 2*N])
ax1.XGrid='on'
ax1.GridAlpha=1;
ax1.FontName='Times New Roman';
ax1.FontSize=12;
%% 

figure
ax2=gca;
Ep1spg=[Ep1(N/2:N,N/2) ; Ep1(N,N/2:N)' ; Ep1(N:-1:N/2,N) ; Ep1(N/2,N:-1:N/2)'];
Ep2spg=[Ep2(N/2:N,N/2) ; Ep2(N,N/2:N)' ; Ep2(N:-1:N/2,N) ; Ep2(N/2,N:-1:N/2)'];
Ep3spg=[Ep3(N/2:N,N/2) ; Ep3(N,N/2:N)' ; Ep3(N:-1:N/2,N) ; Ep3(N/2,N:-1:N/2)'];
Ep4spg=[Ep4(N/2:N,N/2) ; Ep4(N,N/2:N)' ; Ep4(N:-1:N/2,N) ; Ep4(N/2,N:-1:N/2)'];
plot(Ep1spg,'r','LineWidth',0.8)
hold on
plot(Ep2spg,'r','LineWidth',0.8)
plot(Ep3spg,'b','LineWidth',0.8)
plot(Ep4spg,'b','LineWidth',0.8)
% overlap_indicesEp1 = (Ep1spg == Ep3spg);
% plot(highsymm(overlap_indicesEp1),Ep1spg(overlap_indicesEp1), 'k', 'LineWidth', 0.8);
% overlap_indicesEp2 = (Ep2spg == Ep4spg);
% plot(highsymm(overlap_indicesEp2),Ep2spg(overlap_indicesEp2), 'k', 'LineWidth', 0.8);
xticks([0 N/2 N 3*N/2 2*N])
xticklabels({'\Gamma', 'X', 'M', 'Y', '\Gamma'})
xlim([0 2*N])
ax2.XGrid='on'
ax2.GridAlpha=1;
ax2.FontName='Times New Roman';
ax2.FontSize=12;

% title(['FFLO' '\rm{-phase gap structure}'])
% title(['\it{s}' '\rm{-wave gap structure}'])
% title(['\it{p_x}' '\rm{-wave gap structure}'])
% title(['\it{p_y}' '\rm{-wave gap structure}'])
% title('\it{p}\rm{^{+}_{\uparrow\uparrow}-wave gap structure}')
%% 

highsymm=1:2*N+4;
figure
ax1=gca;
EEE1spg=[EEE1(N/2:N,N/2) ; EEE1(N,N/2:N)' ; EEE1(N:-1:N/2,N) ; EEE1(N/2,N:-1:N/2)'];
EEE2spg=[EEE2(N/2:N,N/2) ; EEE2(N,N/2:N)' ; EEE2(N:-1:N/2,N) ; EEE2(N/2,N:-1:N/2)'];
EEE3spg=[EEE3(N/2:N,N/2) ; EEE3(N,N/2:N)' ; EEE3(N:-1:N/2,N) ; EEE3(N/2,N:-1:N/2)'];
EEE4spg=[EEE4(N/2:N,N/2) ; EEE4(N,N/2:N)' ; EEE4(N:-1:N/2,N) ; EEE4(N/2,N:-1:N/2)'];
plot(EEE1spg,'r','LineWidth',0.8)
hold on
plot(EEE2spg,'b','LineWidth',0.8)
plot(EEE3spg,'r','LineWidth',0.8)
plot(EEE4spg,'b','LineWidth',0.8)
% overlap_indicesEE1 = (EE1spg == EE3spg);
% plot(highsymm(overlap_indicesEE1),EE1spg(overlap_indicesEE1), 'k', 'LineWidth', 0.8);
% overlap_indicesEE2 = (EE2spg == EE4spg);
% plot(highsymm(overlap_indicesEE2),EE2spg(overlap_indicesEE2), 'k', 'LineWidth', 0.8);

xticks([0 N/2 N 3*N/2 2*N])
xticklabels({'\Gamma', 'X', 'M', 'Y', '\Gamma'})
xlim([0 2*N])
ax1.XGrid='on'
ax1.GridAlpha=1;
ax1.FontName='Times New Roman';
ax1.FontSize=12;
