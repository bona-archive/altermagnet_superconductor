clear all
load('CustomColormap.mat')
%%%% t,J,mu,del are parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=1;
J=0.4;
mu=-1.8;
m=1;
beta=5000;

%conjugate of del are on third, fourth row
NOGx=200; %%% Number Of Grid x
NOGy=200; %%% Number Of Grid y
NF=1/(NOGx*NOGy); %%% Normalization Factor

kkx=linspace(-pi,pi,NOGx); %%%% k-space %%%%
kky=linspace(-pi,pi,NOGy);

%%%% zero %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Etotpp11=0;
%%%% trial self consistent p-wave interaction params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vp=-9.2;Vp=0;
Dpp11=12.2382; Dpp11=0;% solved by aRuO2AMSC (self consistent input Dpp11)
Dpm11=12.2382; Dpm11=0;% solved by aRuO2AMSC (self consistent input Dpp11)
%%%% phase difference of 11 cooper pair and 22 cooper pair %%%%
phi_width_num=100;
pphi=linspace(-2*pi,2*pi,phi_width_num);
%%%% occupied energy for self-consistent solution array %%%%
EEtotpp11f=zeros(phi_width_num,1);
EEtotpm11f=zeros(phi_width_num,1);
%%%% cell arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rgb_values=[];
rgb_C4T_values=[];
rgb_valuesp=[];
rgb_C4T_valuesp=[];
% for s=1:length(pphi) %%%% p+ TRS broken(pauli(1)) strained system
%     phi=pphi(s);
%     s
% b=0.009
b=0;
for i=1:length(kkx)
    for j=1:length(kky)
        kx=kkx(i);
        ky=kky(j);
        kxm=-kkx(i);
        kym=-kky(j);

        Zeeman=m*pauli(3);
        Zeemanm=m*pauli(3);
        Alter=J*kx*ky*pauli(3);
        Alterm=J*kxm*kym*pauli(3);
        SymmBreaker3xxxy=kx^3*pauli(2)+ky^3*pauli(1);
        SymmBreakerm3xxxy=kxm^3*pauli(2)+kym^3*pauli(1);
        he=(-t*(kx^2+ky^2)-mu)*pauli(0)+Alter+b*SymmBreaker3xxxy;
        hh=(-t*(kxm^2+kym^2)-mu)*pauli(0)+Alterm+b*SymmBreakerm3xxxy;
        Hbdg=blkdiag(he,-transpose(hh));

        [u1,d1]=eig(Hbdg);
        psi3up(i,j)=abs(u1(1,3));
        psi3down(i,j)=abs(u1(2,3));
        psi4up(i,j)=abs(u1(1,4));
        psi4down(i,j)=abs(u1(2,4));
        sz3(i,j)=u1(1:2,3)'*pauli(3)*u1(1:2,3);sz3sphere(NOGx*(i-1)+j)=sz3(i,j);
        sz4(i,j)=u1(1:2,4)'*pauli(3)*u1(1:2,4);sz4sphere(NOGx*(i-1)+j)=sz4(i,j);
        sx3(i,j)=u1(1:2,3)'*pauli(1)*u1(1:2,3);sx3sphere(NOGx*(i-1)+j)=sx3(i,j);
        sx4(i,j)=u1(1:2,4)'*pauli(1)*u1(1:2,4);sx4sphere(NOGx*(i-1)+j)=sx4(i,j);
        sy3(i,j)=u1(1:2,3)'*pauli(2)*u1(1:2,3);sy3sphere(NOGx*(i-1)+j)=sy3(i,j);
        sy4(i,j)=u1(1:2,4)'*pauli(2)*u1(1:2,4);sy4sphere(NOGx*(i-1)+j)=sy4(i,j);
        rgb_values=[rgb_values;(kkx(i)/pi+1)/2 0 (kky(j)/pi+1)/2];
        rgb_C4T_values=[rgb_C4T_values;(kkx(i)/pi+1)/2 0 (kky(j)/pi+1)/2];
        en1(i,j)=d1(1,1);en2(i,j)=d1(2,2);en3(i,j)=d1(3,3);en4(i,j)=d1(4,4);
        H0(i,j,:,:)=(t*(kkx(i)^2+kky(j)^2)-mu)*pauli(0)+J*kkx(i)*kky(j)*pauli(3);
        
        gpp=kx+1i*ky;
        gpm=kx-1i*ky;
        D11=Dpp11*gpp;
        D22=Dpm11*gpm;
        Hdp=Vp*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];
        Hbdg=Hbdg+Hdp;
        [u1,d1]=eig(Hbdg);
        psi3upp(i,j)=abs(u1(1,3));
        psi3downp(i,j)=abs(u1(2,3));
        psi4upp(i,j)=abs(u1(1,4));
        psi4downp(i,j)=abs(u1(2,4));
        sz3p(i,j)=u1(1:2,3)'*pauli(3)*u1(1:2,3);sz3psphere(NOGx*(i-1)+j)=sz3(i,j);
        sz4p(i,j)=u1(1:2,4)'*pauli(3)*u1(1:2,4);sz4psphere(NOGx*(i-1)+j)=sz4(i,j);
        sx3p(i,j)=u1(1:2,3)'*pauli(1)*u1(1:2,3);sx3psphere(NOGx*(i-1)+j)=sx3(i,j);
        sx4p(i,j)=u1(1:2,4)'*pauli(1)*u1(1:2,4);sx4psphere(NOGx*(i-1)+j)=sx4(i,j);
        sy3p(i,j)=u1(1:2,3)'*pauli(2)*u1(1:2,3);sy3psphere(NOGx*(i-1)+j)=sy3(i,j);
        sy4p(i,j)=u1(1:2,4)'*pauli(2)*u1(1:2,4);sy4psphere(NOGx*(i-1)+j)=sy4(i,j);
        rgb_valuesp=[rgb_valuesp;(kkx(i)/pi+1)/2 0 (kky(j)/pi+1)/2];
        rgb_C4T_valuesp=[rgb_C4T_valuesp;(kkx(i)/pi+1)/2 0 (kky(j)/pi+1)/2];
        en1p(i,j)=d1(1,1);en2p(i,j)=d1(2,2);en3p(i,j)=d1(3,3);en4p(i,j)=d1(4,4);
        
    end
end

%%%% before operator
figure;[a b c]=sphere(NOGx);h=surf(a,b,c);set(h,'EdgeColor','none','FaceAlpha','0.7');axis equal;hold on %%블로흐면 나타내기
scatter3(sx3sphere,sy3sphere,sz3sphere,50,rgb_values,'filled')
scatter3(sx4sphere,sy4sphere,sz4sphere,50,rgb_values,'filled')

figure
[xk, yk]=meshgrid(kkx,kky);
scatter3(xk(:),yk(:),zeros(1,NOGx^2),50,rgb_values,'filled')
xlabel('k_x');ylabel('k_y');view(2);
xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});

%%%% C4T operation
figure;[a b c]=sphere(NOGx);h=surf(a,b,c);set(h,'EdgeColor','none','FaceAlpha','0.7');axis equal;hold on %%블로흐면 나타내기
scatter3(sy3sphere,-sx3sphere,-sz3sphere,50,rgb_values,'filled')
scatter3(sy4sphere,-sx4sphere,-sz4sphere,50,rgb_values,'filled');title('Acting C4T operator')

figure
[xk, yk]=meshgrid(kkx,kky);
scatter3(yk(:),-xk(:),zeros(1,NOGx^2),50,rgb_values,'filled');title('Acting C4T operator')
xlabel('k_x');ylabel('k_y');view(2);
xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});

%%%% Fermi Surface, z Spin expectation value distribution in BZ
figure;surf(kkx,kky,en3,Linestyle='none');hold on
xlim([-pi pi]);ylim([-pi pi]);zlim([0.2 1]);ccc=colorbar;clim([0 2]);ccc.Label.String='Energy';
xlabel('k_x');ylabel('k_y');zlabel('Energy');
xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
zticks([0.2 1]);zticklabels({'0', '1'});
surf(kkx,kky,en4,Linestyle='none');
xlim([-pi pi]);ylim([-pi pi]);zlim([0.2 1]);
ccccc=colorbar;clim([0 2]);ccccc.Label.String='Energy';
xlabel('k_x');ylabel('k_y');zlabel('Energy');
hold on;quiver3(kkx,kky,en3,sx3,sy3,sz3,'r');hold on;quiver3(kkx,kky,en4,sx4,sy4,sz4,'r');

figure;surf(kkx,kky,en3,Linestyle='none');hold on
xlim([-pi pi]);ylim([-pi pi]);zlim([0.2 1]);ccc=colorbar;clim([0 2]);ccc.Label.String='Energy';
xlabel('k_x');ylabel('k_y');zlabel('Energy');
xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
zticks([0.2 1]);zticklabels({'0', '1'});
surf(kkx,kky,en4,Linestyle='none');
xlim([-pi pi]);ylim([-pi pi]);zlim([0.2 1]);
ccccc=colorbar;clim([0 2]);ccccc.Label.String='Energy';
xlabel('k_x');ylabel('k_y');zlabel('Energy');
hold on;quiver3(kkx,kky,en3,sx3,sy3,sz3,'r');hold on;quiver3(kkx,kky,en4,sx4,sy4,sz4,'r');
view(90,0)


C3=sz3;C4=sz4;
figure;surf(kkx,kky,en3,C3,Linestyle='none');hold on;
xlim([-pi pi]);ylim([-pi pi]);zlim([0.2 1]);colormap(CustomColormap);cc=colorbar;clim([-1 1]);cc.Label.String='\langle \sigma_z \rangle';
xlabel('k_x');ylabel('k_y');zlabel('Energy');
xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
zticks([0.2 1]);zticklabels({'0', '1'});
surf(kkx,kky,en4,C4,Linestyle='none');
xlim([-pi pi]);ylim([-pi pi]);zlim([0.2 1]);
colormap(CustomColormap);cccc=colorbar;clim([-1 1]);cccc.Label.String='\langle \sigma_z \rangle';
xlabel('k_x');ylabel('k_y');zlabel('Energy');
hold on;quiver3(kkx,kky,en3,sx3,sy3,zeros(NOGx),'r');hold on;quiver3(kkx,kky,en4,sx4,sy4,zeros(NOGx),'r');
%%%% x Spin expectation value distribution in BZ
C3=sx3;C4=sx4;
figure;surf(kkx,kky,en3,C3,Linestyle='none');hold on;
xlim([-pi pi]);ylim([-pi pi]);zlim([0.2 1]);colormap(CustomColormap);cc=colorbar;
clim([-1 1]);cc.Label.String='\langle \sigma_x \rangle';
xlabel('k_x');ylabel('k_y');zlabel('Energy');
xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
zticks([0.2 1]);zticklabels({'0', '1'});
surf(kkx,kky,en4,C4,Linestyle='none');
xlim([-pi pi]);ylim([-pi pi]);zlim([0.2 1]);
colormap(CustomColormap);cccc=colorbar;clim([-1 1]);cccc.Label.String='\langle \sigma_x \rangle';
xlabel('k_x');ylabel('k_y');zlabel('Energy');

%%%% y Spin expectation value distribution in BZ
C3=sy3;C4=sy4;
figure;surf(kkx,kky,en3,C3,Linestyle='none');hold on;
xlim([-pi pi]);ylim([-pi pi]);zlim([0.2 1]);
colormap(CustomColormap);cc=colorbar;clim([-1 1]);cc.Label.String='\langle \sigma_y \rangle';
xlabel('k_x');ylabel('k_y');zlabel('Energy');
xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
zticks([0.2 1]);zticklabels({'0', '1'});
surf(kkx,kky,en4,C4,Linestyle='none');
xlim([-pi pi]);ylim([-pi pi]);zlim([0.2 1]);
colormap(CustomColormap);cccc=colorbar;clim([-1 1]);cccc.Label.String='\langle \sigma_y \rangle';
xlabel('k_x');ylabel('k_y');zlabel('Energy');

%%%% wavefunction
figure
surf(real(psi3up))

%%%% Fermi surface of altermagnetism
figure
surf(kkx,kky,H0(:,:,1,1),'EdgeColor','none','FaceAlpha',1,'FaceColor','r')
hold on
surf(kkx,kky,H0(:,:,2,2),'EdgeColor','none','FaceAlpha',1,'FaceColor','b')
xlim([-pi pi]);ylim([-pi pi]);zlim([0 0.1]);
xlabel('k_x');ylabel('k_y');zlabel('Energy');
xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
zticks([0.2 1]);zticklabels({'0', '1'});

% figure;quiver3(zeros(1,NOGx),zeros(1,NOGy),zeros(NOGx,NOGy),sx3,sy3,sz3,'Color',color3);hold on;quiver3(zeros(1,NOGx),zeros(1,NOGy),zeros(NOGx,NOGy),sx4,sy4,sz4,'Color',color4);
% xlim([-pi pi]);ylim([-pi pi]);xlabel('k_x');ylabel('k_y');
% xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});

%%% p-wave SC
%%%% before operator
figure;[a b c]=sphere(NOGx);h=surf(a,b,c);set(h,'EdgeColor','none','FaceAlpha','0.7');axis equal;hold on %%블로흐면 나타내기
scatter3(sx3psphere,sy3psphere,sz3psphere,50,rgb_valuesp,'filled')
scatter3(sx4psphere,sy4psphere,sz4psphere,50,rgb_valuesp,'filled')

%%%% C4T operation
figure;[a b c]=sphere(NOGx);h=surf(a,b,c);set(h,'EdgeColor','none','FaceAlpha','0.7');axis equal;hold on %%블로흐면 나타내기
scatter3(sy3psphere,-sx3psphere,-sz3psphere,50,rgb_valuesp,'filled')
scatter3(sy4psphere,-sx4psphere,-sz4psphere,50,rgb_valuesp,'filled');title('Acting C4T operator')

figure;surf(kkx,kky,en3p,Linestyle='none');hold on
xlim([-pi pi]);ylim([-pi pi]);%zlim([0.2 1]);
% ccc=colorbar;clim([0 2]);ccc.Label.String='Energy';
xlabel('k_x');ylabel('k_y');zlabel('Energy');
xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
%zticks([0.2 1]);zticklabels({'0', '1'});
surf(kkx,kky,en4p,Linestyle='none');
xlim([-pi pi]);ylim([-pi pi]);%zlim([0.2 1]);
% ccccc=colorbar;clim([0 2]);ccccc.Label.String='Energy';
xlabel('k_x');ylabel('k_y');zlabel('Energy');
hold on;quiver3(kkx,kky,en3p,sx3p,sy3p,sz3p,'r');hold on;quiver3(kkx,kky,en4p,sx4p,sy4p,sz4p,'r');

figure;surf(kkx,kky,en3p,Linestyle='none');hold on
xlim([-pi pi]);ylim([-pi pi]);%zlim([0.2 1]);
ccc=colorbar;clim([0 2]);ccc.Label.String='Energy';
xlabel('k_x');ylabel('k_y');zlabel('Energy');
xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
%zticks([0.2 1]);zticklabels({'0', '1'});
surf(kkx,kky,en4p,Linestyle='none');
xlim([-pi pi]);ylim([-pi pi]);%zlim([0.2 1]);
ccccc=colorbar;clim([0 2]);ccccc.Label.String='Energy';
xlabel('k_x');ylabel('k_y');zlabel('Energy');
hold on;quiver3(kkx,kky,en3p,sx3p,sy3p,sz3p,'r');hold on;quiver3(kkx,kky,en4p,sx4p,sy4p,sz4p,'r');
view(90,0)


C3=sz3p;C4=sz4p;
figure;surf(kkx,kky,en3p,C3,Linestyle='none');hold on;
xlim([-pi pi]);ylim([-pi pi]);%zlim([0.2 1]);
colormap(CustomColormap);
cc=colorbar;clim([-1 1]);cc.Label.String='\langle \sigma_z \rangle';
xlabel('k_x');ylabel('k_y');zlabel('Energy');
xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
%zticks([0.2 1]);zticklabels({'0', '1'});
surf(kkx,kky,en4p,C4,Linestyle='none');
xlim([-pi pi]);ylim([-pi pi]);%zlim([0.2 1]);
colormap(CustomColormap);cccc=colorbar;clim([-1 1]);cccc.Label.String='\langle \sigma_z \rangle';
xlabel('k_x');ylabel('k_y');zlabel('Energy');
hold on;quiver3(kkx,kky,en3p,sx3p,sy3p,zeros(NOGx),'r');hold on;quiver3(kkx,kky,en4p,sx4p,sy4p,zeros(NOGx),'r');
%%%% x Spin expectation value p-wave distribution in BZ
C3=sx3p;C4=sx4p;
figure;surf(kkx,kky,en3p,C3,Linestyle='none');hold on;
xlim([-pi pi]);ylim([-pi pi]);%zlim([0.2 1]);
colormap(CustomColormap);
cc=colorbar;clim([-1 1]);cc.Label.String='\langle \sigma_x \rangle';
xlabel('k_x');ylabel('k_y');zlabel('Energy');
xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
%zticks([0.2 1]);zticklabels({'0', '1'});
surf(kkx,kky,en4p,C4,Linestyle='none');
xlim([-pi pi]);ylim([-pi pi]);%zlim([0.2 1]);
colormap(CustomColormap);cccc=colorbar;clim([-1 1]);cccc.Label.String='\langle \sigma_x \rangle';
xlabel('k_x');ylabel('k_y');zlabel('Energy');

%%%% y Spin expectation p-wave value distribution in BZ
C3=sy3p;C4=sy4p;
figure;surf(kkx,kky,en3p,C3,Linestyle='none');hold on;
xlim([-pi pi]);ylim([-pi pi]);%zlim([0.2 1]);
colormap(CustomColormap);cc=colorbar;clim([-1 1]);cc.Label.String='\langle \sigma_y \rangle';
xlabel('k_x');ylabel('k_y');zlabel('Energy');
xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
%zticks([0.2 1]);zticklabels({'0', '1'});
surf(kkx,kky,en4p,C4,Linestyle='none');
xlim([-pi pi]);ylim([-pi pi]);%zlim([0.2 1]);
colormap(CustomColormap);cccc=colorbar;clim([-1 1]);cccc.Label.String='\langle \sigma_y \rangle';
xlabel('k_x');ylabel('k_y');zlabel('Energy');

%%%% z Spin expectation p-wave value distribution in BZ
C3=sz3p;C4=sz4p;
figure;surf(kkx,kky,en3p,C3,Linestyle='none');hold on;
xlim([-pi pi]);ylim([-pi pi]);%zlim([0.2 1]);
colormap(CustomColormap);cc=colorbar;clim([-1 1]);cc.Label.String='\langle \sigma_z \rangle';
xlabel('k_x');ylabel('k_y');zlabel('Energy');
xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
yticks([-pi -pi/2 0 pi/2 pi]);yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
%zticks([0.2 1]);zticklabels({'0', '1'});
surf(kkx,kky,en4p,C4,Linestyle='none');
xlim([-pi pi]);ylim([-pi pi]);%zlim([0.2 1]);
colormap(CustomColormap);cccc=colorbar;clim([-1 1]);cccc.Label.String='\langle \sigma_z \rangle';
xlabel('k_x');ylabel('k_y');zlabel('Energy');


%%%% wavefunction p-wave
figure
surf(real(psi3upp))

