clear all

%%%% t,J,mu,del are parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=1;
J=1;
mu=0.7;
m=0.3;
beta=5000;

%conjugate of del are on third, fourth row
NOGx=100; %%% Number Of Grid x
NOGy=100; %%% Number Of Grid y
NF=1/(NOGx*NOGy); %%% Normalization Factor

kkx=linspace(-pi,pi,NOGx); %%%% k-space %%%%
kky=linspace(-pi,pi,NOGy);

%%%% number of trial interaction params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_width_num=41;
VVS=linspace(-9.2,-2.2,param_width_num);
VVP=-11.4-VVS;

a=1
Vs=VVS(a);
Vs
Vp=-11.4-Vs;
Vp
legendLabels{a}=sprintf('V_{s}=%.2f', Vs);
legendLabelp{a}=sprintf('V_{p}=%.2f', Vp);

for i=1:length(kkx)
    for j=1:length(kky)
        kx=kkx(i);
        ky=kky(j);
        kxm=-kkx(i);
        kym=-kky(j);

        % SymmBreaker0=m*pauli(1);
        % SymmBreakerm0=m*pauli(1);
        SymmBreaker1=kx*pauli(1)+ky*pauli(2);
        SymmBreakerm1=kxm*pauli(1)+kym*pauli(2);
        % SymmBreaker2=kx^2*pauli(1)+ky^2*pauli(2);
        % SymmBreakerm2=kxm^2*pauli(1)+kym^2*pauli(2);
        % SymmBreaker2cross=kx^2*pauli(2)+ky^2*pauli(1);
        % SymmBreakerm2cross=kxm^2*pauli(2)+kym^2*pauli(1);
        SymmBreaker3xxxx=kx^3*pauli(1);
        SymmBreakerm3xxxx=kxm^3*pauli(1);
        SymmBreaker3yyyy=ky^3*pauli(2);
        SymmBreakerm3yyyy=kym^3*pauli(2);
        SymmBreaker3xxyx=kx^2*ky*pauli(1);
        SymmBreakerm3xxyx=kxm^2*kym*pauli(1);
        SymmBreaker3xxyy=kx^2*ky*pauli(2);
        SymmBreakerm3xxyy=kxm^2*kym*pauli(2);
        SymmBreaker3xyyx=kx*ky^2*pauli(1);
        SymmBreakerm3xyyx=kxm*kym^2*pauli(1);
        SymmBreaker3xyyy=kx*ky^2*pauli(2);
        SymmBreakerm3xyyy=kxm*kym^2*pauli(2);
        SymmBreaker3yyyy=ky^3*pauli(2);
        SymmBreakerm3yyyy=kym^3*pauli(2);

        he=(t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3)+SymmBreaker3xyyy;
        hh=(t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3)+SymmBreakerm3xyyy;
        Hbdg=blkdiag(he,-transpose(hh));

        [u,d]=eig(Hbdg);
        F=diag(1./(exp(diag(d)*beta)+1));
        correls(:,:,i,j)=(u*F*u');

        en1(i,j)=d(1,1);en2(i,j)=d(2,2);en3(i,j)=d(3,3);en4(i,j)=d(4,4);
    end
end


%Ds=1;
%for i=1:length(kkx)
%    for j=1:length(kky)
%        kx=kkx(i);
%        ky=kky(j);
%        kxm=-kkx(i);
%        kym=-kky(j);
%
%        Hds=Vs*[0,0,0,Ds;0,0,-Ds,0;0,-conj(Ds),0,0;conj(Ds),0,0,0];

%        he=(t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3)+m*pauli(2);
%        hh=(t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3)+m*pauli(2);
%        Hbdg=blkdiag(he,-transpose(hh));

%        gs=1;

%        Hbdg=Hbdg+(Hds);

%        [u,d]=eig(Hbdg);
%        F=diag(1./(exp(diag(d)*beta)+1));
%        correls(:,:,i,j)=(u*F*u');

%        en1s(i,j)=d(1,1);en2s(i,j)=d(2,2);en3s(i,j)=d(3,3);en4s(i,j)=d(4,4);
%    end
%end

%Dpp11=1;
%for i=1:length(kkx)
%    for j=1:length(kky)
%        kx=kkx(i);
%        ky=kky(j);
%        kxm=-kkx(i);
%        kym=-kky(j);

%        he=(t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3)+m*pauli(2);
%        hh=(t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3)+m*pauli(2);
%        Hbdg=blkdiag(he,-transpose(hh));

%        gpp=kx+1i*ky;
%        gpm=kx-1i*ky;

%        Dpm11=0;    %D11=Dpp11*gpp solution
%        Dpm22=0;    %D22=Dpp22*gpp solution
%        Dpp22=Dpp11;    % p+ solution
%        D11=Dpp11*gpp+Dpm11*gpm;
%        D22=Dpp22*gpp+Dpm22*gpm;
%        Hdp=Vp*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];

%        Hbdg=Hbdg+(Hdp);

%        [u,d]=eig(Hbdg);
%        F=diag(1./(exp(diag(d)*beta)+1));
%        correlp(:,:,i,j)=(u*F*u');

%        en1p(i,j)=d(1,1);en2p(i,j)=d(2,2);en3p(i,j)=d(3,3);en4p(i,j)=d(4,4);
%    end
%end


%subplot(1,2,1)
%surf(squeeze(kkx),squeeze(kky),en1s(:,:))
%hold on
%surf(squeeze(kkx),squeeze(kky),en2s(:,:))
%surf(squeeze(kkx),squeeze(kky),en3s(:,:))
%surf(squeeze(kkx),squeeze(kky),en4s(:,:))
%subplot(1,2,2)
%surf(squeeze(kkx),squeeze(kky),en1p(:,:))
%hold on
%surf(squeeze(kkx),squeeze(kky),en2p(:,:))
%surf(squeeze(kkx),squeeze(kky),en3p(:,:))
%surf(squeeze(kkx),squeeze(kky),en4p(:,:))

figure
%surf(squeeze(kkx),squeeze(kky),en1(:,:))
hold on
surf(squeeze(kkx),squeeze(kky),en2(:,:))
surf(squeeze(kkx),squeeze(kky),en3(:,:))
%surf(squeeze(kkx),squeeze(kky),en4(:,:))
camlight