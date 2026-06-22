clear all

%%%% t,J,mu,del are parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=1;
J=1;
mu=0.7;
m=0.3;
beta=5000;

%conjugate of del are on third, fourth row
NOGx=15; %%% Number Of Grid x
NOGy=15; %%% Number Of Grid y
NF=1/(NOGx*NOGy); %%% Normalization Factor

kkx=linspace(-pi,pi,NOGx); %%%% k-space %%%%
kky=linspace(-pi,pi,NOGy);

%%%% zero %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Etotpp11=0;
%%%% trial self consistent p-wave interaction params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vp=-2.2;
% Dpp11=2.4116; % solved by aRuO2AMSC (self consistent input Dpp11)
% Dpm11=2.4116; % solved by aRuO2AMSC (self consistent input Dpp11)
Vp=-9.2;
Dpp11=11.8179;
Dpm11=11.8179;
%%%% phase difference of 11 cooper pair and 22 cooper pair %%%%
phi_width_num=100;
pphi=linspace(-2*pi,2*pi,phi_width_num);
%%%% occupied energy for self-consistent solution array %%%%
EEtotpp11f=zeros(phi_width_num,1);
EEtotpm11f=zeros(phi_width_num,1);
%%%% cell arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s=1:length(pphi) %%%% p+ TRS broken(pauli(1)) strained system
    phi=pphi(s);
    s
    for i=1:length(kkx)
        for j=1:length(kky)
            kx=kkx(i);
            ky=kky(j);
            kxm=-kkx(i);
            kym=-kky(j);
            

            SymmBreaker1=m*pauli(1);
            SymmBreakerm1=m*pauli(1);
            SymmBreaker2=kx*pauli(1)+ky*pauli(2);
            SymmBreakerm2=kxm*pauli(1)+kym*pauli(2);
            he=(t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3)+SymmBreaker2;
            hh=(t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3)+SymmBreakerm2;
            Hbdg=blkdiag(he,-transpose(hh));

            gpp=kx+1i*ky;
            gpm=kx-1i*ky;

            Dpm11=0;    %D11=Dpp11*gpp solution
            Dpm22=0;    %D22=Dpp22*gpp solution
            Dpp22=Dpp11;    % p+ solution
            D11=Dpp11*gpp+Dpm11*gpm;
            D22=(Dpp22*gpp+Dpm22*gpm)*exp(1i*phi); %but different phase between wavefunctions
            Hdp=Vp*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];

            Hbdg=Hbdg+(Hdp);

            [u,d]=eig(Hbdg);
            F=diag(1./(exp(diag(d)*beta)+1));
            correlpp(:,:,i,j)=(u*F*u');
            en1pp(i,j)=d(1,1);en2pp(i,j)=d(2,2);en3pp(i,j)=d(3,3);en4pp(i,j)=d(4,4);
        end
    end
    Etotpp11=sum(en1pp(en1pp<0))+sum(en2pp(en2pp<0))+sum(en3pp(en3pp<0))+sum(en4pp(en4pp<0));
    EEtotpp11f(s)=Etotpp11;
    Etotpp11=0;
end


for p=1:length(pphi) % p+ TRS broken(pauli(1)) strained system
    phi=pphi(p);
    p
    for i=1:length(kkx)
        for j=1:length(kky)
            kx=kkx(i);
            ky=kky(j);
            kxm=-kkx(i);
            kym=-kky(j);

            SymmBreaker1=m*pauli(1);
            SymmBreakerm1=m*pauli(1);
            SymmBreaker2=kx*pauli(1)+ky*pauli(2);
            SymmBreakerm2=kxm*pauli(1)+kym*pauli(2);
            he=(t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3)
            hh=(t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3)
            Hbdg=blkdiag(he,-transpose(hh));

            gpp=kx+1i*ky;
            gpm=kx-1i*ky;

            Dpp11=0;    %D11=Dpm11*gpm solution
            Dpp22=0;    %D22=Dpm22*gpm solution
            Dpm22=Dpm11;    % p+ solution
            D11=Dpp11*gpp+Dpm11*gpm;
            D22=(Dpp22*gpp+Dpm22*gpm)*exp(1i*phi);
            Hdp=Vp*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];

            Hbdg=Hbdg+(Hdp);

            [u,d]=eig(Hbdg);
            F=diag(1./(exp(diag(d)*beta)+1));
            correlpm(:,:,i,j)=(u*F*u');

            en1pm(i,j)=d(1,1);en2pm(i,j)=d(2,2);en3pm(i,j)=d(3,3);en4pm(i,j)=d(4,4);
        end
    end
    Etotpm11=sum(en1pm(en1pm<0))+sum(en2pm(en2pm<0))+sum(en3pm(en3pm<0))+sum(en4pm(en4pm<0));
    EEtotpm11f(p)=Etotpm11;
    Etotpm11=0;
end

%%%% Total occupied energy - phi plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(squeeze(pphi),squeeze(EEtotpp11f))
hold on
plot(squeeze(pphi),squeeze(EEtotpm11f), 'o')
xlabel('\phi')
ylabel('E_{tot}')
legend('p^{+} TRS broken occupied energy','p^{-} TRS broken occupied energy')