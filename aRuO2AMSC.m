clear all

%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=1;
J=1;
mu=-1;
m=0.3;
beta=5000;
NOGx=10; % Number Of Grid x
NOGy=10; % Number Of Grid y
NF=1/(NOGx*NOGy); % Normalization Factor
kkx=linspace(-pi,pi,NOGx);
kky=linspace(-pi,pi,NOGy);

%%%% input orderparams limit %%%%
DDslim=5;
DDpp11lim=15;
DDpp22lim=DDpp11lim;
DDpm11lim=15;
DDpm22lim=DDpm11lim;
DDpx11lim=15;
DDpx22lim=DDpx11lim;
DDpy11lim=15;
DDpy22lim=DDpy11lim;
DDmix11lim=15;
DDmix22lim=15;

%%%% number of widths of input orderparams %%%%
width_num=DDslim*100;
width_numpp11=DDpp11lim*100;
width_numpp22=DDpp22lim*100;
width_numpm11=DDpm11lim*100;
width_numpm22=DDpm22lim*100;
width_numpx11=DDpx11lim*100;
width_numpx22=width_numpx11;
width_numpy11=DDpy11lim*100;
width_numpy22=width_numpy11;
width_nummix11=DDmix11lim*100;
width_nummix22=DDmix22lim*100;

%%%% input orderparams %%%%
DDsinput=linspace(0,DDslim,width_num);
DDpp11input=linspace(0,DDpp11lim,width_numpp11);
DDpm11input=DDpp11input;
DDpp22input=linspace(0,DDpp22lim,width_numpp22);
DDpm22input=DDpp22input;
DDpx11input=linspace(0,DDpx11lim,width_numpx11);
DDpx22input=DDpx11input;
DDpy11input=linspace(0,DDpy11lim,width_numpy11);
DDpy22input=DDpy11input;
DDmix11input=linspace(0,DDmix11lim,width_nummix11);
DDmix22input=DDmix11input;

%%%% zero arrays %%%%
%%%% output orderparams receivers %%%%
DDsout=zeros(width_num,1);
DDpp11out=zeros(width_numpp11,1);
DDpp22out=zeros(width_numpp22,1);
DDpm11out=zeros(width_numpp11,1);
DDpm22out=zeros(width_numpp22,1);
DDpx11out=zeros(width_numpx11,1);
DDpx22out=zeros(width_numpx22,1);
DDpy11out=zeros(width_numpy11,1);
DDpy22out=zeros(width_numpy22,1);
DDmix11out=zeros(width_nummix11,1);
DDmix22out=zeros(width_nummix22,1);

%%%% occupied energy receivers for each output orderparams %%%%
EEtots=zeros(width_num,1);
EEtotpp11=zeros(width_numpp11,1);
EEtotpp22=zeros(width_numpp22,1);
EEtotpm11=zeros(width_numpp11,1);
EEtotpm22=zeros(width_numpp22,1);
EEtotpx11=zeros(width_numpx11,1);
EEtotpx22=zeros(width_numpx22,1);
EEtotpy11=zeros(width_numpy11,1);
EEtotpy22=zeros(width_numpy22,1);
EEtotmix11=zeros(width_nummix11,1);
EEtotmix22=zeros(width_nummix22,1);

%%%% output order params 
Dsout=0;
Dpp11out=0;
Dpm11out=0;
Dpp22out=0;
Dpm22out=0;
Dpx11out=0;
Dpy11out=0;
Dpx22out=0;
Dpy22out=0;
Dmix11out=0;
Dmix22out=0;

%%% number of trial interaction params %%%%
param_width_num=10;
L=-2;
R=0;
T=L+R;
VVS=linspace(L,R,param_width_num);
%%%% self-consistent occupied energy receivers %%%%
EEtotsf=zeros(param_width_num,1); 
EEtotpp11f=zeros(param_width_num,1);
EEtotpp22f=zeros(param_width_num,1);
EEtotpm11f=zeros(param_width_num,1);
EEtotpm22f=zeros(param_width_num,1);
EEtotpx11f=zeros(param_width_num,1);
EEtotpx22f=zeros(param_width_num,1);
EEtotpy11f=zeros(param_width_num,1);
EEtotpy22f=zeros(param_width_num,1);
EEtotmix11f=zeros(param_width_num,1);
EEtotmix22f=zeros(param_width_num,1);

%%%% self-consistent input orderparam indexes receivers %%%%
Indexs=zeros(param_width_num,1);
Indexpp11=zeros(param_width_num,1);
Indexpm11=zeros(param_width_num,1);
Indexpx11=zeros(param_width_num,1);
Indexpy11=zeros(param_width_num,1);
Indexmix11=zeros(param_width_num,1);
Indexmix22=zeros(param_width_num,1);

%%% cell arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DDDsout=cell(1,param_width_num);
EEEtots=cell(1,param_width_num);

DDDpp11out=cell(1,param_width_num);
DDDpp22out=cell(1,param_width_num);
EEEtotpp11=cell(1,param_width_num);

DDDpm11out=cell(1,param_width_num);
DDDpm22out=cell(1,param_width_num);
EEEtotpm11=cell(1,param_width_num);

DDDpx11out=cell(1,param_width_num);
DDDpx22out=cell(1,param_width_num);
EEEtotpx11=cell(1,param_width_num);

DDDpy11out=cell(1,param_width_num);
DDDpy22out=cell(1,param_width_num);
EEEtotpy11=cell(1,param_width_num);

DDDmix11out=cell(1,param_width_num);
DDDmix22out=cell(1,param_width_num);
EEEtotmix11=cell(1,param_width_num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon=1e-6; %%%% very small value epsilon for finding self consistent solution %%%%

%%% parameter boxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parboxpp=cell(1,param_width_num);
parboxpm=cell(1,param_width_num);

parboxpx=cell(1,param_width_num);
parboxpy=cell(1,param_width_num);

parboxmix=cell(1,param_width_num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for a=1:length(VVS)
    Vs=VVS(a);
    Vs
    Vp=-2-Vs;
    Vp
   legendLabels{a}=sprintf('V_{s}=%.2f', Vs);
   legendLabelp{a}=sprintf('V_{p}=%.2f', Vp);

   for s=1:length(DDsinput) % s-wave check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Ds=DDsinput(s);
       for i=1:length(kkx)
           for j=1:length(kky)
               kx=kkx(i);
               ky=kky(j);
               kxm=-kkx(i);
               kym=-kky(j);

               Hds=Vs*[0,0,0,Ds;0,0,-Ds,0;0,-conj(Ds),0,0;conj(Ds),0,0,0];

               he=(t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3);
               hh=(t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3);
               Hbdg=blkdiag(he,-transpose(hh));

               gs=1;

               Hbdg=Hbdg+(Hds);

               [u,d]=eig(Hbdg);
               F=diag(1./(exp(diag(d)*beta)+1));
               correls(:,:,i,j)=(u*F*u');

               ds(i,j)=-NF*Vs*conj(gs)*correls(1,4,i,j);
               en1s(i,j)=d(1,1);en2s(i,j)=d(2,2);en3s(i,j)=d(3,3);en4s(i,j)=d(4,4);
               conds(i,j)=NF*Vs*correls(1,4,i,j)*correls(4,1,i,j);
            end
        end
        DDsout(s)=sum(sum(ds));
        EEtots(s)=sum(en1s(en1s<0))+sum(en2s(en2s<0))+sum(en3s(en3s<0))+sum(en4s(en4s<0));
        cconds(s)=sum(sum(conds));
    end
    DDDsout{a}=DDsout;
    EEEtots{a}=EEtots;
    [Nones,I]=min(abs(DDsout./(DDsinput'+epsilon)-1));
    EEtotsf(a)=EEtots(I);
    Indexs(a)=I;
    ccondsf(a)=cconds(I);

    for pp11=1:length(DDpp11input) % (p+11)-wave check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dpp11=DDpp11input(pp11);
        for i=1:length(kkx)
            for j=1:length(kky)
                kx=kkx(i);
                ky=kky(j);
                kxm=-kkx(i);
                kym=-kky(j);

                he=(t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3);
                hh=(t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3);
                Hbdg=blkdiag(he,-transpose(hh));
                
                ges=kx^2+ky^2;
                gpp=kx+1i*ky;
                gpm=kx-1i*ky;
                gd=kx^2-ky^2;

                Dpm11=0;    %D11=Dpp11*gpp solution
                Dpm22=0;    %D22=Dpp22*gpp solution
                Dpp22=Dpp11;    % p+ solution
                D11=Dpp11*gpp+Dpm11*gpm;
                D22=Dpp22*gpp+Dpm22*gpm;
                Hdp=Vp*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];

                Hbdg=Hbdg+(Hdp);

                [u,d]=eig(Hbdg);
                F=diag(1./(exp(diag(d)*beta)+1));
                correlpp(:,:,i,j)=(u*F*u');

                dpp11(i,j)=-NF*Vp*conj(gpp)*correlpp(1,3,i,j);
                dpp22(i,j)=-NF*Vp*conj(gpp)*correlpp(2,4,i,j);
                en1pp(i,j)=d(1,1);en2pp(i,j)=d(2,2);en3pp(i,j)=d(3,3);en4pp(i,j)=d(4,4);
                % condpp11(i,j)=NF*Vp*(conj(ges)*ges+conj(gpp)*gpp+conj(gpm)*gpm+conj(gd)*gd)*correl(1,3,i,j)*correl(3,1,i,j);
                condpp11(i,j)=NF*Vp*(conj(gpp)*gpp+conj(gpm)*gpm)*correlpp(1,3,i,j)*correlpp(3,1,i,j);
            end
        end
        DDpp11out(pp11)=sum(sum(dpp11));
        DDpp22out(pp11)=sum(sum(dpp22));
        EEtotpp11(pp11)=sum(en1pp(en1pp<0))+sum(en2pp(en2pp<0))+sum(en3pp(en3pp<0))+sum(en4pp(en4pp<0));
        ccondpp11(pp11)=sum(sum(condpp11));
    end
    DDDpp11out{a}=DDpp11out;
    DDDpp22out{a}=DDpp22out;
    EEEtotpp11{a}=EEtotpp11;
    [Nonepp11,I]=min(abs(DDpp11out./(DDpp11input'+epsilon)-1));
    EEtotpp11f(a)=EEtotpp11(I);
    Indexpp11(a)=I;
    DDpp11out=zeros(width_numpp11,1);
    DDpp22out=zeros(width_numpp22,1);
    EEtotpp11=zeros(width_numpp11,1);
    parboxpp{a}=[Vp;DDpp11input(I)];
    ccondpp11f(a)=ccondpp11(I);
    for pm11=1:length(DDpm11input) % (p-11)-wave check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dpm11=DDpm11input(pm11);
        for i=1:length(kkx)
            for j=1:length(kky)
                kx=kkx(i);
                ky=kky(j);
                kxm=-kkx(i);
                kym=-kky(j);

                he=(t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3);
                hh=(t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3);
                Hbdg=blkdiag(he,-transpose(hh));
                
                ges=kx^2+ky^2;
                gpp=kx+1i*ky;
                gpm=kx-1i*ky;
                gd=kx^2-ky^2;

                Dpp11=0;    %D11=Dpm11*gpm solution
                Dpp22=0;    %D22=Dpm22*gpm solution
                Dpm22=Dpm11;    % p+ solution
                D11=Dpp11*gpp+Dpm11*gpm;
                D22=Dpp22*gpp+Dpm22*gpm;
                Hdp=Vp*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];

                Hbdg=Hbdg+(Hdp);

                [u,d]=eig(Hbdg);
                F=diag(1./(exp(diag(d)*beta)+1));
                correlpm(:,:,i,j)=(u*F*u');

                dpm11(i,j)=-NF*Vp*conj(gpm)*correlpm(1,3,i,j);
                dpm22(i,j)=-NF*Vp*conj(gpm)*correlpm(2,4,i,j);
                en1pm(i,j)=d(1,1);en2pm(i,j)=d(2,2);en3pm(i,j)=d(3,3);en4pm(i,j)=d(4,4);
                % condpm11(i,j)=NF*Vp*(conj(ges)*ges+conj(gpp)*gpp+conj(gpm)*gpm+conj(gd)*gd)*correl(1,3,i,j)*correl(3,1,i,j);
                condpm11(i,j)=NF*Vp*(conj(gpp)*gpp+conj(gpm)*gpm)*correlpm(1,3,i,j)*correlpm(3,1,i,j);
            end
        end
        DDpm11out(pm11)=sum(sum(dpm11));
        DDpm22out(pm11)=sum(sum(dpm22));
        EEtotpm11(pm11)=sum(en1pm(en1pm<0))+sum(en2pm(en2pm<0))+sum(en3pm(en3pm<0))+sum(en4pm(en4pm<0));
        ccondpm11(pm11)=sum(sum(condpm11));
    end
    DDDpm11out{a}=DDpm11out;
    DDDpm22out{a}=DDpm22out;
    EEEtotpm11{a}=EEtotpm11;
    [Nonepm11,I]=min(abs(DDpm11out./(DDpm11input'+epsilon)-1));
    EEtotpm11f(a)=EEtotpm11(I);
    Indexpm11(a)=I;
    DDpm11out=zeros(width_numpm11,1);
    DDpm22out=zeros(width_numpm22,1);
    EEtotpm11=zeros(width_numpm11,1);
    parboxpm{a}=[Vp;DDpm11input(I)];
    ccondpm11f(a)=ccondpm11(I);

    for px11=1:length(DDpx11input) % (p_x)-wave check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dpp11=DDpx11input(px11);
        for i=1:length(kkx)
            for j=1:length(kky)
                kx=kkx(i);
                ky=kky(j);
                kxm=-kkx(i);
                kym=-kky(j);

                he=(t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3);
                hh=(t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3);
                Hbdg=blkdiag(he,-transpose(hh));
                
                ges=kx^2+ky^2;
                gpx=sqrt(2)*kx;
                gpy=sqrt(2)*ky;
                gpp=kx+1i*ky;
                gpm=kx-1i*ky;
                gd=kx^2-ky^2;

                Dpm11=Dpp11;    %D11=2*Dpp11*kx solution
                Dpp22=Dpp11;    % px solution
                Dpm22=Dpp22;    %D22=2*Dpp22*kx solution
                D11=Dpp11*gpx;%D11=Dpp11*gpp+Dpm11*gpm;
                D22=Dpp11*gpx;%D22=Dpp22*gpp+Dpm22*gpm;
                Hdp=Vp*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];

                Hbdg=Hbdg+(Hdp);

                [u,d]=eig(Hbdg);
                F=diag(1./(exp(diag(d)*beta)+1));
                correlpx(:,:,i,j)=(u*F*u');

                dpx11(i,j)=-NF*Vp*conj(gpx)*correlpx(1,3,i,j);
                dpx22(i,j)=-NF*Vp*conj(gpy)*correlpx(2,4,i,j);
                en1px(i,j)=d(1,1);en2px(i,j)=d(2,2);en3px(i,j)=d(3,3);en4px(i,j)=d(4,4);
                % condpx11(i,j)=NF*Vp*(conj(ges)*ges+conj(gpx)*gpx+conj(gpy)*gpy+conj(gd)*gd)*correlpx(1,3,i,j)*correlpx(3,1,i,j);
                condpx11(i,j)=NF*Vp*(conj(gpx)*gpx+conj(gpy)*gpy)*correlpx(1,3,i,j)*correlpx(3,1,i,j);
            end
        end
        DDpx11out(px11)=sum(sum(dpx11));
        DDpx22out(px11)=sum(sum(dpx22));
        EEtotpx11(px11)=sum(en1px(en1px<0))+sum(en2px(en2px<0))+sum(en3px(en3px<0))+sum(en4px(en4px<0));
        ccondpx11(px11)=sum(sum(condpx11));
    end
    DDDpx11out{a}=DDpx11out;
    DDDpx22out{a}=DDpx22out;
    EEEtotpx11{a}=EEtotpx11;
    [Nonepx11,I]=min(abs(DDpx11out./(DDpx11input'+epsilon)-1));
    EEtotpx11f(a)=EEtotpx11(I);
    Indexpx11(a)=I;
    DDpx11out=zeros(width_numpx11,1);
    DDpx22out=zeros(width_numpx22,1);
    EEtotpx11=zeros(width_numpx11,1);
    parboxpx{a}=[Vp;DDpx11input(I)];
    ccondpx11f(a)=ccondpx11(I);

    for py11=1:length(DDpy11input) % (p_y)-wave check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dpm11=DDpy11input(py11);
        for i=1:length(kkx)
            for j=1:length(kky)
                kx=kkx(i);
                ky=kky(j);
                kxm=-kkx(i);
                kym=-kky(j);

                he=(t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3);
                hh=(t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3);
                Hbdg=blkdiag(he,-transpose(hh));

                ges=kx^2+ky^2;
                gpx=sqrt(2)*kx;
                gpy=sqrt(2)*ky;
                gpp=kx+1i*ky;
                gpm=kx-1i*ky;
                gd=kx^2-ky^2;

                Dpp11=-Dpm11;    %D11=2*i*Dpp11*ky solution
                Dpm22=Dpm11;    % py solution
                Dpp22=-Dpm22;    %D22=2*i*Dpp22*ky solution
                D11=Dpm11*gpy;
                D22=Dpm11*gpy;
                Hdp=Vp*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];

                Hbdg=Hbdg+(Hdp);

                [u,d]=eig(Hbdg);
                F=diag(1./(exp(diag(d)*beta)+1));
                correlpy(:,:,i,j)=(u*F*u');

                dpy11(i,j)=-NF*Vp*conj(gpy)*correlpy(1,3,i,j);
                dpy22(i,j)=-NF*Vp*conj(gpy)*correlpy(2,4,i,j);
                en1py(i,j)=d(1,1);en2py(i,j)=d(2,2);en3py(i,j)=d(3,3);en4py(i,j)=d(4,4);
                % condpy11(i,j)=NF*Vp*(conj(ges)*ges+conj(gpx)*gpx+conj(gpy)*gpy+conj(gd)*gd)*correl(1,3,i,j)*correl(3,1,i,j);
                condpy11(i,j)=NF*Vp*(conj(gpx)*gpx+conj(gpy)*gpy)*correlpy(1,3,i,j)*correlpy(3,1,i,j);
            end
        end
        DDpy11out(py11)=sum(sum(dpy11));
        DDpy22out(py11)=sum(sum(dpy22));
        EEtotpy11(py11)=sum(en1py(en1py<0))+sum(en2py(en2py<0))+sum(en3py(en3py<0))+sum(en4py(en4py<0));
        ccondpy11(py11)=sum(sum(condpy11));
    end
    DDDpy11out{a}=DDpy11out;
    DDDpy22out{a}=DDpy22out;
    EEEtotpy11{a}=EEtotpy11;
    [Nonepy11,I]=min(abs(DDpy11out./(DDpy11input'+epsilon)-1));
    EEtotpy11f(a)=EEtotpy11(I);
    Indexpy11(a)=I;
    DDpy11out=zeros(width_numpy11,1);
    DDpy22out=zeros(width_numpy22,1);
    EEtotpy11=zeros(width_numpy11,1);
    parboxpy{a}=[Vp;DDpy11input(I)];
    ccondpy11f(a)=ccondpy11(I);
    for mix11=1:length(DDmix11input) % (mixed)-wave check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dppmix11=DDmix11input(mix11);
        for i=1:length(kkx)
            for j=1:length(kky)
                kx=kkx(i);
                ky=kky(j);
                kxm=-kkx(i);
                kym=-kky(j);

                he=(t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3);
                hh=(t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3);
                Hbdg=blkdiag(he,-transpose(hh));

                ges=kx^2+ky^2;
                gpp=kx+1i*ky;
                gpm=kx-1i*ky;
                gd=kx^2-ky^2;

                Dpmmix11=0;    %D11=p+ip solution
                Dpmmix22=Dppmix11;    % mixed solution
                Dppmix22=0;    %D22=p-ip solution
                Dmix11=Dppmix11*gpp+Dpmmix11*gpm;
                Dmix22=Dppmix22*gpp+Dpmmix22*gpm;
                Hdp=Vp*[0,0,Dmix11,0;0,0,0,Dmix22;conj(Dmix11),0,0,0;0,conj(Dmix22),0,0];

                Hbdg=Hbdg+(Hdp);

                [u,d]=eig(Hbdg);
                F=diag(1./(exp(diag(d)*beta)+1));
                correlmix(:,:,i,j)=(u*F*u');

                dmix11(i,j)=-NF*Vp*conj(gpp)*correlmix(1,3,i,j);
                dmix22(i,j)=-NF*Vp*conj(gpp)*correlmix(2,4,i,j);
                en1mix(i,j)=d(1,1);en2mix(i,j)=d(2,2);en3mix(i,j)=d(3,3);en4mix(i,j)=d(4,4);
                % condmix11(i,j)=NF*Vp*(conj(ges)*ges+conj(gpp)*gpp+conj(gpm)*gpm+conj(gd)*gd)*correl(1,3,i,j)*correl(3,1,i,j);     
                condmix11(i,j)=NF*Vp*(conj(gpp)*gpp+conj(gpm)*gpm)*correlmix(1,3,i,j)*correlmix(3,1,i,j);   
            end
        end
        DDmix11out(mix11)=sum(sum(dmix11));
        DDmix22out(mix11)=sum(sum(dmix22));
        EEtotmix11(mix11)=sum(en1mix(en1mix<0))+sum(en2mix(en2mix<0))+sum(en3mix(en3mix<0))+sum(en4mix(en4mix<0));
        ccondmix11(mix11)=sum(sum(condmix11));
    end
    DDDmix11out{a}=DDmix11out;
    DDDmix22out{a}=DDmix22out;
    EEEtotmix11{a}=EEtotmix11;
    [Nonemix11,I]=min(abs(DDmix11out./(DDmix11input'+epsilon)-1));
    EEtotmix11f(a)=EEtotmix11(I);
    Indexmix11(a)=I;
    DDmix11out=zeros(width_nummix11,1);
    DDmix22out=zeros(width_nummix22,1);
    EEtotmix11=zeros(width_nummix11,1);
    parboxmix{a}=[Vp;DDmix11input(I)];
    ccondmix11f(a)=ccondmix11(I);
end

%%%% self consistency checking plot %%%%
for v=1:param_width_num
    subplot(2,6,1);
    plot(DDsinput,real(DDDsout{v}))
    hold on
    title('s - wave \Delta_{out}')
    xlabel('\Delta^{s}_{input}')
    ylabel('\Delta^{s}_{output}')
end
grid on
legend(legendLabels)
for v=1:param_width_num
    subplot(2,6,2);
    plot(DDpp11input,real(DDDpp11out{v}))
    hold on
    title('p^{+}_{\uparrow\uparrow} - wave \Delta')
    xlabel('\Delta^{p+}_{\uparrow\uparrow,input}')
    ylabel('\Delta^{p+}_{\uparrow\uparrow,output}')
end
grid on
legend(legendLabelp)
for v=1:param_width_num
    subplot(2,6,3);
    plot(DDpm11input,real(DDDpm11out{v}))
    hold on
    title('p^{-}_{\uparrow\uparrow} - wave \Delta')
    xlabel('\Delta^{p-}_{\uparrow\uparrow,input}')
    ylabel('\Delta^{p-}_{\uparrow\uparrow,output}')
end
grid on
legend(legendLabelp)
for v=1:param_width_num
    subplot(2,6,4);
    plot(DDpx11input,real(DDDpx11out{v}))
    hold on
    title('p^{x}_{\uparrow\uparrow} - wave \Delta')
    xlabel('\Delta^{p^x}_{\uparrow\uparrow,input}')
    ylabel('\Delta^{p^x}_{\uparrow\uparrow,output}')
end
grid on
legend(legendLabelp)
for v=1:param_width_num
    subplot(2,6,5);
    plot(DDpy11input,real(DDDpy11out{v}))
    hold on
    title('p^{y}_{\uparrow\uparrow} - wave \Delta')
    xlabel('\Delta^{p^y}_{\uparrow\uparrow,input}')
    ylabel('\Delta^{p^y}_{\uparrow\uparrow,output}')
end
grid on
legend(legendLabelp)

%%%% self consistent orderparams dots %%%%
for v=1:param_width_num
    subplot(2,6,7);
    plot(DDsinput, real(DDDsout{v}-DDsinput'))
    hold on
    title('s - wave self-consistency check')
    xlabel('\Delta^{s}_{input}')
    ylabel('\Delta^{s}_{output}-\Delta_{input}')
end
grid on
for v=1:param_width_num
    subplot(2,6,7)
    plot(DDsinput(Indexs(v)),real(DDDsout{v}(Indexs(v))-DDsinput(Indexs(v))),'k*')
end

for v=1:param_width_num
    subplot(2,6,8);
    plot(DDpp11input, real(DDDpp11out{v}-DDpp11input'))
    hold on
    title('p^{+}_{\uparrow\uparrow} - wave self-consistency check')
    xlabel('\Delta^{p+}_{\uparrow\uparrow,input}')
    ylabel('\Delta^{p+}_{\uparrow\uparrow,output}-\Delta^{p+}_{\uparrow\uparrow,input}')
end
grid on
for v=1:param_width_num
    subplot(2,6,8)
    plot(DDpp11input(Indexpp11(v)),real(DDDpp11out{v}(Indexpp11(v))-DDpp11input(Indexpp11(v))),'k*')
end

for v=1:param_width_num
    subplot(2,6,9);
    plot(DDpm11input, real(DDDpm11out{v}-DDpm11input'))
    hold on
    title('p^{-}_{\uparrow\uparrow} - wave self-consistency check')
    xlabel('\Delta^{p-}_{\uparrow\uparrow,input}')
    ylabel('\Delta^{p-}_{\uparrow\uparrow,output}-\Delta^{p-}_{\uparrow\uparrow,input}')
end
grid on
for v=1:param_width_num
    subplot(2,6,9)
    plot(DDpm11input(Indexpm11(v)),real(DDDpm11out{v}(Indexpm11(v))-DDpm11input(Indexpm11(v))),'k*')
end

for v=1:param_width_num
    subplot(2,6,10);
    plot(DDpx11input, real(DDDpx11out{v}-DDpx11input'))
    hold on
    title('p^{x}_{\uparrow\uparrow} - wave self-consistency check')
    xlabel('\Delta^{p^x}_{\uparrow\uparrow,input}')
    ylabel('\Delta^{p^x}_{\uparrow\uparrow,output}-\Delta^{p+}_{\uparrow\uparrow,input}')
end
grid on
for v=1:param_width_num
    subplot(2,6,10)
    plot(DDpx11input(Indexpx11(v)),real(DDDpx11out{v}(Indexpx11(v))-DDpx11input(Indexpx11(v))),'k*')
end

for v=1:param_width_num
    subplot(2,6,11);
    plot(DDpy11input, real(DDDpy11out{v}-DDpy11input'))
    hold on
    title('p^{y}_{\uparrow\uparrow} - wave self-consistency check')
    xlabel('\Delta^{p^y}_{\uparrow\uparrow,input}')
    ylabel('\Delta^{p^y}_{\uparrow\uparrow,output}-\Delta^{p^y}_{\uparrow\uparrow,input}')
end
grid on
for v=1:param_width_num
    subplot(2,6,11)
    plot(DDpy11input(Indexpy11(v)),real(DDDpy11out{v}(Indexpy11(v))-DDpy11input(Indexpy11(v))),'k*')
end

%%%% mixed p wave plot
for v=1:param_width_num
    subplot(2,6,6);
    plot(DDmix11input,real(DDDmix11out{v}))
    hold on
    title('p^{mix}_{\uparrow\uparrow} - wave \Delta')
    xlabel('\Delta^{p^{mix}_{\uparrow\uparrow,input}}')
    ylabel('\Delta^{p^{mix}_{\uparrow\uparrow,output}}')
end
grid on
legend(legendLabelp)
for v=1:param_width_num
    subplot(2,6,12);
    plot(DDmix11input, real(DDDmix11out{v}-DDmix11input'))
    hold on
    title('p^{mix}_{\uparrow\uparrow} - wave self-consistency check')
    xlabel('\Delta^{p^{mix}_{\uparrow\uparrow,input}}')
    ylabel('\Delta^{p^{mix}_{\uparrow\uparrow,output}}-\Delta^{p^{mix}_{\uparrow\uparrow,input}}')
end
grid on
for v=1:param_width_num
    subplot(2,6,12)
    plot(DDmix11input(Indexmix11(v)),real(DDDmix11out{v}(Indexmix11(v))-DDmix11input(Indexmix11(v))),'k*')
end
%%%% Total occupied energy - Vs,Vp plot %%%%
figure
ax1=gca;
hold on
plot(VVS,squeeze(EEtotsf)/abs(EEtotsf(1)))
hold(ax1,'on')
plot(ax1,VVS,squeeze(EEtotpp11f)/abs(EEtotsf(1)),'o')
plot(ax1,VVS,squeeze(EEtotpm11f)/abs(EEtotsf(1)))
plot(ax1,VVS,squeeze(EEtotpx11f)/abs(EEtotsf(1)),'o')
plot(ax1,VVS,squeeze(EEtotpy11f)/abs(EEtotsf(1)))
plot(ax1,VVS,squeeze(EEtotmix11f)/abs(EEtotsf(1)), '*')
xlabel(ax1,'V_{s} (V_{p}='+string(T)+'-V_{s})')
ylabel(ax1,'Occupied Energy')
legend(ax1,{'s-wave', 'p^{+}_{\uparrow\uparrow}-wave','p^{-}_{\uparrow\uparrow}-wave','p^{x}_{\uparrow\uparrow}-wave','p^{y}_{\uparrow\uparrow}-wave','p^{mix}_{\uparrow\uparrow} - wave'})
xlim(ax1,[-9.2 -2.2]);xticks(ax1,round(VVS,2));yticks(ax1,[])

ax2=axes('Position', ax1.Position, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none');
xlabel(ax2, 'V_{p}');xlim(ax2,[-9.2 -2.2]);xticks(ax2,round(VVS,2));yticks(ax2,[])
xticklabels(ax2,num2cell(round(-11.4-VVS,2)))
%%%% condensation energy -Vs,Vp plot %%%%
% figure
% hold on
% plot(VVS,squeeze(ccondsf))
% hold on
% plot(VVS,squeeze(abs(ccondpp11f)),'o')
% plot(VVS,squeeze(abs(ccondpm11f)))
% plot(VVS,squeeze(abs(ccondpx11f)),'o')
% plot(VVS,squeeze(abs(ccondpy11f)))
% plot(VVS,squeeze(abs(ccondmix11f)), '*')
% xlabel('V_{s}, (V_{p}='+string(T)+'-V_{s})')
% ylabel('E_{tot}')
% legend('s-wave', 'p^{+}_{\uparrow\uparrow}-wave','p^{-}_{\uparrow\uparrow}-wave','p^{x}_{\uparrow\uparrow}-wave','p^{y}_{\uparrow\uparrow}-wave','p^{mix}_{\uparrow\uparrow} - wave')
%% 
% load('aRuO2AMSCdata.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Spontaneously broken TRS by phase difference between p11 and p22 %%%%
phi_width_num=100;
pphi=linspace(0,2*pi,phi_width_num);
%%%% SymmBreaker params %%%%
beta_width_num=100;
bbeta=linspace(0,1,beta_width_num);
%%%% occupied energy for self-consistent solution array %%%%
phiEEtotpp11f=zeros(phi_width_num,1);
phiEEtotpm11f=zeros(phi_width_num,1);
phiEEtotpx11f=zeros(phi_width_num,1);
phiEEtotpy11f=zeros(phi_width_num,1);
b=1;
%%%% p+ strained system %%%%
for z=1:param_width_num
    parvecpp=parboxpp{z};
    Vp=parvecpp(1);
    Dpp11=parvecpp(2);
    % for BETA=1:length(bbeta)
    %     b=bbeta(BETA);
    for s=1:length(pphi)
        phi=pphi(s);
        for i=1:length(kkx)
            for j=1:length(kky)
                kx=kkx(i);
                ky=kky(j);
                kxm=-kkx(i);
                kym=-kky(j);

                % % % % SymmBreaker0=m*pauli(1);
                % % % % SymmBreakerm0=m*pauli(1);
                SymmBreaker1=kx*pauli(1)+ky*pauli(2);
                SymmBreakerm1=kxm*pauli(1)+kym*pauli(2);
                SymmBreaker2=kx^2*pauli(1)+ky^2*pauli(2);
                SymmBreakerm2=kxm^2*pauli(1)+kym^2*pauli(2);
                % SymmBreaker2cross=kx^2*pauli(2)+ky^2*pauli(1);
                % SymmBreakerm2cross=kxm^2*pauli(2)+kym^2*pauli(1);
                SymmBreaker3xxxx=kx^3*pauli(2)+ky^3*pauli(1);
                SymmBreakerm3xxxx=kxm^3*pauli(2)+kym^3*pauli(1);
                % SymmBreaker3yyyy=ky^3*pauli(2);
                % SymmBreakerm3yyyy=kym^3*pauli(2);
                % SymmBreaker3xxyx=kx^2*ky*pauli(1);
                % SymmBreakerm3xxyx=kxm^2*kym*pauli(1);
                % SymmBreaker3xxyy=kx^2*ky*pauli(2);
                % SymmBreakerm3xxyy=kxm^2*kym*pauli(2);
                % SymmBreaker3xyyx=kx*ky^2*pauli(1);
                % SymmBreakerm3xyyx=kxm*kym^2*pauli(1);
                % SymmBreaker3xyyy=kx*ky^2*pauli(2);
                % SymmBreakerm3xyyy=kxm*kym^2*pauli(2);
                % SymmBreaker3yyyy=ky^3*pauli(2);
                % SymmBreakerm3yyyy=kym^3*pauli(2);

                he=(t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3)+b*SymmBreaker3xxxx;
                hh=(t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3)+b*SymmBreakerm3xxxx;
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
                correlppst(:,:,i,j)=(u*F*u');
                en1ppst(i,j)=d(1,1);en2ppst(i,j)=d(2,2);en3ppst(i,j)=d(3,3);en4ppst(i,j)=d(4,4);
            end
        end
        phiEEtotpp11f(s)=sum(en1ppst(en1ppst<0))+sum(en2ppst(en2ppst<0))+sum(en3ppst(en3ppst<0))+sum(en4ppst(en4ppst<0));
    end
    [Eminpp11,phiminpp11index]=min(phiEEtotpp11f);
    %     pphiminpp11(BETA)=pphi(phiminpp11index);
    % end
    phiEEEtotpp11{z}=phiEEtotpp11f;
end
pppppp=1

%%%% p- strained system
for z=1:param_width_num
    parvecpm=parboxpm{z};
    Vp=parvecpm(1);
    Dpm11=parvecpm(2);
    % for BETA=1:length(bbeta)
    %     b=bbeta(BETA);
    for p=1:length(pphi)
        phi=pphi(p);
        for i=1:length(kkx)
            for j=1:length(kky)
                kx=kkx(i);
                ky=kky(j);
                kxm=-kkx(i);
                kym=-kky(j);
                
                % % % % SymmBreaker0=m*pauli(1);
                % % % % SymmBreakerm0=m*pauli(1);
                SymmBreaker1=kx*pauli(1)+ky*pauli(2);
                SymmBreakerm1=kxm*pauli(1)+kym*pauli(2);
                SymmBreaker2=kx^2*pauli(1)+ky^2*pauli(2);
                SymmBreakerm2=kxm^2*pauli(1)+kym^2*pauli(2);
                % SymmBreaker2cross=kx^2*pauli(2)+ky^2*pauli(1);
                % SymmBreakerm2cross=kxm^2*pauli(2)+kym^2*pauli(1);
                SymmBreaker3xxxx=kx^3*pauli(2)+ky^3*pauli(1);
                SymmBreakerm3xxxx=kxm^3*pauli(2)+kym^3*pauli(1);
                % SymmBreaker3yyyy=ky^3*pauli(2);
                % SymmBreakerm3yyyy=kym^3*pauli(2);
                % SymmBreaker3xxyx=kx^2*ky*pauli(1);
                % SymmBreakerm3xxyx=kxm^2*kym*pauli(1);
                % SymmBreaker3xxyy=kx^2*ky*pauli(2);
                % SymmBreakerm3xxyy=kxm^2*kym*pauli(2);
                % SymmBreaker3xyyx=kx*ky^2*pauli(1);
                % SymmBreakerm3xyyx=kxm*kym^2*pauli(1);
                % SymmBreaker3xyyy=kx*ky^2*pauli(2);
                % SymmBreakerm3xyyy=kxm*kym^2*pauli(2);
                % SymmBreaker3yyyy=ky^3*pauli(2);
                % SymmBreakerm3yyyy=kym^3*pauli(2);

                he=(t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3)+b*SymmBreaker3xxxx;
                hh=(t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3)+b*SymmBreakerm3xxxx;
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
                correlpmst(:,:,i,j)=(u*F*u');

                en1pmst(i,j)=d(1,1);en2pmst(i,j)=d(2,2);en3pmst(i,j)=d(3,3);en4pmst(i,j)=d(4,4);
            end
        end
        phiEEtotpm11f(p)=sum(en1pmst(en1pmst<0))+sum(en2pmst(en2pmst<0))+sum(en3pmst(en3pmst<0))+sum(en4pmst(en4pmst<0));
    end
    [Eminpm11,phiminpm11index]=min(phiEEtotpm11f);
    % pphiminpm11(BETA)=pphi(phiminpm11index);
    % end
    phiEEEtotpm11{z}=phiEEtotpm11f;
end
pmmmmm=1

for z=1:param_width_num
    parvecpx=parboxpx{z};
    Vp=parvecpx(1);
    Dpp11=parvecpx(2);
%     % for BETA=1:length(bbeta)
%     %     b=bbeta(BETA);    
    for s=1:length(pphi) %%%% p_x topologically TRS broken(pauli(1)) strained system
        phi=pphi(s);
        for i=1:length(kkx)
            for j=1:length(kky)
                kx=kkx(i);
                ky=kky(j);
                kxm=-kkx(i);
                kym=-kky(j);

                % % % % SymmBreaker0=m*pauli(1);
                % % % % SymmBreakerm0=m*pauli(1);
                SymmBreaker1=kx*pauli(1)+ky*pauli(2);
                SymmBreakerm1=kxm*pauli(1)+kym*pauli(2);
                SymmBreaker2=kx^2*pauli(1)+ky^2*pauli(2);
                SymmBreakerm2=kxm^2*pauli(1)+kym^2*pauli(2);
                % SymmBreaker2cross=kx^2*pauli(2)+ky^2*pauli(1);
                % SymmBreakerm2cross=kxm^2*pauli(2)+kym^2*pauli(1);
                SymmBreaker3xxxx=kx^3*pauli(2)+ky^3*pauli(1);
                SymmBreakerm3xxxx=kxm^3*pauli(2)+kym^3*pauli(1);
                % SymmBreaker3yyyy=ky^3*pauli(2);
                % SymmBreakerm3yyyy=kym^3*pauli(2);
                % SymmBreaker3xxyx=kx^2*ky*pauli(1);
                % SymmBreakerm3xxyx=kxm^2*kym*pauli(1);
                % SymmBreaker3xxyy=kx^2*ky*pauli(2);
                % SymmBreakerm3xxyy=kxm^2*kym*pauli(2);
                % SymmBreaker3xyyx=kx*ky^2*pauli(1);
                % SymmBreakerm3xyyx=kxm*kym^2*pauli(1);
                % SymmBreaker3xyyy=kx*ky^2*pauli(2);
                % SymmBreakerm3xyyy=kxm*kym^2*pauli(2);
                % SymmBreaker3yyyy=ky^3*pauli(2);
                % SymmBreakerm3yyyy=kym^3*pauli(2);

                he=(t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3)+b*SymmBreaker3xxxx;
                hh=(t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3)+b*SymmBreakerm3xxxx;
                Hbdg=blkdiag(he,-transpose(hh));

                gpp=kx+1i*ky;
                gpm=kx-1i*ky;

                Dpm11=Dpp11;    %D11=2*Dpp11*kx solution
                Dpp22=Dpp11;    % px solution
                Dpm22=Dpp22;    %D22=2*Dpp22*kx solution
                D11=Dpp11*gpp+Dpm11*gpm;
                D22=(Dpp22*gpp+Dpm22*gpm)*exp(1i*phi); %but different phase between wavefunctions
                Hdp=Vp*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];

                Hbdg=Hbdg+(Hdp);

                [u,d]=eig(Hbdg);
                F=diag(1./(exp(diag(d)*beta)+1));
                correlpxst(:,:,i,j)=(u*F*u');
                en1pxst(i,j)=d(1,1);en2pxst(i,j)=d(2,2);en3pxst(i,j)=d(3,3);en4pxst(i,j)=d(4,4);
            end
        end
        phiEEtotpx11f(s)=sum(en1pxst(en1pxst<0))+sum(en2pxst(en2pxst<0))+sum(en3pxst(en3pxst<0))+sum(en4pxst(en4pxst<0));
    end
    [Eminpx11,phiminpx11index]=min(phiEEtotpx11f);
    % pphiminpx11(BETA)=pphi(phiminpx11index);
    % end
    phiEEEtotpx11{z}=phiEEtotpx11f;    
end
pxxxxx=1

for z=1:param_width_num
    parvecpy=parboxpy{z};
    Vp=parvecpy(1);
    Dpm11=parvecpy(2);
%     % for BETA=1:length(bbeta)
%     %     b=bbeta(BETA);
    for p=1:length(pphi) % p_y topologically TRS broken(pauli(1)) strained system
        phi=pphi(p);
        for i=1:length(kkx)
            for j=1:length(kky)
                kx=kkx(i);
                ky=kky(j);
                kxm=-kkx(i);
                kym=-kky(j);

                 % % % % SymmBreaker0=m*pauli(1);
                % % % % SymmBreakerm0=m*pauli(1);
                SymmBreaker1=kx*pauli(1)+ky*pauli(2);
                SymmBreakerm1=kxm*pauli(1)+kym*pauli(2);
                SymmBreaker2=kx^2*pauli(1)+ky^2*pauli(2);
                SymmBreakerm2=kxm^2*pauli(1)+kym^2*pauli(2);
                % SymmBreaker2cross=kx^2*pauli(2)+ky^2*pauli(1);
                % SymmBreakerm2cross=kxm^2*pauli(2)+kym^2*pauli(1);
                SymmBreaker3xxxx=kx^3*pauli(2)+ky^3*pauli(1);
                SymmBreakerm3xxxx=kxm^3*pauli(2)+kym^3*pauli(1);
                % SymmBreaker3yyyy=ky^3*pauli(2);
                % SymmBreakerm3yyyy=kym^3*pauli(2);
                % SymmBreaker3xxyx=kx^2*ky*pauli(1);
                % SymmBreakerm3xxyx=kxm^2*kym*pauli(1);
                % SymmBreaker3xxyy=kx^2*ky*pauli(2);
                % SymmBreakerm3xxyy=kxm^2*kym*pauli(2);
                % SymmBreaker3xyyx=kx*ky^2*pauli(1);
                % SymmBreakerm3xyyx=kxm*kym^2*pauli(1);
                % SymmBreaker3xyyy=kx*ky^2*pauli(2);
                % SymmBreakerm3xyyy=kxm*kym^2*pauli(2);
                % SymmBreaker3yyyy=ky^3*pauli(2);
                % SymmBreakerm3yyyy=kym^3*pauli(2);

                he=(t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3)+b*SymmBreaker3xxxx;
                hh=(t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3)+b*SymmBreakerm3xxxx;
                Hbdg=blkdiag(he,-transpose(hh));

                gpp=kx+1i*ky;
                gpm=kx-1i*ky;

                Dpp11=-Dpm11;    %D11=2*i*Dpp11*ky solution
                Dpm22=Dpm11;    % py solution
                Dpp22=-Dpm22;    %D22=2*i*Dpp22*ky solution
                D11=(Dpp11*gpp+Dpm11*gpm)*1i;
                D22=(Dpp22*gpp+Dpm22*gpm)*exp(1i*phi)*1i;
                Hdp=Vp*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];

                Hbdg=Hbdg+(Hdp);

                [u,d]=eig(Hbdg);
                F=diag(1./(exp(diag(d)*beta)+1));
                correlpyst(:,:,i,j)=(u*F*u');

                en1pyst(i,j)=d(1,1);en2pyst(i,j)=d(2,2);en3pyst(i,j)=d(3,3);en4pyst(i,j)=d(4,4);
            end
        end
        phiEEtotpy11f(p)=sum(en1pyst(en1pyst<0))+sum(en2pyst(en2pyst<0))+sum(en3pyst(en3pyst<0))+sum(en4pyst(en4pyst<0));
    end
        [Eminpy11,phiminpy11index]=min(phiEEtotpy11f);
    %     pphiminpy11(BETA)=pphi(phiminpy11index);
    % end
    phiEEEtotpy11{z}=phiEEtotpy11f;    
end
pyyyyy=1

for z=1:param_width_num
    parvecmix=parboxmix{z};
    Vp=parvecmix(1);
    Dpp11=parvecmix(2);
%     % for BETA=1:length(bbeta)
%     %     b=bbeta(BETA);
    for s=1:length(pphi) %%%% mixed topologically TRS broken(pauli(1)) strained system
        phi=pphi(s);
        for i=1:length(kkx)
            for j=1:length(kky)
                kx=kkx(i);
                ky=kky(j);
                kxm=-kkx(i);
                kym=-kky(j);

                % % % % SymmBreaker0=m*pauli(1);
                % % % % SymmBreakerm0=m*pauli(1);
                SymmBreaker1=kx*pauli(1)+ky*pauli(2);
                SymmBreakerm1=kxm*pauli(1)+kym*pauli(2);
                SymmBreaker2=kx^2*pauli(1)+ky^2*pauli(2);
                SymmBreakerm2=kxm^2*pauli(1)+kym^2*pauli(2);
                % SymmBreaker2cross=kx^2*pauli(2)+ky^2*pauli(1);
                % SymmBreakerm2cross=kxm^2*pauli(2)+kym^2*pauli(1);
                SymmBreaker3xxxx=kx^3*pauli(2)+ky^3*pauli(1);
                SymmBreakerm3xxxx=kxm^3*pauli(2)+kym^3*pauli(1);
                % SymmBreaker3yyyy=ky^3*pauli(2);
                % SymmBreakerm3yyyy=kym^3*pauli(2);
                % SymmBreaker3xxyx=kx^2*ky*pauli(1);
                % SymmBreakerm3xxyx=kxm^2*kym*pauli(1);
                % SymmBreaker3xxyy=kx^2*ky*pauli(2);
                % SymmBreakerm3xxyy=kxm^2*kym*pauli(2);
                % SymmBreaker3xyyx=kx*ky^2*pauli(1);
                % SymmBreakerm3xyyx=kxm*kym^2*pauli(1);
                % SymmBreaker3xyyy=kx*ky^2*pauli(2);
                % SymmBreakerm3xyyy=kxm*kym^2*pauli(2);
                % SymmBreaker3yyyy=ky^3*pauli(2);
                % SymmBreakerm3yyyy=kym^3*pauli(2);

                he=(t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3)+b*SymmBreaker3xxxx;
                hh=(t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3)+b*SymmBreakerm3xxxx;
                Hbdg=blkdiag(he,-transpose(hh));

                gpp=kx+1i*ky;
                gpm=kx-1i*ky;

                Dpm11=Dpp11;    %D11=2*Dpp11*kx solution
                Dpp22=Dpp11;    % px solution
                Dpm22=Dpp22;    %D22=2*Dpp22*kx solution
                D11=Dpp11*gpp+Dpm11*gpm;
                D22=(Dpp22*gpp+Dpm22*gpm)*exp(1i*phi); %but different phase between wavefunctions
                Hdp=Vp*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];

                Hbdg=Hbdg+(Hdp);

                [u,d]=eig(Hbdg);
                F=diag(1./(exp(diag(d)*beta)+1));
                correlmixst(:,:,i,j)=(u*F*u');
                en1mixst(i,j)=d(1,1);en2mixst(i,j)=d(2,2);en3mixst(i,j)=d(3,3);en4mixst(i,j)=d(4,4);
            end
        end
        phiEEtotmix11f(s)=sum(en1mixst(en1mixst<0))+sum(en2mixst(en2mixst<0))+sum(en3mixst(en3mixst<0))+sum(en4mixst(en4mixst<0));
    end
        [Eminmix11,phiminmix11index]=min(phiEEtotmix11f);
        % pphiminmix11(BETA)=pphi(phiminmix11index);
    % end
    phiEEEtotmix11{z}=phiEEtotmix11f;
end
pmixmix=1

%%%% Total occupied energy - phi plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(squeeze(pphi),squeeze(phiEEtotpp11f*NF))
hold on
plot(squeeze(pphi),squeeze(phiEEtotpm11f*NF), 'o')
xlabel('\phi');ylabel('Occupied Energy');xticks([0 pi/2 pi 3*pi/2 2*pi]);xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});xlim([0 2*pi]);yticklabels({});
legend('p^{+} C_4𝒯 broken occupied energy','p^{-} C_4𝒯 broken occupied energy')

figure
plot(squeeze(pphi),squeeze(phiEEtotpx11f*NF))
hold on
plot(squeeze(pphi),squeeze(phiEEtotpy11f*NF), 'o')
xlabel('\phi');ylabel('Occupied Energy');xticks([0 pi/2 pi 3*pi/2 2*pi]);xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});yticklabels({});xlim([0 2*pi])
legend('p^{x} C_4𝒯 broken occupied energy','p^{y} C_4𝒯 broken occupied energy')

figure
plot(squeeze(pphi),squeeze(phiEEtotpx11f*NF))
xlabel('\phi');ylabel('Occupied Energy');xticks([0 pi/2 pi 3*pi/2 2*pi]);xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});yticklabels({});xlim([0 2*pi])
legend('mixed state C_4𝒯 broken occupied energy')

%%%% beta vs phimin phase diagram %%%%
% plot(squeeze(bbeta),squeeze(pphiminpp11))
% figure
% plot(squeeze(bbeta),squeeze(pphiminpm11))
% figure
% plot(squeeze(bbeta),squeeze(pphiminpx11))
% figure
% plot(squeeze(bbeta),squeeze(pphiminpy11))
% figure
% plot(squeeze(bbeta),squeeze(pphiminmix11))