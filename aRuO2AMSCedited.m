clear all
%% parameters
t=1;
J=0.4;
mu=-1.8;
m=0;
beta=5000;
NOGx=10; % Number Of Grid x
NOGy=10; % Number Of Grid y
NF=1/(NOGx*NOGy); % Normalization Factor
kkx=linspace(-pi,pi,NOGx+1);
kkx(end)=[];
kky=linspace(-pi,pi,NOGy+1);
kky(end)=[];
%%%% input orderparams limit %%%%
% DDslim=5;
% DDpp11lim=15;
% DDpm11lim=15;
% DDpx11lim=15;
% DDpy11lim=15;
% DDmix11lim=15;

DDslim=5;
DDpp11lim=5;
DDpm11lim=5;
DDpx11lim=5;
DDpy11lim=5;
DDmix11lim=5;
%%%% number of widths of input orderparams %%%%
slice=200;
width_num=DDslim*slice;
width_numpp11=DDpp11lim*slice;
width_numpm11=DDpm11lim*slice;
width_numpx11=DDpx11lim*slice;
width_numpy11=DDpy11lim*slice;
width_nummix11=DDmix11lim*slice;
%%%% input orderparams %%%%
DDsinput=linspace(0,DDslim,width_num);
% DDsinput=linspace(0,0.281,20);
DDpp11input=linspace(0,DDpp11lim,width_numpp11);
DDpm11input=linspace(0,DDpm11lim,width_numpm11);
DDpx11input=linspace(0,DDpx11lim,width_numpx11);
DDpy11input=linspace(0,DDpy11lim,width_numpy11);
DDmix11input=linspace(0,DDmix11lim,width_nummix11);
%%%% zero arrays prelocation %%%%
%%%% output orderparams receivers %%%%
DDsout=zeros(width_num,1);
DDpp11out=zeros(width_numpp11,1);
DDpm11out=zeros(width_numpp11,1);
DDpx11out=zeros(width_numpx11,1);
DDpy11out=zeros(width_numpy11,1);
DDmix11out=zeros(width_nummix11,1);
%%%% occupied energy receivers for each output orderparams %%%%
EEtots=zeros(width_num,1);
EEtotpp11=zeros(width_numpp11,1);
EEtotpm11=zeros(width_numpp11,1);
EEtotpx11=zeros(width_numpx11,1);
EEtotpy11=zeros(width_numpy11,1);
EEtotmix11=zeros(width_nummix11,1);
%%%% output order params 
Dsout=0;
Dpp11out=0;
Dpm11out=0;
Dpx11out=0;
Dpy11out=0;
Dmix11out=0;
%%% number of trial interaction params %%%%
param_width_num=1;
% L=-9.2;
L=-2;
% R=-2.2;
R=0;
T=L+R;
% VVS=linspace(0,-3,param_width_num);
% VVP=flip(linspace(0,-3,param_width_num));
VVS=linspace(0,-2,11);
VVP=linspace(-2,0,11);
%%%% self-consistent occupied energy receivers %%%%
EEtotsf=zeros(param_width_num,1); 
EEtotpp11f=zeros(param_width_num,1);
EEtotpm11f=zeros(param_width_num,1);
EEtotpx11f=zeros(param_width_num,1);
EEtotpy11f=zeros(param_width_num,1);
EEtotmix11f=zeros(param_width_num,1);
%%%% self-consistent input orderparam indexes receivers %%%%
Indexs=zeros(param_width_num,1);
Indexpp11=zeros(param_width_num,1);
Indexpm11=zeros(param_width_num,1);
Indexpx11=zeros(param_width_num,1);
Indexpy11=zeros(param_width_num,1);
Indexmix11=zeros(param_width_num,1);
%%% cell arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DDDsout=cell(1,param_width_num);
EEEtots=cell(1,param_width_num);
DDDpp11out=cell(1,param_width_num);
EEEtotpp11=cell(1,param_width_num);
DDDpm11out=cell(1,param_width_num);
EEEtotpm11=cell(1,param_width_num);
DDDpx11out=cell(1,param_width_num);
EEEtotpx11=cell(1,param_width_num);
DDDpy11out=cell(1,param_width_num);
EEEtotpy11=cell(1,param_width_num);
DDDmix11out=cell(1,param_width_num);
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
%% gap_eq_solver
q=1;
for a=1:length(VVS)
    Vs=VVS(a);
    Vs
    % Vp=-11.4-Vs;
    % Vp=T-Vs;
    Vp=VVP(a);
    Vp

    legendLabels{a}=sprintf('V_{s}=%.2f', Vs);
    legendLabelp{a}=sprintf('V_{p}=%.2f', Vp);

    [DDsout,EEtots,cconds,ccorrels,en1s,en2s,en3s,en4s,ddeltas]=gap_eq_solver('s',NF,t,mu,J,beta,Vs,kkx,kky,DDsinput); %s-wave gap_eq_solver
    DDDsout{a}=DDsout;
    EEEtots{a}=EEtots;
    [Nones,I]=min(abs(DDsout./(DDsinput+epsilon)-1));
    if (DDsout(2)-DDsout(1))/(DDsinput(2)-DDsinput(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotsf(a)=EEtots(I);
    Indexs(a)=I;
    ccondsf(a)=cconds(I);
    correls{a}=ccorrels{I};
    ASDFs{a}=ddeltas;
    % deltas(a,:,:)=ddeltas{I};

    [DDpp11out,EEtotpp11,ccondpp,ccorrelpp,en1pp,en2pp,en3pp,en4pp,ddeltapp]=gap_eq_solver('pp11',NF,t,mu,J,beta,Vp,kkx,kky,DDpp11input);  %pp-wave gap_eq_solver
    DDDpp11out{a}=DDpp11out;
    EEEtotpp11{a}=EEtotpp11;
    [Nonepp11,I]=min(abs(DDpp11out./(DDpp11input+epsilon)-1));
    if (DDpp11out(2)-DDpp11out(1))/(DDpp11input(2)-DDpp11input(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotpp11f(a)=EEtotpp11(I);
    Indexpp11(a)=I;
    DDpp11out=zeros(width_numpp11,1);
    EEtotpp11=zeros(width_numpp11,1);
    parboxpp{a}=[Vp;DDpp11input(I)];
    ccondpp11f(a)=ccondpp(I);
    % correlpp(a,:,:,:,:)=ccorrelpp{I};
    ASDFpp{a}=ddeltapp;
    % deltapp(a,:,:)=ddeltapp{I};
    
    [DDpm11out,EEtotpm11,ccondpm,ccorrelpm,en1pm,en2pm,en3pm,en4pm,ddeltapm]=gap_eq_solver('pm11',NF,t,mu,J,beta,Vp,kkx,kky,DDpm11input); %pm-wave gap_eq_solver
    DDDpm11out{a}=DDpm11out;
    EEEtotpm11{a}=EEtotpm11;
    [Nonepm11,I]=min(abs(DDpm11out./(DDpm11input+epsilon)-1));
    if (DDpm11out(2)-DDpm11out(1))/(DDpm11input(2)-DDpm11input(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotpm11f(a)=EEtotpm11(I);
    Indexpm11(a)=I;
    DDpm11out=zeros(width_numpm11,1);
    EEtotpm11=zeros(width_numpm11,1);
    parboxpm{a}=[Vp;DDpm11input(I)];
    ccondpm11f(a)=ccondpm(I);
    % correlpm(a,:,:,:,:)=ccorrelpm{I};
    % ASDFpm{a}=ddeltapm;
    % deltapm(a,:,:)=ddeltapm{I};

    [DDpx11out,EEtotpx11,ccondpx,ccorrelpx,en1px,en2px,en3px,en4px,ddeltapx]=gap_eq_solver('px11',NF,t,mu,J,beta,Vp,kkx,kky,DDpx11input); %px-wave gap_eq_solver
    DDDpx11out{a}=DDpx11out;
    EEEtotpx11{a}=EEtotpx11;
    [Nonepx11,I]=min(abs(DDpx11out./(DDpx11input+epsilon)-1));
    if (DDpx11out(2)-DDpx11out(1))/(DDpx11input(2)-DDpx11input(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotpx11f(a)=EEtotpx11(I);
    Indexpx11(a)=I;
    DDpx11out=zeros(width_numpx11,1);
    EEtotpx11=zeros(width_numpx11,1);
    parboxpx{a}=[Vp;DDpx11input(I)];
    ccondpx11f(a)=ccondpx(I);
    % correlpx(a,:,:,:,:)=ccorrelpx{I};
    ASDFpx{a}=ddeltapx;
    % deltapx(a,:,:)=ddeltapx{I};
    
    [DDpy11out,EEtotpy11,ccondpy,ccorrelpy,en1py,en2py,en3py,en4py,ddeltapy]=gap_eq_solver('py11',NF,t,mu,J,beta,Vp,kkx,kky,DDpy11input); %py-wave gap_eq_solver
    DDDpy11out{a}=DDpy11out;
    EEEtotpy11{a}=EEtotpy11;
    [Nonepy11,I]=min(abs(DDpy11out./(DDpy11input+epsilon)-1));
    if (DDpy11out(2)-DDpy11out(1))/(DDpy11input(2)-DDpy11input(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotpy11f(a)=EEtotpy11(I);
    Indexpy11(a)=I;
    DDpy11out=zeros(width_numpy11,1);
    EEtotpy11=zeros(width_numpy11,1);
    parboxpy{a}=[Vp;DDpy11input(I)];
    ccondpy11f(a)=ccondpy(I);
    % correlpy(a,:,:,:,:)=ccorrelpy{I};
    % deltapy(a,:,:)=ddeltapy{I};

    
    [DDmix11out,EEtotmix11,ccondmix,ccorrelmix,en1mix,en2mix,en3mix,en4mix,ddeltamix]=gap_eq_solver('mix11',NF,t,mu,J,beta,Vp,kkx,kky,DDmix11input); %pmix-wave gap_eq_solver
    DDDmix11out{a}=DDmix11out;
    EEEtotmix11{a}=EEtotmix11;
    [Nonemix11,I]=min(abs(DDmix11out./(DDmix11input+epsilon)-1));
    if (DDmix11out(2)-DDmix11out(1))/(DDmix11input(2)-DDmix11input(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotmix11f(a)=EEtotmix11(I);
    Indexmix11(a)=I;
    DDmix11out=zeros(width_nummix11,1);
    EEtotmix11=zeros(width_nummix11,1);
    parboxmix{a}=[Vp;DDmix11input(I)];
    ccondmix11f(a)=ccondmix(I);
    % correlmix(a,:,:,:,:)=ccorrelmix{I};
    % deltamix(a,:,:)=ddeltamix{I};

    for k=1:length(DDsinput)
        kkxq=linspace(-pi,pi,NOGx+1);kkxq(end)=[];
        kkyq=linspace(-pi,pi,NOGx+1);kkyq(end)=[];
        qqx=linspace(-pi,pi,NOGx+1);qqx(end)=[];
        qqy=linspace(-pi,pi,NOGy+1);qqy(end)=[];
        for i=1:length(qqx)
            % for j=1:length(qqy)
            [correlf,en1f,en2f,en3f,en4f,deltaf]=Hbdg_solver('s',NF,t,mu,J,beta,Vs,kkx,kky,DDsinput(k),true,qqx(i),0);
            EEf(i)=sum(en1f(en1f<0))+sum(en2f(en2f<0))+sum(en3f(en3f<0))+sum(en4f(en4f<0));
            % end
        end
        Ef=fix((EEf*10^(-3))*10^(13))/10^13;
        [xx, yy]=find(Ef==min(Ef));
        qqq=qqx(yy(1));
        qqqq(a,k)=qqq;
        [correlff,en1ff,en2ff,en3ff,en4ff,deltaff]=Hbdg_solver('s',NF,t,mu,J,beta,Vs,kkxq,kkyq,DDsinput(k),true,qqq,0);
        DDffout(k)=sum(sum(deltaff));
        EEtotff(k)=sum(en1ff(en1ff<0))+sum(en2ff(en2ff<0))+sum(en3ff(en3ff<0))+sum(en4ff(en4ff<0));
        ccorrelff{k}=correlff;
    end
    DDDffout{a}=DDffout;
    EEEtotff{a}=EEtotff;
    [Noneff,I]=min(abs(DDffout./(DDsinput+epsilon)-1));
    if (DDffout(2)-DDffout(1))/(DDsinput(2)-DDsinput(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotfff(a)=EEtotff(I);
    ASDFff{a}=deltaff;
    Indexff(a)=I;
    % correlff(a,:,:,:,:)=ccorrelff{I};
    % clear ccorrelff
end


%% excluding I=2(debugging process)
for a=1:length(VVS)
    DDsout=DDDsout{a};    
    [Nones,I]=min(abs(DDsout./(DDsinput+epsilon)-1));
    if (DDsout(2)-DDsout(1))/(DDsinput(2)-DDsinput(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtots=EEEtots{a};
    EEtotsf(a)=EEtots(I);
    Indexs(a)=I;

    DDpp11out=DDDpp11out{a};
    [Nonepp11,I]=min(abs(DDpp11out./(DDpp11input+epsilon)-1));
    if (DDpp11out(2)-DDpp11out(1))/(DDpp11input(2)-DDpp11input(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotpp11=EEEtotpp11{a};
    EEtotpp11f(a)=EEtotpp11(I);
    Indexpp11(a)=I;

    DDpm11out=DDDpm11out{a};
    [Nonepm11,I]=min(abs(DDpm11out./(DDpm11input+epsilon)-1));
    if (DDpm11out(2)-DDpm11out(1))/(DDpm11input(2)-DDpm11input(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotpm11=EEEtotpm11{a};
    EEtotpm11f(a)=EEtotpm11(I);
    Indexpm11(a)=I;

    DDmix11out=DDDmix11out{a};    
    [Nonemix11,I]=min(abs(DDmix11out./(DDmix11input+epsilon)-1));
    if (DDmix11out(2)-DDmix11out(1))/(DDmix11input(2)-DDmix11input(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotmix11=EEEtotmix11{a};
    EEtotmix11f(a)=EEtotmix11(I);
    Indexmix11(a)=I;

    DDpx11out=DDDpx11out{a};    
    [Nonepx11,I]=min(abs(DDpx11out./(DDpx11input+epsilon)-1));
    if (DDpx11out(2)-DDpx11out(1))/(DDpx11input(2)-DDpx11input(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotpx11=EEEtotpx11{a};
    EEtotpx11f(a)=EEtotpx11(I);
    Indexpx11(a)=I;

    DDpy11out=DDDpy11out{a};
    [Nonepy11,I]=min(abs(DDpy11out./(DDpy11input+epsilon)-1));
    if (DDpy11out(2)-DDpy11out(1))/(DDpy11input(2)-DDpy11input(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotpy11=EEEtotpy11{a};
    EEtotpy11f(a)=EEtotpy11(I);
    Indexpy11(a)=I;

    DDffout=DDDffout{a};    
    [Noneff,I]=min(abs(DDffout./(DDsinput+epsilon)-1));
    if (DDffout(2)-DDffout(1))/(DDsinput(2)-DDsinput(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotff=EEEtotff{a};
    EEtotfff(a)=EEtotff(I);
    Indexff(a)=I;
end
%% self consistency checking plot %%%%
figure
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
% for v=1:param_width_num
%     subplot(2,6,3);
%     plot(DDpm11input,real(DDDpm11out{v}))
%     hold on
%     title('p^{-}_{\uparrow\uparrow} - wave \Delta')
%     xlabel('\Delta^{p-}_{\uparrow\uparrow,input}')
%     ylabel('\Delta^{p-}_{\uparrow\uparrow,output}')
% end
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
% for v=1:param_width_num
%     subplot(2,6,5);
%     plot(DDpy11input,real(DDDpy11out{v}))
%     hold on
%     title('p^{y}_{\uparrow\uparrow} - wave \Delta')
%     xlabel('\Delta^{p^y}_{\uparrow\uparrow,input}')
%     ylabel('\Delta^{p^y}_{\uparrow\uparrow,output}')
% end
% grid on
legend(legendLabelp)
%% self consistent orderparams dots
for v=1:param_width_num
    subplot(2,6,7);
    plot(DDsinput, real(DDDsout{v}-DDsinput))
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
    plot(DDpp11input, real(DDDpp11out{v}-DDpp11input))
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

% for v=1:param_width_num
%     subplot(2,6,9);
%     plot(DDpm11input, real(DDDpm11out{v}-DDpm11input))
%     hold on
%     title('p^{-}_{\uparrow\uparrow} - wave self-consistency check')
%     xlabel('\Delta^{p-}_{\uparrow\uparrow,input}')
%     ylabel('\Delta^{p-}_{\uparrow\uparrow,output}-\Delta^{p-}_{\uparrow\uparrow,input}')
% end
% grid on
% for v=1:param_width_num
%     subplot(2,6,9)
%     plot(DDpm11input(Indexpm11(v)),real(DDDpm11out{v}(Indexpm11(v))-DDpm11input(Indexpm11(v))),'k*')
% end

for v=1:param_width_num
    subplot(2,6,10);
    plot(DDpx11input, real(DDDpx11out{v}-DDpx11input))
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

% for v=1:param_width_num
%     subplot(2,6,11);
%     plot(DDpy11input, real(DDDpy11out{v}-DDpy11input))
%     hold on
%     title('p^{y}_{\uparrow\uparrow} - wave self-consistency check')
%     xlabel('\Delta^{p^y}_{\uparrow\uparrow,input}')
%     ylabel('\Delta^{p^y}_{\uparrow\uparrow,output}-\Delta^{p^y}_{\uparrow\uparrow,input}')
% end
% grid on
% for v=1:param_width_num
%     subplot(2,6,11)
%     plot(DDpy11input(Indexpy11(v)),real(DDDpy11out{v}(Indexpy11(v))-DDpy11input(Indexpy11(v))),'k*')
% end

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
    plot(DDmix11input, real(DDDmix11out{v}-DDmix11input))
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

 %% FFLO state

figure
for v=1:param_width_num
    subplot(2,1,1);
    plot(DDsinput,real(DDDffout{v}))
    hold on
    title('FFLO state \Delta')
    xlabel('\Delta_{FFLO,input}')
    ylabel('\Delta_{FFLO,output}')
end
grid on
legend(legendLabels)

for v=1:param_width_num
    subplot(2,1,2);
    plot(DDsinput, real(DDDffout{v}-DDsinput))
    hold on
    title('FFLO self-consistency check')
    xlabel('\Delta_{FFLO,input}')
    ylabel('\Delta_{FFLO,output}-\Delta_{FFLO,input}')
end
grid on
for v=1:param_width_num
    subplot(2,1,2)
    plot(DDsinput(Indexff(v)),real(DDDffout{v}(Indexff(v))-DDsinput(Indexff(v))),'k*')
end
%% FFLO p?
% figure
% for v=1:param_width_num
%     subplot(2,1,1);
%     plot(DDpp11input,real(DDDffpout{v}))
%     hold on
%     title('FFLOp state \Delta')
%     xlabel('\Delta_{FFLOp,input}')
%     ylabel('\Delta_{FFLOp,output}')
% end
% grid on
% legend(legendLabels)
% 
% for v=1:param_width_num
%     subplot(2,1,2);
%     plot(DDpp11input, real(DDDffpout{v}-DDpp11input))
%     hold on
%     title('FFLOp self-consistency check')
%     xlabel('\Delta_{FFLOp,input}')
%     ylabel('\Delta_{FFLOp,output}-\Delta_{FFLOp,input}')
% end
% grid on
% for v=1:param_width_num
%     subplot(2,1,2)
%     plot(DDpp11input(Indexffp(v)),real(DDDffpout{v}(Indexffp(v))-DDpp11input(Indexffp(v))),'k*')
% end
%% finished..
figure
tt2=tiledlayout(1,1);
% tt1=tiledlayout(1,2,'TileSpacing','tight');
% tt2=tiledlayout(tt1,'flow','TileSpacing','Compact');
% tt3=tiledlayout(tt1,'flow','TileSpacing','Compact');
% nexttile(tt2);
ax2=axes(tt2);
% title('(a)','FontSize',12)
% ax2.TitleHorizontalAlignment = 'left'; 
xlabel(ax2, 'Nearest-neighbor interaction strength','FontName','Helvetica', 'FontSize',20)
ax2.XAxisLocation='top';
ax2.YAxisLocation='right';
% ax2.XColor='#E56F15';
% ax2.YColor='#E56F15';
% ax2.XColor='#D62D20';
ax2.YColor='k';
% xlim(ax2, [L R]);
% xticks(ax2, [-2 -1.5 -1 -0.5 0]);
yticks(ax2,[]);
% xticklabels(ax2, {'0', '-0.5', '-1.0', '-1.5', '-2'})
ax2.XAxis.FontSize=20;
ax2.XAxis.FontName='Helvetica';
% ax2.TickLength=[0.02,0.025];
ax2.XMinorTick='on';
ax2.TickDir='out';
axis square
ax2.LineWidth=1;

ax1=axes(tt2);
% plot(ax1, VVS, squeeze(EEtotsf)/abs(squeeze(EEtotsf(1))),'LineStyle','none','LineWidth',1,'Marker','o','MarkerSize',12,'MarkerEdgeColor','#0057E7','MarkerFaceColor','none')
% hold on
% plot(ax1, VVS, squeeze(EEtotfff)/abs(squeeze(EEtotfff(1))),'LineStyle','none','LineWidth',1,'Marker','^','MarkerSize',12,'MarkerEdgeColor','#0057E7','MarkerFaceColor','none')
hold on
plot(ax1, VVS, squeeze(EEtotsf),'LineStyle','none','LineWidth',1,'Marker','o','MarkerSize',12,'MarkerEdgeColor','#0057E7','MarkerFaceColor','none')
plot(ax1, VVS, squeeze(EEtotfff),'LineStyle','none','LineWidth',1,'Marker','^','MarkerSize',12,'MarkerEdgeColor','#0057E7','MarkerFaceColor','none')
plot(ax1, VVS, squeeze(EEtotpx11f),'LineStyle','none','LineWidth',1,'Marker','square','MarkerSize',12,'MarkerEdgeColor','#D62D20','MarkerFaceColor','none')
plot(ax1, VVS, squeeze(EEtotpp11f),'LineStyle','none','LineWidth',1, 'Marker','o', 'MarkerSize',12,'MarkerEdgeColor','#D62D20','MarkerFaceColor','none')
plot(ax1, VVS, squeeze(EEtotmix11f), 'LineStyle','none','LineWidth',1,'Marker','x','MarkerSize',12,'MarkerEdgeColor','#D62D20','MarkerFaceColor','none')


legend(ax1, {'\it{s}\rm{-wave} (q=0)','\it{s}\rm{-wave} (FFLO)','\it{p}\rm{_{x}-wave}', 'chiral \it{p}\rm{-wave}', 'helical \it{p}\rm{-wave}'},'Position',[.16 .15, .211 .2],'FontName','Helvetica', 'FontSize',8)
% 
% axs1=axes('Position',[.20 .4 .261 .25]);
% box off
% plot(axs1,VVS(1:4),squeeze(EEtotsf(1:4)/abs(squeeze(EEtotsf(end)))),'LineStyle','none','LineWidth',1,'Marker','o','MarkerSize',8,'MarkerEdgeColor','#0057E7','MarkerFaceColor','none')
% hold on
% plot(axs1,VVS(1:4), squeeze(EEtotfff(1:4))/abs(squeeze(EEtotfff(end))),'LineStyle','none','LineWidth',1,'Marker','^','MarkerSize',8,'MarkerEdgeColor','#0057E7','MarkerFaceColor','none')
% plot(axs1,VVS(1:4), squeeze(EEtotpx11f(1:4))/abs(squeeze(EEtotpx11f(1))), 'LineStyle','none','LineWidth',1,'Marker','square','MarkerSize',8,'MarkerEdgeColor','#D62D20','MarkerFaceColor','none')
% plot(axs1,VVS(1:4), squeeze(EEtotpp11f(1:4))/abs(squeeze(EEtotpp11f(1))),'LineStyle','none','LineWidth',1, 'Marker','o', 'MarkerSize',8,'MarkerEdgeColor','#D62D20','MarkerFaceColor','none')
% plot(axs1,VVS(1:4), squeeze(EEtotmix11f(1:4))/abs(squeeze(EEtotmix11f(1))), 'LineStyle','none','LineWidth',1,'Marker','x','MarkerSize',8,'MarkerEdgeColor','#D62D20','MarkerFaceColor','none')
% axs1.YAxisLocation='left';
% axs1.XMinorTick='on';
% axs1.YMinorTick='on';
% axs1.LineWidth=1;
% % yticklabels(axs1, {'-0.50', '-0.45', '-0.40'})
% 

ax1.Color='none';
% ax1.XColor='#0057E7';
xlabel(ax1, 'Onsite interaction strength','FontName','Helvetica', 'FontSize',20)
ylabel(ax1, 'Ground state energy','FontName','Helvetica', 'FontSize',20)
% xlim(ax1, [L R]);
% ylim(ax1, [-9 -0.75])
% xticks(ax1, [-2 -1.5 -1 -0.5 0]);
% yticks(ax1, [-9 -7 -5 -3 -1])
% xticklabels(ax1, {'-2', '-1.5', '-1.0', '-0.5', '0'})
% yticklabels(ax1, {'-9.0', '-7.0', '-5.0', '-3.0', '-1.0'})
axis square
ax1.LineWidth=1;
ax1.XAxis.FontSize=20;
ax1.XAxis.FontName='Helvetica';
ax1.YAxis.FontSize=20;
ax1.YAxis.FontName='Helvetica';
ax1.TickLength=[0.02,0.025];
ax1.XMinorTick='on';
ax1.YMinorTick='on';
ax1.Box='off';
ax2.Box='off';



% saveas(gcf, 'Vs vs Vp phase diagram AM2.png');

%% 

%subplot of Onsite interaction
% figure;
% ax4 = gca;
% hold on
% plot(VVS(1:6), squeeze(EEtotsf(1:6)*NF/abs(squeeze(EEtotsf(end)))),'LineStyle','none','LineWidth',1.5,'Marker','o','MarkerSize',8,'MarkerEdgeColor','#0057E7','MarkerFaceColor','none')
% hold(ax4, 'on')
% plot(ax4, VVS(1:6), squeeze(EEtotfff(1:6)*NF)/abs(squeeze(EEtotfff(end))),'LineStyle','none','LineWidth',1.5,'Marker','^','MarkerSize',8,'MarkerEdgeColor','#0057E7','MarkerFaceColor','none')
% xticks(ax4, [-2 -1.5 -1]);
% xlim([-2 -1])
% yticks(ax4, [-2.5 0]*10^(-3))
% xticklabels(ax4, [-2 -1.5 -1])
% yticklabels(ax4, [-0.0025 0])

%% condensation energy -Vs,Vp plot %%%%
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
%% load previous data
% load('aRuO2AMSCdata.mat')
%% Spontaneously broken TRS by phase difference between p11 and p22 %%%%
% phi_width_num=100;
% pphi=linspace(0,2*pi,phi_width_num);
% beta_width_num=10;
% bbeta=linspace(0.01,0.1,beta_width_num);
% phiEEtotpp11f=zeros(phi_width_num,1);
% phiEEtotpm11f=zeros(phi_width_num,1);
% phiEEtotpx11f=zeros(phi_width_num,1);
% phiEEtotpy11f=zeros(phi_width_num,1);
% b=1;strain='SymmB3cross';
% %%%% p+ strained system %%%%
% for z=1:param_width_num
%     parvecpp=parboxpp{z};
%     Vp=parvecpp(1);
%     Dpp11=parvecpp(2);
%     for BETA=1:length(bbeta)
%         b=bbeta(BETA);
%         [phiEEtotpp11f,correlppst,en1ppst,en2ppst,en3ppst,en4ppst]=SymmB_gap_eq_solver('pp11',strain,t,mu,J,beta,Vp,Dpp11,b,kkx,kky,pphi);
%         [Eminpp11,phiminpp11index]=min(phiEEtotpp11f);
%         pphiminpp11(BETA)=pphi(phiminpp11index);
%     end
%     phiEEEtotpp11{z}=phiEEtotpp11f;
% end
% pppppp=1
% %%%% p- strained system
% for z=1:param_width_num
%     parvecpm=parboxpm{z};
%     Vp=parvecpm(1);
%     Dpm11=parvecpm(2);
%     % for BETA=1:length(bbeta)
%     %     b=bbeta(BETA);
%     [phiEEtotpm11f,correlpmst,en1pmst,en2pmst,en3pmst,en4pmst]=SymmB_gap_eq_solver('pm11',strain,t,mu,J,beta,Vp,Dpm11,b,kkx,kky,pphi);
%     [Eminpm11,phiminpm11index]=min(phiEEtotpm11f);
%     % pphiminpm11(BETA)=pphi(phiminpm11index);
%     % end
%     phiEEEtotpm11{z}=phiEEtotpm11f;
% end
% pmmmmm=1
% %%%% px strained system
% for z=1:param_width_num
%     parvecpx=parboxpx{z};
%     Vp=parvecpx(1);
%     Dpx11=parvecpx(2);
%     %     % for BETA=1:length(bbeta)
%     %     %     b=bbeta(BETA);
%     [phiEEtotpx11f,correlpxst,en1pxst,en2pxst,en3pxst,en4pxst]=SymmB_gap_eq_solver('px11',strain,t,mu,J,beta,Vp,Dpx11,b,kkx,kky,pphi);
%     [Eminpx11,phiminpx11index]=min(phiEEtotpx11f);
%     % pphiminpx11(BETA)=pphi(phiminpx11index);
%     % end
%     phiEEEtotpx11{z}=phiEEtotpx11f;
% end
% pxxxxx=1
% %%%% py strained system
% for z=1:param_width_num
%     parvecpy=parboxpy{z};
%     Vp=parvecpy(1);
%     Dpy11=parvecpy(2);
%     %     % for BETA=1:length(bbeta)
%     %     %     b=bbeta(BETA);
%     [phiEEtotpy11f,correlpyst,en1pyst,en2pyst,en3pyst,en4pyst]=SymmB_gap_eq_solver('py11',strain,t,mu,J,beta,Vp,Dpy11,b,kkx,kky,pphi);
%     [Eminpy11,phiminpy11index]=min(phiEEtotpy11f);
%     %     pphiminpy11(BETA)=pphi(phiminpy11index);
%     % end
%     phiEEEtotpy11{z}=phiEEtotpy11f;
% end
% pyyyyy=1
% %%%% mixed strained system
% for z=1:param_width_num
%     parvecmix=parboxmix{z};
%     Vp=parvecmix(1);
%     Dmix11=parvecmix(2);
%     %     % for BETA=1:length(bbeta)
%     %     %     b=bbeta(BETA);
%     [phiEEtotmix11f,correlmixst,en1mixst,en2mixst,en3mixst,en4mixst]=SymmB_gap_eq_solver('mix11',strain,t,mu,J,beta,Vp,Dmix11,b,kkx,kky,pphi);
%     [Eminmix11,phiminmix11index]=min(phiEEtotmix11f);
%     % pphiminmix11(BETA)=pphi(phiminmix11index);
%     % end
%     phiEEEtotmix11{z}=phiEEtotmix11f;
% end
% pmixmix=1
%% Total occupied energy - phi plot
% figure
% plot(squeeze(pphi),squeeze(phiEEtotpp11f*NF))
% hold on
% plot(squeeze(pphi),squeeze(phiEEtotpm11f*NF), 'o')
% xlabel('\phi');ylabel('Occupied Energy');xticks([0 pi/2 pi 3*pi/2 2*pi]);xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});xlim([0 2*pi]);yticklabels({});
% legend('C_4𝒯 broken p^{+}-wave occupied energy','C_4𝒯 broken p^{-}-wave occupied energy')
% 
% figure
% plot(squeeze(pphi),squeeze(phiEEtotpx11f*NF))
% hold on
% plot(squeeze(pphi),squeeze(phiEEtotpy11f*NF), 'o')
% xlabel('\phi');ylabel('Occupied Energy');xticks([0 pi/2 pi 3*pi/2 2*pi]);xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});yticklabels({});xlim([0 2*pi])
% legend('C_4𝒯 broken p^{x}-wave occupied energy','C_4𝒯 broken p^{y}-wave occupied energy')
% 
% figure
% plot(squeeze(pphi),squeeze(phiEEtotpx11f*NF))
% xlabel('\phi');ylabel('Occupied Energy');xticks([0 pi/2 pi 3*pi/2 2*pi]);xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});yticklabels({});xlim([0 2*pi])
% legend('C_4𝒯 broken mixed p-wave occupied energy')
%% beta vs phimin phase diagram %%%%
% figure
% plot(squeeze(bbeta),squeeze(pphiminpp11))
% % figure
% % plot(squeeze(bbeta),squeeze(pphiminpm11))
% % figure
% % plot(squeeze(bbeta),squeeze(pphiminpx11))
% % figure
% % plot(squeeze(bbeta),squeeze(pphiminpy11))
% % figure
% plot(squeeze(bbeta),squeeze(pphiminmix11))