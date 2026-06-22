clear all
%% parameters
t=1;
J_ran=linspace(0,1,11);
mu=0;
m=0;
U=-1;
% V=-1.5;
tl = sprintf('t=%.3f, U=%.3f, mu=%.3f', t, U, mu);
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

DDslim=4;
DDpp11lim=4;
DDpm11lim=4;
DDpx11lim=4;
DDpy11lim=4;
DDmix11lim=4;
%%%% number of widths of input orderparams %%%%
slice=100;
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
VVS=linspace(U,U,1);
VVP=linspace(-1,-1,1);
%%%% self-consistent occupied energy receivers %%%%
EEtotsf=zeros(length(VVP), length(J_ran)); 
EEtotpp11f=zeros(length(VVP),length(J_ran));
EEtotpm11f=zeros(length(VVP),length(J_ran));
EEtotpx11f=zeros(length(VVP),length(J_ran));
EEtotpy11f=zeros(length(VVP),length(J_ran));
EEtotmix11f=zeros(length(VVP),length(J_ran));
%%%% self-consistent input orderparam indexes receivers %%%%
Indexs=zeros(length(VVP),length(J_ran));
Indexpp11=zeros(length(VVP),length(J_ran));
Indexpm11=zeros(length(VVP),length(J_ran));
Indexpx11=zeros(length(VVP),length(J_ran));
Indexpy11=zeros(length(VVP),length(J_ran));
Indexmix11=zeros(length(VVP),length(J_ran));
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
% EEf=zeros(1,NOGx);
% qqqq=zeros(NOGx,length(DDsinput));
% DDffout=zeros(length(VVP),length(DDsinput));
% % EEtotff=zeros(length(VVP),length(DDsinput));
% EEtotfff=zeros(length(VVP),length(J_ran));
% Indexff=zeros(length(VVP),length(J_ran));
BFSff=zeros(1,width_num);
%% gap_eq_solver
for a=1:length(VVP)
    for j=1:length(J_ran)
    Vs=VVS(1);
    Vp=VVP(a)
    J=J_ran(j)

    legendLabels{j}=sprintf('J=%.2f', J);
    legendLabelp{j}=sprintf('J=%.2f', J);

    [DDsout,EEtots,cconds,ccorrels,en1s,en2s,en3s,en4s,ddeltas]=gap_eq_solver('s',NF,t,mu,J,beta,Vs,kkx,kky,DDsinput); %s-wave gap_eq_solver
    % DDDsout{j}=DDsout;
    % EEEtots{j}=EEtots;    
    [Nones,I]=min(abs(DDsout./(DDsinput+epsilon)-1));
    if (DDsout(2)-DDsout(1))/(DDsinput(2)-DDsinput(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotsf(j)=EEtots(I);
    Delta_s(j)=DDsout(I);      
    Indexs(j)=I;
    [DDpp11out,EEtotpp11,ccondpp,ccorrelpp,en1pp,en2pp,en3pp,en4pp,ddeltapp]=gap_eq_solver('pp11',NF,t,mu,J,beta,Vp,kkx,kky,DDpp11input);  %pp-wave gap_eq_solver
    % DDDpp11out{j}=DDpp11out;
    % EEEtotpp11{j}=EEtotpp11;
    [Nonepp11,I]=min(abs(DDpp11out./(DDpp11input+epsilon)-1));
    if (DDpp11out(2)-DDpp11out(1))/(DDpp11input(2)-DDpp11input(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotpp11f(j)=EEtotpp11(I);
    Indexpp11(j)=I;
    Delta_pp(j)=DDpp11out(I);      
    % parboxpp{j}=[Vp;DDpp11input(I)];

    [DDpx11out,EEtotpx11,ccondpx,ccorrelpx,en1px,en2px,en3px,en4px,ddeltapx]=gap_eq_solver('px11',NF,t,mu,J,beta,Vp,kkx,kky,DDpx11input); %px-wave gap_eq_solver
    % DDDpx11out{j}=DDpx11out;
    % EEEtotpx11{j}=EEtotpx11;    
    [Nonepx11,I]=min(abs(DDpx11out./(DDpx11input+epsilon)-1));
    if (DDpx11out(2)-DDpx11out(1))/(DDpx11input(2)-DDpx11input(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotpx11f(j)=EEtotpx11(I);
    Indexpx11(j)=I;
    Delta_px(j)=DDpx11out(I);      
    % parboxpx{j}=[Vp;DDpx11input(I)];

    [DDmix11out,EEtotmix11,ccondmix,ccorrelmix,en1mix,en2mix,en3mix,en4mix,ddeltamix]=gap_eq_solver('mix11',NF,t,mu,J,beta,Vp,kkx,kky,DDmix11input); %pmix-wave gap_eq_solver
    % DDDmix11out{j}=DDmix11out;
    % EEEtotmix11{j}=EEtotmix11;
    [Nonemix11,I]=min(abs(DDmix11out./(DDmix11input+epsilon)-1));
    if (DDmix11out(2)-DDmix11out(1))/(DDmix11input(2)-DDmix11input(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotmix11f(j)=EEtotmix11(I);
    Indexmix11(j)=I;
    Delta_mix(j)=DDmix11out(I);    
    % parboxmix{j}=[Vp;DDmix11input(I)];

    for k=1:length(DDsinput)
        kkxq=linspace(-pi,pi,NOGx+1);kkxq(end)=[];
        kkyq=linspace(-pi,pi,NOGx+1);kkyq(end)=[];
        qqx=linspace(-pi,pi,NOGx+1);qqx(end)=[];
        qqy=linspace(-pi,pi,NOGy+1);qqy(end)=[];
        for i=1:length(qqx)
            [correlf,en1f,en2f,en3f,en4f,deltaf]=Hbdg_solver('s',NF,t,mu,J,beta,Vs,kkx,kky,DDsinput(k),true,qqx(i),0);
            EEf(i)=sum(en1f(en1f<0))+sum(en2f(en2f<0))+sum(en3f(en3f<0))+sum(en4f(en4f<0));
        end
        Ef=fix((EEf*10^(-3))*10^(13))/10^13;
        [xx, yy]=find(Ef==min(Ef));
        qqq=qqx(yy(1));
        qqqq(j,k)=qqq;
        [correlff,en1ff,en2ff,en3ff,en4ff,deltaff]=Hbdg_solver('s',NF,t,mu,J,beta,Vs,kkxq,kkyq,DDsinput(k),true,qqq,0);
        DDffout(k)=sum(sum(deltaff));
        EEtotff(k)=sum(en1ff(en1ff<0))+sum(en2ff(en2ff<0))+sum(en3ff(en3ff<0))+sum(en4ff(en4ff<0));
        % if sum(en1(en1<0))==0 | sum(en2(en2<0))==0 | sum(en3(en3<0))==0 | sum(en4(en4<0)) ==0
        %     BFSff=true;
        % else
        %     BFSff=false;
        % end
    end
    % DDDffout{j}=DDffout;
    % EEEtotff{j}=EEtotff;    
    [Noneff,I]=min(abs(DDffout./(DDsinput+epsilon)-1));
    if (DDffout(2)-DDffout(1))/(DDsinput(2)-DDsinput(1)) < 1
        I=1;
    else
        I=I;
    end
    EEtotfff(j)=EEtotff(I);
    Indexff(j)=I;
    Delta_ff(j)=DDffout(I);
    end
end

%% 
figure
hold on
yyaxis left
ylabel('Analytically solved')
plot(Delta_s_lst)
yyaxis right
ylabel('Numerically ED')
plot(Delta_s)
%% 
% allMatrices = cat(3, EEtotsf, EEtotpp11f);
allMatrices = cat(3, EEtotsf, EEtotpp11f,EEtotfff);
[~, minIndex] = min(allMatrices, [], 3);
figure
surf(J_ran,VVP,transpose(minIndex))
xlabel('J')
ylabel('V')
title(tl)
view(2)
%% 

% figure
% surf(VVP,J_ran,EEtotpp11f-EEtotfff);xlabel('V');ylabel('J');title('E0p-E0ff');view(2)