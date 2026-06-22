clear all
t=1;
J=1;
mu=-1;
del=0.1;
epsilon=1e-6;
N=500;
kx=linspace(-pi,pi,N);kx(end)=[];
ky=linspace(-pi,pi,N);ky(end)=[];

for l=1:length(kx)
    q=kx(l);
    for i=1:length(kx)
        for j=1:length(ky)
            kkx=kx(i);
            kky=ky(j);
            A=-t*(cos(kkx+q)+cos(kky+q))-mu;
            B=J*sin(kkx+q)*sin(kky+q);
            C=del;
            AA=-t*(cos(-kkx+q)+cos(-kky+q))-mu;
            BB=J*sin(-kkx+q)*sin(-kky+q);
            HFFLO=[A+B 0 0 -C;0 A-B C 0;0 C -AA-BB 0;-C 0 0 -AA+BB];
            HFFLO1=[A+B 0 0 -C;0 A-B C 0;0 C -AA-BB 0;-C 0 0 -AA+BB];
            HFFLO2=blkdiag(HFFLO,HFFLO1);
            % [u,d]=eig(HFFLO);
            % E1(i)=d(1,1);E2(i)=d(2,2);E3(i)=d(3,3);E4(i)=d(4,4);
            [u,d]=eig(HFFLO2);
            EEE1(i,j,l)=d(1,1);EEE2(i,j,l)=d(2,2);EEE3(i,j,l)=d(3,3);EEE4(i,j,l)=d(4,4);EEE5(i,j,l)=d(5,5);EEE6(i,j,l)=d(6,6);EEE7(i,j,l)=d(7,7);EE8(i,j,l)=d(8,8);
            EE1(i,j,l)=1/2*((-AA-BB+A-B)+sqrt((-AA-BB-A+B)^2+4*C^2));
            EE2(i,j,l)=1/2*((-AA-BB+A-B)-sqrt((-AA-BB-A+B)^2+4*C^2));
            EE3(i,j,l)=1/2*((-AA+BB+A+B)+sqrt((-AA+BB-A-B)^2+4*C^2));
            EE4(i,j,l)=1/2*((-AA+BB+A+B)-sqrt((-AA+BB-A-B)^2+4*C^2));
            e1(i,j)=1/2*((-AA-BB+A-B)+sqrt((-AA-BB-A+B)^2+4*C^2));
            e2(i,j)=1/2*((-AA-BB+A-B)-sqrt((-AA-BB-A+B)^2+4*C^2));
            e3(i,j)=1/2*((-AA+BB+A+B)+sqrt((-AA+BB-A-B)^2+4*C^2));
            e4(i,j)=1/2*((-AA+BB+A+B)-sqrt((-AA+BB-A-B)^2+4*C^2));
        end
    end
    localmin1(l)=min(min(EE1(:,:,l)));
    localmin3(l)=min(min(EE3(:,:,l)));
    localmax2(l)=max(max(EE2(:,:,l)));
    localmax4(l)=max(max(EE4(:,:,l)));
    E(l)=sum(e1(e1<0))+sum(e2(e2<0))+sum(e3(e3<0))+sum(e4(e4<0));
end
A1=abs(localmin1);
A3=abs(localmin3);
A2=abs(localmax2);
A4=abs(localmax4);
TF1=islocalmin(A1);I1=find(TF1==1);
TF3=islocalmin(A3);I3=find(TF3==1);
TF2=islocalmin(A2);I2=find(TF2==1);
TF4=islocalmin(A4);I4=find(TF4==1);
figure
plot(kx,A1,kx(TF1),A1(TF1),'r*')
xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2','\pi'});
ylabel('gap \Delta')

figure
plot(kx,A2,kx(TF2),A2(TF2),'r*')
xticks([-pi -pi/2 0 pi/2 pi]);xticklabels({'-\pi', '-\pi/2', '0', '\pi/2','\pi'});
ylabel('gap \Delta')
%% 


ind=3;
figure
plot(EE1(:,I1(ind)))
hold on
plot(EE2(:,I2(ind)))
plot(EE3(:,I3(ind)))
plot(EE4(:,I4(ind)))
plot(localmin1(I1(ind))*ones(1,N),'k')
plot(localmin3(I3(ind))*ones(1,N),'k')
plot(localmax2(I2(ind))*ones(1,N),'k')
plot(localmax4(I4(ind))*ones(1,N),'k')