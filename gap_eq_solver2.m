function [Etotf, delta,correl,en1,en2,en3,en4]=gap_eq_solver2(type,NF,t,mu,J,beta,V,kkx,kky,Del)
D=Del;
for s=1:10000
    % D=DDinput(s);
    for i=1:length(kkx)
        for j=1:length(kky)
            kx=kkx(i);
            ky=kky(j);
            kxm=-kkx(i);
            kym=-kky(j);
            % he=(-t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3);
            % hh=(-t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3);
            he=(-t*(cos(kx)+cos(ky))-mu)*pauli(0)+J*sin(kx)*sin(ky)*pauli(3);
            hh=(-t*(cos(kxm)+cos(kym))-mu)*pauli(0)+J*sin(kxm)*sin(kym)*pauli(3);

            Hbdg=blkdiag(he,-transpose(hh));
            if strcmp(type,'s')
                gs=1;

                Hds=(-1)*[0,0,0,D;0,0,-D,0;0,-conj(D),0,0;conj(D),0,0,0];
                Hbdg=Hbdg+(Hds);

                [u,d]=eig(Hbdg);
                F=diag(1./(exp(diag(d)*beta)+1));
                correl(:,:,i,j)=(u*F*u');

                delta(i,j)=-NF*V*conj(gs)*correl(1,4,i,j);
                en1(i,j)=d(1,1);en2(i,j)=d(2,2);en3(i,j)=d(3,3);en4(i,j)=d(4,4);
                % cond(i,j)=NF*V*correl(1,4,i,j)*correl(4,1,i,j);
            else
                ges=kx^2+ky^2;
                gpp=sin(kx)+1i*sin(ky);
                gpm=sin(kx)-1i*sin(ky);
                gpx=sqrt(2)*sin(kx);
                gpy=sqrt(2)*sin(ky);
                gd=kx^2-ky^2;

                % ges=kx^2+ky^2;
                % gpp=kx+1i*ky;
                % gpm=kx-1i*ky;
                % gpx=sqrt(2)*kx;
                % gpy=sqrt(2)*ky;
                % gd=kx^2-ky^2;

                if strcmp(type,'pp11')
                    D11=D*gpp+0*gpm;
                    D22=D*gpp+0*gpm;
                    Hdp=(-1)*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];
                elseif strcmp(type,'pm11')
                    D11=0*gpp+D*gpm;
                    D22=0*gpp+D*gpm;
                    Hdp=(-1)*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];
                elseif strcmp(type,'px11')
                    D11=D*gpx;
                    D22=D*gpx;
                    Hdp=(-1)*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];
                elseif strcmp(type,'py11')
                    D11=D*gpy;
                    D22=D*gpy;
                    Hdp=(-1)*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];
                elseif strcmp(type,'mix11')
                    D11=D*gpp+0*gpm;
                    D22=0*gpp+D*gpm;
                    Hdp=(-1)*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];
                end
                Hbdg=Hbdg+(Hdp);

                [u,d]=eig(Hbdg);
                F=diag(1./(exp(diag(d)*beta)+1));
                correl(:,:,i,j)=(u*F*u');

                if strcmp(type,'pp11')
                    delta(i,j)=-NF*(V)*conj(gpp)*correl(1,3,i,j);
                    en1(i,j)=d(1,1);en2(i,j)=d(2,2);en3(i,j)=d(3,3);en4(i,j)=d(4,4);
                    cond(i,j)=NF*(V)*(conj(gpp)*gpp+conj(gpm)*gpm)*correl(1,3,i,j)*correl(3,1,i,j);
                elseif strcmp(type,'pm11') 
                    delta(i,j)=-NF*(V)*conj(gpm)*correl(1,3,i,j);
                    en1(i,j)=d(1,1);en2(i,j)=d(2,2);en3(i,j)=d(3,3);en4(i,j)=d(4,4);
                    cond(i,j)=NF*(V)*(conj(gpp)*gpp+conj(gpm)*gpm)*correl(1,3,i,j)*correl(3,1,i,j);
                elseif strcmp(type,'px11')
                    delta(i,j)=-NF*(V)*conj(gpx)*correl(1,3,i,j);
                    en1(i,j)=d(1,1);en2(i,j)=d(2,2);en3(i,j)=d(3,3);en4(i,j)=d(4,4);
                    cond(i,j)=NF*(V)*(conj(gpx)*gpx+conj(gpy)*gpy)*correl(1,3,i,j)*correl(3,1,i,j);
                elseif strcmp(type,'py11')
                    delta(i,j)=-NF*(V)*conj(gpy)*correl(1,3,i,j);
                    en1(i,j)=d(1,1);en2(i,j)=d(2,2);en3(i,j)=d(3,3);en4(i,j)=d(4,4);
                    cond(i,j)=NF*V*(conj(gpx)*gpx+conj(gpy)*gpy)*correl(1,3,i,j)*correl(3,1,i,j);                    
                elseif strcmp(type,'mix11')
                    delta(i,j)=-NF*(V)*conj(gpp)*correl(1,3,i,j);
                    en1(i,j)=d(1,1);en2(i,j)=d(2,2);en3(i,j)=d(3,3);en4(i,j)=d(4,4);
                    cond(i,j)=NF*(V)*(conj(gpp)*gpp+conj(gpm)*gpm)*correl(1,3,i,j)*correl(3,1,i,j);
                end
            end
        end
    end
    D=sum(sum(delta));
end
Etotf=sum(en1(en1<0))+sum(en2(en2<0))+sum(en3(en3<0))+sum(en4(en4<0));
end