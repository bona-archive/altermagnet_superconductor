function [correl,en1,en2,en3,en4,delta]=Hbdg_solver(type,NF,t,mu,J,beta,V,kkx,kky,D,booleanInput,qx,qy)
for i=1:length(kkx)
    for j=1:length(kky)
        if booleanInput
            % circkx=circshift(kkx,qx);
            % kx=circkx(i);
            % circky=circshift(kky,qx);
            % ky=circky(j);
            % circkxm=circshift(-kkx,-qy);
            % kxm=circkxm(i);
            % circkym=circshift(-kky,-qy);            
            % kym=circkym(j);
            kx=kkx(i)+qx;
            ky=kky(j)+qy;
            kxm=-kkx(i)+qx;
            kym=-kky(j)+qy;
        else
            kx=kkx(i);
            ky=kky(j);
            kxm=-kkx(i);
            kym=-kky(j);
        end
        he=(-t*(cos(kx)+cos(ky))-mu)*pauli(0)+J*sin(kx)*sin(ky)*pauli(3);
        hh=(-t*(cos(kxm)+cos(kym))-mu)*pauli(0)+J*sin(kxm)*sin(kym)*pauli(3);

        Hbdg=blkdiag(he,-transpose(hh));
        if strcmp(type,'s')
            gs=1;

            Hds=[0,0,0,D;0,0,-D,0;0,-conj(D),0,0;conj(D),0,0,0];
            Hbdg=Hbdg+(Hds);

            [u,d]=eig(Hbdg);
            F=diag(1./(exp(diag(d)*beta)+1));
            correl(:,:,i,j)=(u*F*u');

            delta(i,j)=-NF*V*conj(gs)*correl(1,4,i,j);
            en1(i,j)=d(1,1);en2(i,j)=d(2,2);en3(i,j)=d(3,3);en4(i,j)=d(4,4);
            cond(i,j)=NF*V*correl(1,4,i,j)*correl(4,1,i,j);
        else
            gpp=sin(kx)+1i*sin(ky);
            gpm=sin(kx)-1i*sin(ky);
            gpx=sqrt(2)*sin(kx);
            gpy=sqrt(2)*sin(ky);
            ges=cos(kx) + cos(ky);
            gd=cos(kx) - cos(ky);

            % ges=kx^2+ky^2;
            % gpp=sin(kx)+1i*sin(ky);
            % gpm=sin(kx)-1i*sin(ky);
            % gpx=sqrt(2)*sin(kx);
            % gpy=sqrt(2)*sin(ky);
            % gd=kx^2-ky^2;

            if strcmp(type,'pp11')
                D11=D*gpp+0*gpm;
                D22=D*gpp+0*gpm;
                Hdp=[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];
            elseif strcmp(type,'pm11')
                D11=0*gpp+D*gpm;
                D22=0*gpp+D*gpm;
                Hdp=[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];
            elseif strcmp(type,'px11')
                D11=D*gpx;
                D22=D*gpx;
                Hdp=[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];
            elseif strcmp(type,'py11')
                D11=D*gpy;
                D22=D*gpy;
                Hdp=[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];
            elseif strcmp(type,'mix11')
                D11=D*gpp+0*gpm;
                D22=0*gpp+D*gpm;
                Hdp=[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];
            elseif strcmp(type,'es')
                D=D*ges;
                Hdp=[0,0,0,D;0,0,-D,0;0,-conj(D),0,0;conj(D),0,0,0];
            elseif strcmp(type,'d')
                D=D*gd;
                Hdp=[0,0,0,D;0,0,-D,0;0,-conj(D),0,0;conj(D),0,0,0];
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
                cond(i,j)=NF*(V)*(conj(gpx)*gpx+conj(gpy)*gpy)*correl(1,3,i,j)*correl(3,1,i,j);
            elseif strcmp(type,'mix11')
                delta(i,j)=-NF*(V)*conj(gpp)*correl(1,3,i,j);
                en1(i,j)=d(1,1);en2(i,j)=d(2,2);en3(i,j)=d(3,3);en4(i,j)=d(4,4);
                cond(i,j)=NF*(V)*(conj(gpp)*gpp+conj(gpm)*gpm)*correl(1,3,i,j)*correl(3,1,i,j);
            elseif strcmp(type,'es')
                delta(i,j)=-NF*V*conj(ges)*correl(1,4,i,j);
                en1(i,j)=d(1,1);en2(i,j)=d(2,2);en3(i,j)=d(3,3);en4(i,j)=d(4,4);
                cond(i,j)=NF*V*correl(1,4,i,j)*correl(4,1,i,j);
            elseif strcmp(type,'d')
                delta(i,j)=-NF*V*conj(gd)*correl(1,4,i,j);
                en1(i,j)=d(1,1);en2(i,j)=d(2,2);en3(i,j)=d(3,3);en4(i,j)=d(4,4);
                cond(i,j)=NF*V*correl(1,4,i,j)*correl(4,1,i,j);
            end
        end
    end
end
end