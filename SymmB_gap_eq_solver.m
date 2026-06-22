function [phiEEtotf,correl,en1,en2,en3,en4]=SymmB_gap_eq_solver(type,SymmBtype,t,mu,J,beta,V,D,b,kkx,kky,pphi)
for s=1:length(pphi)
    phi=pphi(s);
    for i=1:length(kkx)
        for j=1:length(kky)
            kx=kkx(i);
            ky=kky(j);
            kxm=-kkx(i);
            kym=-kky(j);

            if strcmp(SymmBtype,'Zeeman')
                SymmB=m*pauli(1);
                SymmBm=m*pauli(1);
            elseif strcmp(SymmBtype,'SymmB1')
                SymmB=kx*pauli(1)+ky*pauli(2);
                SymmBm=kxm*pauli(1)+kym*pauli(2);
            elseif strcmp(SymmBtype,'SymmB2')
                SymmB=kx^2*pauli(1)+ky^2*pauli(2);
                SymmBm=kxm^2*pauli(1)+kym^2*pauli(2);
                % elseif strcmp(SymmBtype,'SymmB2cross') %TRS broken
                % SymmB=kx^2*pauli(2)+ky^2*pauli(1);
                % SymmBm=kxm^2*pauli(2)+kym^2*pauli(1);
            elseif strcmp(SymmBtype,'SymmB3cross')
                SymmB=kx^3*pauli(2)+ky^3*pauli(1);
                SymmBm=kxm^3*pauli(2)+kym^3*pauli(1);
                % elseif strcmp(SymmBtype,'SymmB')
                % % SymmBreaker3xxyx=kx^2*ky*pauli(1);
                % % SymmBreakerm3xxyx=kxm^2*kym*pauli(1);
                % elseif strcmp(SymmBtype,'SymmB')
                % % SymmBreaker3xxyy=kx^2*ky*pauli(2);
                % % SymmBreakerm3xxyy=kxm^2*kym*pauli(2);
                % elseif strcmp(SymmBtype,'SymmB')
                % % SymmBreaker3xyyx=kx*ky^2*pauli(1);
                % % SymmBreakerm3xyyx=kxm*kym^2*pauli(1);
                % elseif strcmp(SymmBtype,'SymmB')
                % % SymmBreaker3xyyy=kx*ky^2*pauli(2);
                % % SymmBreakerm3xyyy=kxm*kym^2*pauli(2);
                % elseif strcmp(SymmBtype,'SymmB')
                % % SymmBreaker3yyyy=ky^3*pauli(2);
                % % SymmBreakerm3yyyy=kym^3*pauli(2);
            end
            he=(-t*(kx^2+ky^2)-mu)*pauli(0)+J*kx*ky*pauli(3)+b*SymmB;
            hh=(-t*(kxm^2+kym^2)-mu)*pauli(0)+J*kxm*kym*pauli(3)+b*SymmBm;
            Hbdg=blkdiag(he,-transpose(hh));

            ges=kx^2+ky^2;
            gpp=kx+1i*ky;
            gpm=kx-1i*ky;
            gpx=sqrt(2)*kx;
            gpy=sqrt(2)*ky;
            gd=kx^2-ky^2;
            if strcmp(type,'pp11')
                D11=D*gpp+0*gpm;
                D22=(D*gpp+0*gpm)*exp(1i*phi);
                Hdp=V*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];
            elseif strcmp(type,'pm11')
                D11=0*gpp+D*gpm;
                D22=(0*gpp+D*gpm)*exp(1i*phi);
                Hdp=V*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];
            elseif strcmp(type,'px11')
                D11=D*gpx+0*gpy;
                D22=(D*gpx+0*gpy)*exp(1i*phi);
                Hdp=V*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];
            elseif strcmp(type,'py11')
                D11=0*gpx+D*gpy;
                D22=(0*gpx+D*gpy)*exp(1i*phi);
                Hdp=V*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];
            elseif strcmp(type,'mix11')
                D11=D*gpp+0*gpm;
                D22=(0*gpp+D*gpm)*exp(1i*phi);
                Hdp=V*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];
            end

            Hdp=V*[0,0,D11,0;0,0,0,D22;conj(D11),0,0,0;0,conj(D22),0,0];

            Hbdg=Hbdg+(Hdp);

            [u,d]=eig(Hbdg);
            F=diag(1./(exp(diag(d)*beta)+1));
            correl(:,:,i,j)=(u*F*u');
            en1(i,j)=d(1,1);en2(i,j)=d(2,2);en3(i,j)=d(3,3);en4(i,j)=d(4,4);
        end
    end
    phiEEtotf(s)=sum(en1(en1<0))+sum(en2(en2<0))+sum(en3(en3<0))+sum(en4(en4<0));
end
end