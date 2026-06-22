function Jk=Alterdispersion(kkx,kky,J)
Jk=zeros(length(kkx),length(kky));
for i=1:length(kkx)
    for j=1:length(kky)
        kx=kkx(i);
        ky=kky(j);
        Jk(i,j)=J*sin(kx)*sin(ky);
    end
end
end