figure(5);clf;
nf=0.59;
plot(cutx+0.5,nf*(cutm./max(cutm)))
grid on
fprintf('IntE = %1.3e, ',trapz(cutx,nf*(cutm./max(cutm))))
k=2*pi/0.633;
r=0.5E4;z=6E4; 
roz= 1.02;

a=roz*cutx;
intad=(2*besselj(1,a)./a).^2;
intad(cutx<=-3.8255) = 0;
intad(cutx>=3.8255) = 0;

hold on
plot(cutx,intad)  
fprintf('IntT = %1.3e\n',trapz(cutx,intad))
legend('Experimental','Theoretical')
