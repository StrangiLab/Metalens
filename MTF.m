clearvars;
m = 50;
n = 50;
cpt = [m/2,n/2];
dgrid = 0.25;
[X,Y] = meshgrid(0:dgrid:m,0:dgrid:n);
R = sqrt((X-cpt(1)).^2+(Y-cpt(2)).^2);
roz= 1.00; %gives the correct FWHM of airy disk
a=roz.*R;
intad=(2*besselj(1,a)./a).^2;% tint = trapz(cutx,intad);
intad(isnan(intad)) = 1.2;
intad = intad./sum(intad,'all');
[~,cline] = max(sum(intad));
psfTheo = intad(cline,:);
mtfTheo = abs(fftshift(fft(psfTheo)));
[~,cline] = max(mtfTheo);
mtfTheo = mtfTheo(cline+1:cline+40)/max(mtfTheo);

freq = linspace(0,1,31);
defLimit = (2/pi).*(acos(freq)-abs(freq).*sqrt(1-freq.^2));
clf;
plot(mtfTheo)
hold on
plot(defLimit,'k-')

emptyData = importdata('uninfiltrated_ffp.txt');
empty = zeros(99,1000);
empty(1:50,:) = flip(emptyData);
empty(50:end,:) = (emptyData);
[~,a] = max(empty(50,:));
psfEmp = empty(:,a);
mtfEmp = abs(fftshift(fft(psfEmp)));
mtfEmp = mtfEmp/max(mtfEmp);


isodata = importdata('results-1.58ind_0.25h_FFP.txt');
isotropic = zeros(199,1334);
isotropic(1:100,:) = flip(isodata);
isotropic(100:end,:) = (isodata);
[~,a] = max(isotropic(100,:));
areaIso = sum(isotropic(100-50:100+50,a));
psfIso = isotropic(101-50:101+50,a);
mtfIso = abs(fftshift(fft(psfIso)));
mtfIso = mtfIso(51:51+40)/max(mtfIso);

nematicData = importdata('results1-3750cpt_0.25pinf_FFP.txt');
nematic = zeros(99,1000);
nematic(1:50,:) = flip(nematicData);
nematic(50:end,:) = (nematicData);
[~,a] = max(nematic(50,:));
psfNem = nematic(:,a);
mtfNem = abs(fftshift(fft(psfNem)));
mtfNem = mtfNem/max(mtfNem);






figure(2);clf;
hold on
% plot(-24:25,psfIso(100-24:100+25)/max(psfIso))
% plot(-24:25,psfNem(50-24:50+25)/max(psfNem))
% plot(-24:25,psfEmp(50-24:50+25)/max(psfEmp))
plot(-24:25,psfTheo(50-24:50+25))


figure(3);clf;
hold on
plot(0:32,mtfTheo,'k-','linewidth',2)
plot(0:16,mtfIso)
grid on
