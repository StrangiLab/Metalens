m = 50;n = 50;dgrid = 0.1;
cpt = [m/2,n/2];
[X,Y] = meshgrid(0:dgrid:m,0:dgrid:n);
R = sqrt((X-cpt(1)).^2+(Y-cpt(2)).^2);
roz= 0.94; a=roz.*R;%gives the correct FWHM of airy disk
intad=(2*besselj(1,a)./a).^2;
intad(isnan(intad)) = 1.0;
intad = intad./sum(intad,'all');
[~,cline] = max(sum(intad)); psfTheo = intad(cline,:);

otfTheo = psf2otf(intad);
mtfTheo = abs(otfTheo);
mtfTheo(1,1) = 1.025;
mtfTheo = mtfTheo/max(mtfTheo,[],'all');

emptyData = importdata('uninfiltrated_ffp.txt');
empty = zeros(99,1000);
empty(1:50,:) = flip(emptyData);
empty(50:end,:) = (emptyData);
[~,a] = max(empty(50,:));
psfEmp = empty(:,a);
otfEmp = psf2otf(psfEmp);

mtfEmp = abs(otfEmp);
mtfEmp = mtfEmp/max(mtfEmp);


isodata = importdata('results-1.58ind_0.25h_FFP.txt');
isotropic = zeros(199,1334);
isotropic(1:100,:) = flip(isodata);
isotropic(100:end,:) = (isodata);
[~,a] = max(isotropic(100,:));
areaIso = sum(isotropic(100-50:100+50,a));
psfIso = isotropic(101-50:101+50,a);
% psfIso = interp1(-50:50,psfIso,-50:1.5:50);

otfIso = psf2otf(psfIso);
mtfIso = abs(otfIso);
mtfIso = mtfIso/max(mtfIso);

nematicData = importdata('results1-3750cpt_0.25pinf_FFP.txt');
nematic = zeros(99,1000);
nematic(1:50,:) = flip(nematicData);
nematic(50:end,:) = (nematicData);
[~,a] = max(nematic(50,:));
psfNem = nematic(:,a);

otfNem = psf2otf(psfNem);
mtfNem = abs(otfNem);
mtfNem = mtfNem/max(mtfNem);



figure(1);clf;
hold on
plot(linspace(0,400,20),mtfTheo(1,1:20),'k-','linewidth',2)
plot((1:10:400),mtfEmp(1:40),'bo','linewidth',2)
plot((1:10:600)*0.66,mtfIso(1:60),'ro','linewidth',2)
plot((1:10:400),mtfNem(1:40),'go','linewidth',2)

% freq = linspace(0,1,32);
% defLimit = (2/pi).*(acos(freq)-abs(freq).*sqrt(1-freq.^2));
% hold on
% plot(linspace(0,320,32),defLimit)

figure(2);clf
plot(mtfIso)


