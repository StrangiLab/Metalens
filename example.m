clear all

m = 50;
n = 50;
cpt = [m/2,n/2];
dgrid = 0.25;
[X,Y] = meshgrid(0:dgrid:m,0:dgrid:n);
R = sqrt((X-cpt(1)).^2+(Y-cpt(2)).^2);
roz= pi*10000/0.633; %gives the correct FWHM of airy disk
a=roz.*R./50000;
intad=(2*besselj(1,a)./a).^2;
intad(isnan(intad)) = 1.0;
intad = intad./sum(intad,'all');
[~,cline] = max(sum(intad));
psfTheo = intad(cline,:);

x = 0:dgrid:m;
y = psfTheo;

% x=[-25.132,-18.849,-12.566,-6.283,-3,0,3,6.283,12.566,18.849,25.132];
% y=[0.21446,0.21446,0.21446,0.21446,44,44.15754,44,0.21446,0.21446,0.21446,0.21446];
LSF = [x,y]; % This is more like a line spread function
OTF = fftshift(fft(psfTheo)); % OTF
MTF = abs(OTF); % absolute value of OTF
MTF = MTF./max(MTF); % normalize
% correct sampling frequency for conversion on frequency bins to frequency
Size = length(MTF);
spacing = 0.001*dgrid;% spacing between data points 
fsx = abs(1/spacing); % turn into sampling frequency
a = linspace(-Size/2,Size/2,Size); % form scale for conversion based on frequency bins
conversionx = fsx./Size; % conversion factor for frequency bin units to frequency (unit^-1)
Psi = a.*conversionx; % frequency (unit^-1)
plot(Psi,MTF)
xlim([0,400])
xlabel('cycles/unit')
ylabel('MTF')
