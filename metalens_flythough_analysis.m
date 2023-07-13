cuts = 10;             %number of cuts
mag = 1*200/200;       %magnification
pp = 1.00;               %pixel pitch
cut_um = 50;            %length of the line cuts [um]
leng = round(cut_um*mag/pp);%pixles

filename = '110_OFF1_rescale.png';
try
    I = im2double(imread(filename));
    try
        I = rgb2gray(I);
    catch
    end
catch
end



figure(1);clf;
I = (I./max(I,[],'all')); %normalize peak
imshow(I)
[m,n] = size(I);
I(I<0.04) = 0; %noise threshold
intExp = sum(sum(I)); %integral experimental

lengout=leng;

%not the best way to do this if not pseudo-diffraction limited
%should be gaussian fitting but takes a long time for 2d
%if similar to gaussian MPE same as fitting
[~,x] = max(mean(I,1));
[~,y] = max(mean(I,2));

% [x,y] = getpts;
cpt = [x,y];
figure(1);
hold on
plot(cpt(1),cpt(2),'r+','LineWidth',2)
theta = linspace(0,2*pi,cuts+1);
theta = theta(1:end-1);
pts = 7*lengout;%sampling density for interpolation
line = zeros(2,pts);

[X,Y] = meshgrid(1:m,1:n);
R = sqrt((X-cpt(1)).^2+(Y-cpt(2)).^2);
roz= 1.00; %gives the correct FWHM of airy disk
a=roz.*R;
intad=(2*besselj(1,a)./a).^2;% tint = trapz(cutx,intad);
intad(cpt(2),cpt(1)) = 1;
intTheo = sum(sum(intad));
sf = (intTheo/intExp);%normalize to the theoretical intensity, same spacing ect.
I = I.*sf;

line(1,:) = linspace(-lengout,lengout,pts);
line(2,:) = zeros(1,pts);
cut = zeros(cuts,pts);
srs = zeros(1,cuts);
for i = 1:cuts
    %line cuts
    rline = [cos(theta(i)),-sin(theta(i));sin(theta(i)),cos(theta(i))]*line;
    rline(1,:) = rline(1,:) + cpt(1);
    rline(2,:) = rline(2,:) + cpt(2);
    for k = 1:pts
        %interp since camera res is low
        cut(i,k) = interp2(1:n,1:m,I,rline(1,k),rline(2,k));
    end
    %cool plotting, not necessary
    if mod(i,5)==0
        figure(1);
        plot(rline(1,:),rline(2,:))
    end
    cutx = linspace(-lengout,lengout,pts).*pp/mag;
    %SR (if cpt not excatly centered max not necessarialy on axis)
    srs(i) = max(cut(i,:));
    figure(3);clf;
    plot(cutx,cut(i,:))
    grid on
    hold on
    a=1.0*cutx;
    theo=(2*besselj(1,a)./a).^2;
    plot(cutx,theo)  
end

cutm = mean(cut);
cutx = linspace(-lengout,lengout,pts)*pp/mag;

nf = mean(srs);

figure(2);clf;
plot(cutx,nf*(cutm-0.5*max(cutm)))
grid on

figure(3);clf;
plot(cutx,(cutm))
grid on
fprintf('IntE = %1.3e, ',trapz(cutx,nf*(cutm./max(cutm))))
k=2*pi/0.633;r=0.5E4;z=6E4; 
a=1.00.*cutx;
intad=(2*besselj(1,a)./a).^2;
hold on
plot(cutx,intad)  
fprintf('IntT = %1.3e\n',trapz(cutx,intad))
legend('Experimental','Theoretical')
saveas(gcf,strcat(filename(1:end-4),'-Comp.jpg'),'jpg')



figure(4);clf;
histogram(srs)
title(strcat('Distribution of SRs: ',num2str(mean(srs)),'+-',num2str(std(srs))))
saveas(gcf,strcat(filename(1:end-4),'-SR.jpg'),'jpg')

