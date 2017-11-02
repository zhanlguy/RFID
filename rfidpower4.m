clear all
close all
clc

%% Parameters
SNR = 30;
SNRlinear = 10^(SNR/10);
snr = 1/SNRlinear;
lane = 3;
noise = snr;
c  = 3e8;
f  = 9.6e6;
lamd = c/f;

%% Set the map
stepx = .01;
stepy = .01;
igrid = 0:stepx:9;
jgrid = -5:stepy:5;
gridsize = [length(jgrid) length(igrid)];
xI = ones(gridsize(1),1) * igrid;
yI = (ones(gridsize(2),1) * jgrid)';

%  Reciever co-ordinates
Rx1 = lane + 0j;
Rx2 = 2*lane + 0j;
phian1 =[];
phian2 =[];
%% Attenuation is equal to free-space plus ground-reflection
for b = 1:length(jgrid)
    for a = 1:length(igrid)
        x  = igrid(a);
        y  = jgrid(b);
        p  = x + 1i*y;
        D1(b,a) = norm(p-Rx1);
        D2(b,a) = norm(p-Rx2);
        signal1(b,a) = 1./(eps+D1(b,a).^2);
        signal2(b,a) = 1./(eps+D2(b,a).^2);
        Freeloss1(b,a) = D1(b,a);
        Freeloss2(b,a) = D2(b,a);
        phian1(b,a) = atan((lane-x)/y);
        phian2(b,a) = atan((lane*2-x)/y);
    end
end


beamwidth =4;
[rada, radp] = genRadiationPattern4(beamwidth);
ai2=[];
ai1=[];
for ww = 1:length(jgrid)
    for qq= 1:length(igrid)
        ai1{ww,qq} = find(rada>phian1(ww,qq));
        ai2{ww,qq} = find(rada>phian2(ww,qq));
    end
end


emptyIndex1 = cellfun(@isempty,ai1);       %# Find indices of empty cells
ai1(emptyIndex1) = {1};                    %# Fill empty cells with 1
emptyIndex2 = cellfun(@isempty,ai2);       %# Find indices of empty cells
ai2(emptyIndex2) = {1};                    %# Fill empty cells with 1

PL1=[];
PL2=[];

phase1 = exp((-1i*pi*2.*Freeloss1)./lamd);
phase2 = exp((-1i*pi*2.*Freeloss2)./lamd);


for o = 1:b
    for p= 1:a
        ai11 = cell2mat(ai1(o,p));
        ai22 = cell2mat(ai2(o,p));
        PL1{o,p}= radp(ai11(1));
        PL2{o,p}= radp(ai22(1));
    end
end

PLMAP1 = cell2mat(PL1).*phase1.*(signal1);
PLMAP2 = cell2mat(PL2).*phase2.*(signal2);

PLMAP11 = phase1.*(signal1);
PLMAP22 = phase2.*(signal2);


MAP=[];
for ii = 1:b
    for jj = 1:a
        if D1(ii,jj) < D2(ii,jj)
            DMap(ii,jj) = PLMAP1(ii,jj)./PLMAP2(ii,jj)+noise;
            OMap(ii,jj) = PLMAP11(ii,jj)./PLMAP22(ii,jj)+noise;
        else
            DMap(ii,jj) = PLMAP2(ii,jj)./PLMAP1(ii,jj)+noise;
            OMap(ii,jj) = PLMAP22(ii,jj)./PLMAP11(ii,jj)+noise;
        end
    end
end


figure(1);
pcolor(xI, yI,10*log10(abs(PLMAP1).^2)); shading flat;
colorbar
figure(2);
pcolor(xI, yI,10*log10(abs(PLMAP2).^2)); shading flat;
colorbar
figure(3);
pcolor(xI, yI,10*log10(abs(DMap).^2)); shading flat;
colorbar
figure(4);
pcolor(xI, yI,10*log10(abs(PLMAP11).^2)); shading flat;
colorbar
figure(5);
pcolor(xI, yI,10*log10(abs(PLMAP22).^2)); shading flat;
colorbar
figure(6);
pcolor(xI, yI,10*log10(abs(OMAP).^2)); shading flat;
colorbar