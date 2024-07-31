clc; close all; clear;

fname = checkDir('\Users\rdatta\Dropbox (MIT)\PUFFIN\Codes\RadTrans Analysis\Data\exp_spectra.csv');
data = dlmread(fname,',');

wl = data(:,1);
idx = wl > 400; 
wl = wl(idx);

Iout = data(:,4);
Iout = Iout(idx);
% Iout = (Iout - min(Iout)) ./ (max(Iout) - min(Iout)); 


f1 = figure();
plot(wl, Iout);
xlabel('nm');
ylabel('Intensity');


% attenutaion data


fname = checkDir('\Users\rdatta\Dropbox (MIT)\PUFFIN\Codes\RadTrans Analysis\Data\attenuation.csv');
attn = dlmread(fname,',');
dist = 85.95 * 1e-3; % km
figure
plot(attn(:,1), attn(:,2) * dist,marker='o');

xlabel('nm');
ylabel('dB');

dB = interp1(attn(:,1),attn(:,2)*dist,wl);

Iin = Iout .* 10.^(dB/10); 
Iin = Iin ./ max(Iin);
% Iin = (Iin - min(Iin)) ./ (max(Iin) - min(Iin));

figure(f1);
hold on;
plot(wl, Iin,linewidth=3);


% fit a plankian to the spectrum

figure
plot(wl,Iin,linewidth=3,DisplayName='Data',color='r'); hold on;

addpath('/Users/rishabhdatta/Dropbox (MIT)/PUFFIN/Codes/Other/PlasmaFormulary/')
load physicalConstants-SI.mat T_eV
F = FundamentalPlasma();
T = [1, 5, 10, 50];
for ii =1:numel(T)
bb = F.getPlankianDistribution(T(ii) * T_eV,wl*1e-9);
bb = bb ./ max(bb);
plot(wl,bb,linewidth=3,color=sqclr('b',3+ii),...
    DisplayName=['blackbody fit T [eV] = ', num2str(T(ii))]); hold on;
end
ylabel('Intensity');
xlabel('nm');
formatPlots();

