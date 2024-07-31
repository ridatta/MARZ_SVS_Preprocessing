clc; close all; clear;

fname = checkDir('\Users\rdatta\Dropbox (MIT)\PUFFIN\Data\MARZ\SVS/z3697 svs5/z3697svs5_shot.hdf');
saveDir = checkDir('/Users/Rishabh/Dropbox (MIT)/PUFFIN/Data/MARZ/SVS/');
data = hdfread(fname,'fore');

figure
imagesc(data);
set(gca,'Ydir','Normal');
colormap(hot)
title('SVS5');

ylabel('wavelength (a.u.)');
xlabel('time (a.u.)');
formatPlots(600,16);
set(gcf,'Position', [262   118   820   556]);

exportgraphics(gcf,[saveDir, 'SVS5_raw.tiff']);

%%
clc; close all; 
% Calibration

f1_id = 'z3697svs5_HgNe';
f1 = checkDir(['\Users\rdatta\Dropbox (MIT)\PUFFIN\Data\MARZ\SVS/z3697 svs5/' f1_id '.hdf']);
f1 = hdfread(f1,'fore');
f2_id = 'z3697svs5_543nm';
f2 = checkDir(['\Users\rdatta\Dropbox (MIT)\PUFFIN\Data\MARZ\SVS/z3697 svs5/' f2_id '.hdf']);
f2 = hdfread(f2,'fore');
f3_id = 'z3697svs5_458nm';
f3 = checkDir(['\Users\rdatta\Dropbox (MIT)\PUFFIN\Data\MARZ\SVS/z3697 svs5/' f3_id '.hdf']);
f3 = hdfread(f3,'fore');

fig = figure;
h1 = subplot(1,3,1);
imagesc(f1); set(gca,'Ydir','Normal');
colormap(hot); title(strrep(f1_id,'_','\_'));
h2 = subplot(1,3,2);
imagesc(f2); set(gca,'Ydir','Normal');
colormap(hot); title(strrep(f2_id,'_','\_')); yticklabels({});
h3 = subplot(1,3,3);
imagesc(f3); set(gca,'Ydir','Normal');
colormap(hot); title(strrep(f3_id,'_','\_')); yticklabels({});
colorbar();

set(gcf,'Position',[0 0 1200 300])

% wavelength cal 
temp = double(f3); miny = 400;
temp = temp(miny:end,1:500);
[idr,idc] = find(temp > 4e3);
y1 = round(mean(idr)) + miny;

% figure
% imagesc(temp); set(gca,'Ydir','Normal'); colormap(hot);

temp = double(f2);
temp = temp(miny:end,1:500);
[idr,idc] = find(temp > 4e3);
y2 = round(mean(idr)) + miny;

axes(h1);  set(gca,'TickDir','out');
axes(h3); text(250,y1,f3_id(end-4:end),'Color','w','Fontsize',14); set(gca,'TickDir','out'); 
axes(h2); text(250,y2,f2_id(end-4:end),'Color','w','Fontsize',14); set(gca,'TickDir','out');

exportgraphics(gcf,[saveDir, 'SVS5_wl_calib.tiff']);

% Time calibration

temp = double(f2);

temp = temp(140,:);
[m,n] = size(temp);
[pks,idx] = findpeaks(temp,'MinPeakDistance',50,'MinPeakHeight',1e3);
figure
subplot(2,1,1);
plot(1:n,temp); hold on;
plot(idx,pks,'og'); 
xlabel('x (px)'); xlim([0,n]);
ylabel('counts (a.u.)');

subplot(2,1,2);
t = linspace(2800,3200,numel(pks));
plot(t,idx,'o'); hold on;
P = polyfit(t,idx,1);
tq = 2800:3200;
plot(tq,polyval(P,tq),'-','color','r'); 
xlabel('t (ns)');
ylabel('x (px)');
exportgraphics(gcf,[saveDir, 'SVS5_time_calib.tiff']);

% time cal functions
t2px = @(t) polyval(P,t); % return ns 2 px
px2t = @(x) (x - P(2)) / P(1); % return ns


% get wl calibration linear lam[nm] = m[nm/px]*y[px] + c[nm]
[m,n] = size(f3); 
ypx = 0:m; xpx = 0:n;
lam2 = 543; lam1 = 458; % nm
m = (lam2 - lam1) / (y2 - y1); % nm/px
c = lam1 - m* y1; % nm
px2wl = @(ypx) m * ypx + c;
wl2px = @(lam) (lam - c) / m;


saveDir = checkDir('/Users/Rishabh/Dropbox (MIT)/PUFFIN/Data/MARZ/SVS/');
save([saveDir, 'calibration.mat'],'t2px', 'px2t', 'px2wl', 'wl2px');

% get lineouts

[m,n] = size(data);
t = 1:m; wl = 1:n; % px
t = px2t(t); wl = px2wl(wl); % convert
[tt,wl] = meshgrid(t,wl); 

Lq = linspace(350,700,round(wl2px(700)-wl2px(350))); 
tq = 2950 .* ones(size(Lq)); 
tstp = (t2px(3000)-t2px(2999)) / 2; 
N = 2;

out1 =  getavgLineout(tt,wl,double(f1),tq,Lq,tstp,N);


figure
plot(Lq,out1,'Linewidth',2,'DisplayName',strrep(f1_id,'_','\_')); hold on;

out2 =  getavgLineout(tt,wl,double(f2),tq,Lq,tstp,N);

plot(Lq,out2,'Linewidth',2,'DisplayName',strrep(f2_id,'_','\_')); hold on;

out3 =  getavgLineout(tt,wl,double(f3),tq,Lq,tstp,N);



plot(Lq,out3,'Linewidth',2,'DisplayName',strrep(f3_id,'_','\_')); hold on;


grid on; xlabel('Wavelength (nm)'); ylabel('Counts (a.u.)');
title(['Instrument Response t = ' num2str(tq(1)) ' \pm ' num2str(round(tstp*N)) ' ns']);
set(gcf,'Position',[0 0 1200 600]); legend('location','northwest');
formatPlots()
set(gca, 'YScale', 'log')


% magnifyPlot([570,700],[0 5000],[0.59 0.45 0.3 0.35],2,1);
%  grid on;    
saveDir = checkDir('/Users/Rishabh/Dropbox (MIT)/PUFFIN/Data/MARZ/SVS/');


 saveas(gcf,[saveDir, 'svs5-instrumentbroadening' num2str(tq(1)) '.png']);
 
 % FWHM of the 543 qnm and 458 nm peaks
 figure
 
 subplot(2,2,1)

 [pks,locs] = findpeaks(out3,'MinPeakProminence',50);
 
 idx = Lq >= (Lq(locs) - 15) & Lq <= (Lq(locs) + 15);
 xi = Lq(idx); yi = out3(idx); yi = yi - mean(yi(1:5));
 plot(xi,yi,'-o','Color',[0.9290 0.6940 0.1250]); hold on;
plot(Lq(locs),pks,'X','MarkerSize',10); 

wid = fwhm(xi,yi);

plot([-wid/2 wid/2] + Lq(locs),[pks, pks]*0.5,'r-'); 
text( Lq(locs)+wid/2,pks*0.5,['FWHM = ' num2str(wid,3) ' nm'],'color','r'); 
 
 grid on; xlabel('Wavelength (nm)'); ylabel('Counts (a.u.)');

ylim([0,2.2e4]);
  subplot(2,2,2)
  

 [pks,locs] = findpeaks(out2,'MinPeakProminence',50);

 
 idx = Lq >= (Lq(locs) - 10) & Lq <= (Lq(locs) + 10);
 xi = Lq(idx); yi = out2(idx);yi = yi - mean(yi(1:5));
 plot(xi,yi,'-o','Color',[0.8500 0.3250 0.0980]); hold on;
plot(Lq(locs),pks,'X','MarkerSize',10); 

wid = fwhm(xi,yi);

plot([-wid/2 wid/2] + Lq(locs),[pks, pks]*0.5,'r-'); 
text( Lq(locs)+wid/2,pks*0.5,['FWHM = ' num2str(wid,3) ' nm'],'color','r'); 
 
 grid on; xlabel('Wavelength (nm)'); ylabel('Counts (a.u.)');
ylim([0,2.2e4]);
  xlim([535,555])
%  saveas(gcf,[saveDir, 'svs5-instrumentbroadening_FWHM' num2str(tq(1)) '.png']);
 
  % FWHM of the peaks form HgNe lamp

 subplot(2,2,3);
 idx = Lq > 420 & Lq < 450;
 xi = Lq(idx); yi = out1(idx);yi = yi - mean(yi(1:5));
 
  [pks,locs] = findpeaks(yi,'MinPeakProminence',50);

 plot(xi,yi,'-o'); hold on;
plot(xi(locs),pks,'X','MarkerSize',10); 

wid = fwhm(xi,yi);

plot([-wid/2 wid/2] + xi(locs),[pks, pks]*0.5,'r-'); 
text( xi(locs)+wid/2,pks*0.5,['FWHM = ' num2str(wid,3) ' nm'],'color','r'); 
 
 grid on; xlabel('Wavelength (nm)'); ylabel('Counts (a.u.)');
  ylim([-100,4500]); xlim([420,450]);
 subplot(2,2,4);
 idx = Lq > 535 & Lq < 555;
 xi = Lq(idx); yi = out1(idx);yi = yi - mean(yi(1:5));
 
  [pks,locs] = findpeaks(yi,'MinPeakProminence',50);

 plot(xi,yi,'-o'); hold on;
plot(xi(locs),pks,'X','MarkerSize',10); 

wid = fwhm(xi,yi);

plot([-wid/2 wid/2] + xi(locs),[pks, pks]*0.5,'r-'); 
text( xi(locs)+wid/2,pks*0.5,['FWHM = ' num2str(wid,3) ' nm'],'color','r'); 
 
 grid on; xlabel('Wavelength (nm)'); ylabel('Counts (a.u.)');
 ylim([-100,4500]);

set(gcf,'position',[0 0 900 600]);
  saveas(gcf,[saveDir, 'svs5-instrumentbroadening_FWHM-2' num2str(tq(1)) '.png']);
  
  
  figure
   idx = Lq > 385 & Lq < 415;
 xi = Lq(idx); yi = out1(idx); yi = yi - mean(yi(1:5));
 
  [pks,locs] = findpeaks(yi,'MinPeakProminence',50);
  [pks,id] = max(pks); locs = locs(id);

 plot(xi,yi,'-o'); hold on;
plot(xi(locs),pks,'X','MarkerSize',10); 

wid = fwhm(xi,yi);

plot([-wid/2 wid/2] + xi(locs),[pks, pks]*0.5,'r-'); 
text( xi(locs)+wid/2,pks*0.5,['FWHM = ' num2str(wid,3) ' nm'],'color','r'); 
 
 grid on; xlabel('Wavelength (nm)'); ylabel('Counts (a.u.)');
 

%%
clc; close all;
% Instrument response
saveDir = checkDir('/Users/Rishabh/Dropbox (MIT)/PUFFIN/Data/MARZ/SVS/');
load([saveDir, 'calibration.mat'],'t2px', 'px2t', 'px2wl', 'wl2px');

close all;
f_id = {'z3697svs5_W_ND1_1','z3697svs5_W_ND1_2','z3697svs5_W_ND1_3'};
for ii = 1:length(f_id)
f{ii} = checkDir(['\Users\rdatta\Dropbox (MIT)\PUFFIN\Data\MARZ\SVS/z3697 svs5/' f_id{ii} '.hdf']);
f{ii} = hdfread(f{ii},'fore');
end

Tid = 2950;
% Plot
figure
for ii = 1:length(f_id)
subplot(1,length(f_id),ii);
imagesc(f{ii}); set(gca,'Ydir','Normal');
colormap(hot); title(strrep(f_id{ii},'_','\_'));
xlabels = 2800:50:3200;
xticks(t2px(xlabels'));
xticklabels(xlabels');
xlabel('t (ns)');
ylabels = 350:50:800;
yticks(wl2px(ylabels'));
yticklabels(ylabels');
ylabel('wavelength (nm)');
set(gca,'TickDir','out');
hold on;
plot(t2px([Tid,Tid]),wl2px([350,700]),'--','Color',sqclr('b',3),'Linewidth',2);
end
set(gcf,'Position',[0 100 1600 400]);
exportgraphics(gcf,[saveDir, 'SVS5_background.tiff']);


% Lineouts
[m,n] = size(data);
t = 1:m; wl = 1:n; % px
t = px2t(t); wl = px2wl(wl); % convert
[tt,wl] = meshgrid(t,wl); 


Lq = linspace(350,800,round(wl2px(800)-wl2px(350))); 
tq = Tid .* ones(size(Lq)); 
figure
tot = zeros(size(Lq));
for ii = 1:length(f_id)
tstp = (t2px(3000)-t2px(2999)) / 2; 
N = 2;
out =  getavgLineout(tt,wl,double(f{ii}),tq,Lq,tstp,N);
plot(Lq,out,'Linewidth',2,'DisplayName',strrep(f_id{ii},'_','\_')); hold on;
tot = tot + out;
end

plot(Lq,tot/length(f_id),'Linewidth',2,'DisplayName','Average','Color',[1 0 1 0.25]); hold on;

title(['t = ' num2str(Tid) ' \pm ' num2str(round(tstp*N)) ' ns']); 
xlabel('Wavelength (nm)'); ylabel('Counts (a.u.)'); 
grid on; legend('location','northwest');
formatPlots(700);
set(gcf,'Position',[0 0 1200 600]);
exportgraphics(gcf,[saveDir, 'SVS5_background_lineout.tiff']);


%% lineouts
clc; close all;
% load shot data
fid = 'z3697svs5_shot';
fname = checkDir(['\Users\rdatta\Dropbox (MIT)\PUFFIN\Data\MARZ\SVS/z3697 svs5/' fid '.hdf']);
data = hdfread(fname,'fore');
data = double(data);

f1 = figure();
imagesc(data); set(gca,'Ydir','Normal');
colormap(hot); title(strrep(fid,'_','\_'));
xlabels = 2800:50:3200;
xticks(t2px(xlabels'));
xticklabels(xlabels');
xlabel('t (ns)');
ylabels = 350:50:800;
yticks(wl2px(ylabels'));
yticklabels(ylabels');
ylabel('wavelength (nm)');
set(gca,'TickDir','out');


% lineouts

tid = [2950];
figure(f1);
hold on;
plot(t2px([Tid,Tid]),wl2px([350,700]),'--','Color',sqclr('b',3),'Linewidth',2);

f2 = figure();
Lq = linspace(350,700,round(wl2px(700)-wl2px(350))); 
for ii = 1:numel(tid)
    tq = tid(ii) .* ones(size(Lq)); 
    tstp = (t2px(3000)-t2px(2999)) / 2; N = 2;
    out =  getavgLineout(tt,wl,data,tq,Lq,tstp,N);
    figure(f2);
%     plot(Lq,out/5e3,'Linewidth',2,'DisplayName',[num2str(tid(ii)) ' ns'],'Color',sqclr('b',ii*2-1)); hold on; 
    plot(Lq,out/5e3,'Linewidth',2,'DisplayName',[num2str(tid(ii)) ' ns'],'Color','k'); hold on;
    figure(f1);
    plot(tq,Lq,'--',...
        'HandleVisibility','off','Color',sqclr('b',ii*2-1));
end

figure(f1);
set(gca,'TickDir','out');
exportgraphics(gcf,[saveDir, 'SVS5_calibrated.tiff']);

figure(f2);
title(['t = ' num2str(Tid) ' \pm ' num2str(round(tstp*N)) ' ns']); 
xlabel('Wavelength (nm)'); ylabel('Counts (a.u.)'); 
grid on; legend('location','northwest');
formatPlots(700);
set(gcf,'Position',[0 0 1200 600]);
exportgraphics(gcf,[saveDir, 'SVS5_lines.tiff']);

saveDir = checkDir('/Users/Rishabh/Dropbox (MIT)/PUFFIN/Data/MARZ/SVS/');
save([saveDir, 'svs5-' num2str(tid) '.mat'],'Lq','out');

%% 
clc; close all;

% identifying lines

tid = 2950;
saveDir = checkDir('/Users/Rishabh/Dropbox (MIT)/PUFFIN/Data/MARZ/SVS/');
load([saveDir, 'svs5-' num2str(tid) '.mat']);


figure
subplot(3,1,1)
plot(Lq,out,'k','Linewidth',2); 
 ylabel('Counts (a.u.)'); 
xlim([350,700]);grid on; % set(gca,'Position',[0 0 1200  400])
title(['t = ' num2str(Tid) ' \pm ' num2str(round(tstp*N)) ' ns']); 
set(gca,'xticklabels',[]);

AlII = importfile3('Al-II.csv');
wl = AlII(:,1); rel_int = AlII(:,2);
idx = strcmpi(wl,""); wl(idx) = []; rel_int(idx) = [];
rel_int = erase(rel_int,'hfs'); 
wl = double(wl); rel_int = double(rel_int); 


subplot(3,1,2)
bar(wl,rel_int,5,'FaceColor','r','EdgeColor','r',...
    'DisplayName','Al-II','linewidth',2); hold on;
xlim([350,700]);grid on;
ylabel('Rel. Intensity (a.u.)'); 

AlIII = importfile4('Al-III.csv');
wl = string(AlIII(:,1)); rel_int = string(AlIII(:,2));
idx = strcmpi(wl,"              "); wl(idx) = []; rel_int(idx) = [];
rel_int = erase(rel_int,'d*'); rel_int = erase(rel_int,'d'); 
wl = double(wl); rel_int = double(rel_int); 
idx = isnan(rel_int); rel_int(idx) = []; wl(idx) = [];
[wl,id] = unique(wl);rel_int= rel_int(id);
set(gca,'xticklabels',[]);
legend();
set(gca, 'YScale', 'log')

subplot(3,1,3)
bar(wl,rel_int,5,'FaceColor','b','EdgeColor','b',...
    'DisplayName','Al-III','linewidth',2);
xlim([350,700]);grid on;
ylabel('Rel. Intensity (a.u.)'); 
xlabel('Wavelength (nm)');

set(gcf,'Position',[0 0 1100 900]);
legend();
set(gca, 'YScale', 'log')
saveas(gcf,[saveDir, 'svs5-' num2str(tid) '.png']);



%%


load([saveDir, 'svs5-' num2str(tid) '.mat']);
figure
plot(Lq,out,'k','Linewidth',2); 
 ylabel('Counts (a.u.)'); 
xlim([350,700]);grid on; % set(gca,'Position',[0 0 1200  400])
title(['t = ' num2str(Tid) ' \pm ' num2str(round(tstp*N)) ' ns']); 


