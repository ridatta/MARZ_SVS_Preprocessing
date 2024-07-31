clc; close all; clear;

fname = checkDir('\Users\rdatta\Dropbox (MIT)\PUFFIN group\Data\Z\MARZ\Shots\z3697\SVS\z3697 svs4/z3697svs4_shot.hdf');
saveDir = checkDir('/Users/Rishabh/Dropbox (MIT)/PUFFIN/Data/MARZ/SVS/');
data = hdfread(fname,'fore');

figure
imagesc(data);
set(gca,'Ydir','Normal');
colormap(hot)
title('SVS4');
ylabel('wavelength (a.u.)');
xlabel('time (a.u.)');
formatPlots(600,16);
set(gcf,'Position', [262   118   820   556]);

exportgraphics(gcf,[saveDir, 'SVS4_raw.tiff']);

%%
clc; close all; 
% Calibration

f1_id = 'z3697svs4_HgNe';
f1 = ['/Users/Rishabh/Dropbox (MIT)/PUFFIN group/Data/Z/z3697/shot_data/SVS/z3697 svs4/' f1_id '.hdf'];
f1 = hdfread(f1,'fore');
f2_id = 'z3697svs4_543nm';
f2 = ['/Users/Rishabh/Dropbox (MIT)/PUFFIN group/Data/Z/z3697/shot_data/SVS/z3697 svs4/' f2_id '.hdf'];
f2 = hdfread(f2,'fore');
f3_id = 'z3697svs4_458nm';
f3 = ['/Users/Rishabh/Dropbox (MIT)/PUFFIN group/Data/Z/z3697/shot_data/SVS/z3697 svs4/' f3_id '.hdf'];
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

exportgraphics(gcf,[saveDir, 'SVS4_wl_calib.tiff']);

% Time calibration

temp = double(f2);

temp = temp(240,:);
[m,n] = size(temp);
[pks,idx] = findpeaks(temp,'MinPeakDistance',50,'MinPeakHeight',1e3);
figure
subplot(2,1,1);
plot(1:n,temp); hold on;
plot(idx,pks,'og'); 
xlabel('x (px)'); xlim([0,n]);
ylabel('counts (a.u.)');

subplot(2,1,2);
t = linspace(2900,2900+numel(pks)*20,numel(pks));
plot(t,idx,'o'); hold on;
P = polyfit(t,idx,1);
tq = 2900:3300;
plot(tq,polyval(P,tq),'-','color','r'); 
xlabel('t (ns)');
ylabel('x (px)');
exportgraphics(gcf,[saveDir, 'SVS4_time_calib.tiff']);

% time cal functions
t2px = @(t) polyval(P,t); % return px
px2t = @(x) (x - P(2)) / P(1); % return ns

% get wl calibration linear lam[nm] = m[nm/px]*y[px] + c[nm]
[m,n] = size(f3); 
ypx = 0:m; xpx = 0:n;
lam2 = 543; lam1 = 458; % nm
m = (lam2 - lam1) / (y2 - y1); % nm/px
c = lam1 - m* y1; % nm
px2wl = @(ypx) m * ypx + c;
wl2px = @(lam) (lam - c) / m;


% % plot
% fname = '/Users/Rishabh/Dropbox (MIT)/PUFFIN group/Data/Z/z3697/shot_data/SVS/z3697 svs4/z3697svs4_shot.hdf';
% data = hdfread(fname,'fore');
% 
% % Plot
% [m,n] = size(data);
% t = 1:n; wl = 1:n; % px
% t = px2t(t); wl = px2wl(wl);
% [tt,wl] = meshgrid(t,wl); 
% figure
% scatter(tt(:),wl(:),20,data(:),'s','filled'); colormap('hot');
% xlabel('t (ns)'); 
% ylabel('Wavelength (nm)'); 
% xlim([2900,3300]);
% ylim([300,850]);
% formatPlots();


%% lineouts
clc; close all;
% load shot data
fname = '/Users/Rishabh/Dropbox (MIT)/PUFFIN group/Data/Z/z3697/shot_data/SVS/z3697 svs4/z3697svs4_shot.hdf';
data = hdfread(fname,'fore');
data = double(data);

% Plot
f1 = figure();
[m,n] = size(data);
t = 1:n; wl = 1:n; % px
t = px2t(t); wl = px2wl(wl); % convert
[tt,wl] = meshgrid(t,wl); 
scatter(tt(:),wl(:),20,data(:),'s','filled','HandleVisibility','off'); colormap('hot'); hold on;
xlabel('t (ns)'); ylabel('Wavelength (nm)'); 
xlim([2900,3300]); ylim([250,850]);
formatPlots();



% lineouts
f2 = figure();
tid = [3000,3050,3100];
Lq = linspace(350,700,round(wl2px(700)-wl2px(350))); 
for ii = 1:numel(tid)
    tq = tid(ii) .* ones(size(Lq)); 
    tstp = 1; N = 5;
    out =  getavgLineout(tt,wl,data,tq,Lq,1,5);
    figure(f2);
    plot(Lq,out/3e3,'Linewidth',2,'DisplayName',[num2str(tid(ii)) ' ns'],'Color',sqclr('b',ii*2-1)); hold on; 
    figure(f1);
    plot(tq,Lq,'--',...
        'HandleVisibility','off','Color',sqclr('b',ii*2-1));
end

figure(f1);
set(gca,'TickDir','out');
exportgraphics(gcf,[saveDir, 'SVS4_calibrated.tiff']);

figure(f2);
xlabel('Wavelength (nm)'); ylabel('Counts (a.u.)'); 
formatPlots(900,16); grid on;

AlII = importfile('Al-II.csv', 7, 152);
AlII_wl = AlII(:,1); AlII_rel_int = AlII(:,7);
AlII_rel_int= AlII_rel_int(~isnan(AlII_wl));
AlII_wl= AlII_wl(~isnan(AlII_wl));
AlII_wl = AlII_wl(~isnan(AlII_rel_int));
AlII_rel_int= AlII_rel_int(~isnan(AlII_rel_int));
[AlII_rel_int, idx]= sort(AlII_rel_int);
AlII_wl = AlII_wl(idx);

figure(f2);
bar(AlII_wl,AlII_rel_int/2e8,5,'FaceColor','r','EdgeColor','r',...
    'DisplayName','Al-II','linewidth',2);

AlIII = importfile1('Al-III.csv', 7, 152);
AlIII_wl = AlIII(:,1); AlIII_rel_int = AlIII(:,8);
AlIII_rel_int= AlIII_rel_int(~isnan(AlIII_wl));
AlIII_wl= AlIII_wl(~isnan(AlIII_wl));
AlIII_wl = AlIII_wl(~isnan(AlIII_rel_int));
AlIII_rel_int= AlIII_rel_int(~isnan(AlIII_rel_int));
[AlIII_rel_int, idx]= sort(AlIII_rel_int);
AlIII_wl = AlIII_wl(idx);
[AlIII_wl,id] = unique(AlIII_wl);
AlIII_rel_int = AlIII_rel_int(id);

figure(f2);
bar(AlIII_wl,AlIII_rel_int/2e8,5,'FaceColor','g','EdgeColor','g',...
    'DisplayName','Al-III','linewidth',2);
set(gcf,'Position', [0  0   900   500]); ylim([0,1]);
exportgraphics(gcf,[saveDir, 'SVS4_lineouts.tiff']);

function out = getavgLineout(tt,wl,data,tq,Lq,tstp,N)
    out = 0;
    for ii = 0:N
        out = out + interp2(tt,wl,data,tq-N*tstp+ii*tstp,Lq); 
    end
    out = out / (2 * N);
end








