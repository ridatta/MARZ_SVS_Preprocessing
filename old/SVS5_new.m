% SVS5 Processing

clc; close all; clear;

% (1) Show raw speactral Image

inum = 5;
fid = 'z3697svs5_shot';
inDir = checkDir(['\Users\rdatta\Dropbox (MIT)\PUFFIN\Data\MARZ\SVS/z3697 svs' num2str(inum) '/']);
saveDir = checkDir('/Users/Rishabh/Dropbox (MIT)/PUFFIN/Data/MARZ/SVS/');
data = hdfread([inDir, fid, '.hdf'],'fore');

figure
imagesc(data);
set(gca,'Ydir','Normal');
colormap(hot)
title(strrep(fid,'_','\_'));

ylabel('wavelength (a.u.)');
xlabel('time (a.u.)');
set(gcf,'Position', [262   118   820   556]);

exportgraphics(gcf,[saveDir, fid '_raw.tiff']);

%% Alignmnet test and rectification

