clc; close all; clear;

obj = SVSAnalysis() 

x = linspace(0,1,100);

y = @(x) -(x - 1) .* x;
y2 = y(0.8*x - 0.1);

figure
plot(x,y(x),DisplayName='target'); hold on;
plot(x,y2,DisplayName='Measured');

ylim([0,0.4]); legend();


[a,b,corr] = obj.find_shift_scale(x,y(x),y2)


plot(x,corr,'--',DisplayName='Corr.');
plot(x,obj.sweep_correct(x,y2,b,a),'--',DisplayName='Corr. 2');
 

