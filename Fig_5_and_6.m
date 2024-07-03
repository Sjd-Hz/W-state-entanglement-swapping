clear;
clc;
close all;

[r, q] = meshgrid(0:0.01:1, 0:0.01:1);

% Compute fidw and fider
fidw = 1./(1 + r.*(1 - q));
fider = 1./(1 + r);
p=((1-q).^2.*(r-q.*r+1))./(1+r); %Weak measurements success probability
%p=((1-q).^2.*(r-q.*r+1)).*(1-r)./2; %Entanglement swapping success probability


figure(1);
% Set font and color settings
set(gcf, 'DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Times New Roman');
surf(r, q, fidw, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
hold on;
surf(r, q, fider, 'EdgeColor', 'none', 'FaceAlpha', 1,  'FaceColor', [0.5, 0.5, 0.5]);
xlabel('r', 'FontWeight', 'bold');
ylabel('q', 'FontWeight', 'bold');
zlabel('Fidelity', 'FontWeight', 'bold');
axis tight;


figure(2);
set(gcf, 'DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Times New Roman');
surf(r, q, p, 'EdgeColor', 'none');
colormap("cool"); 
xlabel('r', 'FontWeight', 'bold');
ylabel('q', 'FontWeight', 'bold');
zlabel('Success probability', 'FontWeight', 'bold');
%zlabel('Entanglement swapping success probability', 'FontWeight', 'bold');
axis tight;