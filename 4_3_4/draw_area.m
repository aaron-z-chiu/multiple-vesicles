close all; clear ; clc;

y1 = load('surface_areac1.m');                     

dt = 2e-5;
nPoints = length(y1);
t = (0:nPoints-1)*dt;
figure;
plot(t, y1, 'b-', ...
    'Marker', 'o', ...
    'MarkerFaceColor', 'b', ...   
    'MarkerEdgeColor', 'b', ...   
    'LineWidth', 1.2, ... 
    'MarkerIndices', 1:1000:length(y1));
hold on;
y2 = load('surface_areac2.m');
plot(t, y2, 'r-', ...
    'Marker', '^', ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'r', ...
    'LineWidth', 1.2, ... 
    'MarkerIndices', 1:1000:length(y2));
hold off;

xlabel('Time (s)');
ylabel('Surface Area');
grid on;

legend("Vesicle 1","Vesicle 2");

set(gca,'fontsize',14)
text('Interpreter','latex','String','$x$','Position',[0.891962081682177 -0.0577590895778421 7.105427357601e-15],'FontSize',19);
text('Interpreter','latex','String','$y$','Position',[-0.0749639883567335 0.915105165952411 7.105427357601e-15],'FontSize',19);

tt = sprintf('area.eps');
print('-deps', tt);