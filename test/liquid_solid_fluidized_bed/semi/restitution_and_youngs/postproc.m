restitution = 0.05:0.05:0.95;
y = -log(restitution) ./ sqrt(log(restitution) .* log(restitution) .+ pi * pi);

linienstaerke = 1;
MarkerGroesse = 6;

figure(2)
h = plot(restitution, y, '-');
set(h,'LineWidth', linienstaerke, 'MarkerSize', MarkerGroesse);
set(gca, 'FontSize', 12)
set(gca, 'FontSize', 12)
print('restitution.png')
