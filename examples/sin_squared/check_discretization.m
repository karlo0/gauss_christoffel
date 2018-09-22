clear all; clf; 
BNAME = 'sin_squared'
EXTNAME = cell(1,3);
EXTNAME{1} = '_sample_1000.txt'
EXTNAME{2} = '_downsample_N100.txt'
EXTNAME{3} = '_downsample_N5.txt'

NAME = cell(1,3);
X = cell(1,length(NAME));
for k=1:length(NAME)
    NAME{k} = [BNAME EXTNAME{k}]; 
    X{k} = load(NAME{k});
end

f1 = figure(1);
hh = cell(1,length(X));
for k=1:length(X)
    hold on;
    hh{k} = plot(X{k}(:,1)', X{k}(:,2)');
    hold off;
end

set(hh{1}, 'Color', 'black');
set(hh{2}, 'Color', 'black', 'LineStyle', 'none', 'marker', '*');
set(hh{3}, 'Color', 'black', 'LineStyle', 'none', 'marker', 'x');

xlabel(gca, 'x & x-nodes')
ylabel(gca, 'f(x)& f-weights')

legend('f(x)', 'GC, N=30', 'GC, N=5', 'Location','northeast');

print(f1, [BNAME '.pdf'], '-dpdf');
close(f1);
