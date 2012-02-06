function plotData(model)

[~, cl] = max(model.post);
K = size(model.post, 1);

c = combnk(1:3:12, 2);
colors = lines;
figure(101), clf
hdl = zeros(1, K);
for j = 1:size(c, 1)
    subplot(2, 3, j), hold all
    for k = 1:K
        plot(model.Ytrain(c(j, 1), cl == k), model.Ytrain(c(j, 2), cl == k), '.', 'markersize', 1, 'color', colors(k, :))
        hdl(k) = plot(model.mu(c(j, 1), :, k), model.mu(c(j, 2), :, k), '*-', 'color', colors(k, :));
    end
    axis tight
end
legend(hdl, arrayfun(@(x) sprintf('%d', x), 1:K, 'UniformOutput', false))


figure(102), clf
for k = 1:K
    subplot(ceil(K / 4), 4, k), hold on
    [~, chan] = max(mean(model.mu(:, :, k), 2));
    ndx = find(cl ~= k);
    plot(ndx, model.Ytrain(chan, ndx), '.', 'markersize', 1, 'color', 0.5 * ones(3, 1))
    ndx = find(cl == k);
    plot(ndx, model.Ytrain(chan, ndx), '.', 'markersize', 6, 'color', colors(k, :))
    plot(model.mu(chan, :, k), '-k', 'linewidth', 2)
    ylabel(sprintf('Ch %d', chan))
    axis tight
end

c = [2 3; 5 6; 8 9; 11 12];
 
 
