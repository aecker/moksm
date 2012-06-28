units = 5;       % Number of units
D = 3;           % Dimensionality of fake features
mean_rate = 2;   % Mean firing rate (Hz)
T = 600000;      % Total time (seconds)
dr = 1e-2;       % Amount to (on average) move mean per second

features = zeros(0,D);
spike_times = [];
spike_ids = [];
for i = 1:units
    isi = [];
    while sum(isi) < T
        isi = [isi; exprnd(ones(1000,1) * 1 / mean_rate)];
    end
    t = cumsum(isi);
    t(t > T) = [];
    
    mu(i,:) = randn(1,D) * 5;
    mix(:,:,i) = randn(D,D);
    
    drift = [zeros(1,D); dr * cumsum(bsxfun(@times, randn(length(t)-1,D), diff(t)))];
    features = [features; bsxfun(@plus, (randn(length(t),D) + drift) * mix(:,:,i), mu(i,:))];
    spike_times = [spike_times; t];
    spike_ids = [spike_ids; i * ones(length(t),1)];
end

% Reorganize data by time
[spike_times, ids] = sort(spike_times);
features = features(ids,:);
spike_ids = spike_ids(ids,:);

m = MoKsm(features,spike_times);
m.params.dTmu = 100;
m = fitModel(m);
plot(m)

figure
subplot(121)
cla
hold on
subplot(122)
cla
hold on
c = hsv(units);
for i = 1:units
    subplot(121)
    idx = find(spike_ids == i);
    plot(spike_times(idx(1:10:end)), features(idx(1:10:end), 1), '.', 'MarkerSize',5, 'color', c(i,:));
    subplot(122)
    plot(features(idx(1:10:end), 1), features(idx(1:10:end), 2), '.', 'MarkerSize',5, 'color', c(i,:));
end
% Do afterwards so on top of all spikes
c = jet(units);
for i = 1:units
    subplot(121)
    plot(m.model.mu_t, squeeze(m.model.mu(1,:,i)), 'k','LineWidth',3,'color',c(i,:))
    subplot(122)
    plot(squeeze(m.model.mu(1,:,i)),squeeze(m.model.mu(2,:,i)), 'k','LineWidth',3,'color',c(i,:))
end

