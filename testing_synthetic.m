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
m.params.driftRate = 2e-4;
m = fitModel(m);
plot(m)

