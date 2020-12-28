function labels = TrunkThis(features,thresholds,indices)


sel = (features(:,indices(1)) > thresholds(1) | ... %linearity 1
    features(:,indices(2)) > thresholds(2)) & ... %linearity 0.5
    features(:,indices(3)) > thresholds(3) & ... %verticality 0.5
    features(:,indices(4)) > thresholds(4); % density 0.25
labels(sel) = 1;
labels(~sel) = 2;
