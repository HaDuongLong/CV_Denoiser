%% merge Q and P datasets
data_merge = [Q1n P1n];

%% Split predictor and target randomly into n parts
n_parts = 5;
data_parts = splitdata(data_merge, n_parts);
data_parts_Q = data_parts(:,1:size(Q1n,2),:);
data_parts_P = data_parts(:,size(Q1n,2)+1:size(Q1n,2)+size(P1n,2),:);

%% plot to check
idx3 = 1;
figure(1)
plot(data_parts_Q(idx3,:,idx3));
hold on
plot(data_parts_P(idx3,:,idx3));
hold off

%% return to original dimension for comparison
data_Q1n = permute(data_parts_Q,[1 3 2]);
data_Q1n = reshape(data_Q1n, [(n_parts)*size(data_parts_Q,1) size(data_parts_Q,2)]);

data_P1n = permute(data_parts_P,[1 3 2]);
data_P1n = reshape(data_P1n, [(n_parts)*size(data_parts_P,1) size(data_parts_P,2)]);

%% plot to check
idx4 = 1;
figure(2)
plot(data_Q1n(idx4,:));
hold on
plot(data_P1n(idx4,:));
hold off

%% K-fold MAPE
MAPE_each = zeros(1,n_parts);
MAPE_smooth = zeros(1,n_parts);
MAPE_smooth2 = zeros(1,n_parts);
for i = 1:n_parts
    validatePredictor = reshape(data_parts_Q(:,:,i), [size(data_parts_Q,1) size(data_parts_Q,2)]);
    validateTarget = reshape(data_parts_P(:,:,i), [size(data_parts_P,1) size(data_parts_P,2)]);

    trainPredictor= data_parts_Q;
    trainPredictor(:,:,i) = [];
    trainPredictor = permute(trainPredictor,[1 3 2]);
    trainPredictor = reshape(trainPredictor(:,:,:),[(n_parts-1)*size(data_parts_Q,1) size(data_parts_Q,2)]);


    trainTarget = data_parts_P;
    trainTarget(:,:,i) = [];
    trainTarget = permute(trainTarget,[1 3 2]);
    trainTarget = reshape(trainTarget(:,:,:),[(n_parts-1)*size(data_parts_P,1) size(data_parts_P,2)]);

    [numSegments, numFeatures]=size(P1n(1,:));
% ML training
layers = [
    featureInputLayer(numFeatures)
    fullyConnectedLayer(11200)
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(11200)
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(numFeatures)
    regressionLayer
    ];

miniBatchSize = 20;
options = trainingOptions("sgdm", ... %rmsprop sgdm adam
    MaxEpochs= 10, ...
    InitialLearnRate=1e-5,...
    MiniBatchSize=miniBatchSize, ...
    Shuffle="every-epoch", ...
    Plots="training-progress", ...
    Verbose=false, ...
    ValidationFrequency=floor(numFeatures/miniBatchSize), ...
    LearnRateSchedule="piecewise", ...
    LearnRateDropFactor=0.9, ...
    LearnRateDropPeriod=1, ...
    ValidationData={validatePredictor,validateTarget});


denoiseNetFullyConnected = trainNetwork(trainPredictor,trainTarget,layers,options);

% predict with trained ML
predictor_r=predict(denoiseNetFullyConnected, data_Q1n);
predictor_n1 = zeros(size(predictor_r,1),size(predictor_r,2));
    for j = 1:size(predictor,1)
        predictor_n1(j,:) = smooth(predictor_r(j,:), 0.023);
    end
predictor_n = smoothdata(predictor_r,2,"SmoothingFactor", 0.007);
MAPE_each(1,i) = mape(predictor_r, data_P1n, 'all');
MAPE_smooth(1,i) = mape(predictor_n, data_P1n, 'all');
MAPE_smooth2(1,i) = mape(predictor_n1, data_P1n, 'all');
end
%% avg MAPE
avg_mape_raw = mean(MAPE_each);
avg_mape_smooth = mean(MAPE_smooth);
avg_mape_smooth2 = mean(MAPE_smooth2);
%%
idx5 = 10900;
figure(3)
plot(data_Q1n(idx5,:));
hold on
plot(data_P1n(idx5,:),'k.');
plot(predictor_r(idx5,:));
plot(predictor_n(idx5,:));
plot(predictor_n1(idx5,:),'g');
hold off
%% Split function
function [data_parts] = splitdata(odata, n)
    data_shuffled = odata(randperm(size(odata,1)),:);
    part_size = floor(length(data_shuffled) / n);

    data_parts = zeros(part_size, size(odata,2), n);
    for i = 1:n
        start_index = (i-1) * part_size + 1;
        end_index = i * part_size;
        data_parts(:,:,i) = data_shuffled(start_index:end_index,:);
    end
end