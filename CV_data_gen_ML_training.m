%% parameters and CV dataset generation
%% scanrate
%v(:,1) = logspace(-3 ,1, 50);
exponent = -4:1:1;
b = 10.^(exponent(:));
n_factor = 1:1:9;
v2 = (b.*(n_factor))';
v = reshape(v2,[numel(v2) 1]);

% k0
%k0_vary = logspace(-3, 1, 10); 
k0_vary = [1e-3, 5e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1, 1, 5, 10];
% CPE
CPE = linspace(0.6e-1, 1e-1, 5);

n_exp = linspace(0.8, 0.98, 5);
set_size = length(v)*length(k0_vary)*length(CPE)*length(n_exp);
% preallocation of data set
L = 1599;  % [=] number of iterations per t_k (pg 790). should be odd number
L1 = L + 1;
data_cur = zeros(set_size,L1);
data_cur_noise = zeros(set_size,L1);
data_eta = zeros(set_size,L1);
ft_i = zeros(set_size,L1);
ft_i_noise = zeros(set_size,L1);
ft_v = zeros(set_size,L1);
P1 = zeros(set_size,L1);
P1n = zeros(set_size,L1);
Q1 = zeros(set_size,L1);
Q1n = zeros(set_size,L1);

%
for i=1:length(v)
    v1 = v(i);
    for m = 1:length(k0_vary)
        for m1 = 1:length(CPE)
            for m2 = 1:length(n_exp)

% INDEPENDENT VARIABLES %%
C      = 1e-3;    % [=] mol/L, initial concentration of O. Default = 1.0
D      = 0.76E-5;   % [=] cm^2/s, O & R diffusion coefficient. Default = 1E-5
etai   = +0.4;   % [=] V, initial overpotential (relative to redox potential). Default = +0.2
etaf   = -0.4;   % [=] V, final overpotential (relative to redox potential). Default = -0.2
%v1      = v(i);   % [=] V/s, sweep rate. Default = 1E-3
n      = 1.0;    % [=] number of electrons transfered. Default = 1
alpha  = 0.5;    % [=] dimensionless charge-transfer coefficient. Default = 0.5
k0     = k0_vary(m);   % [=] cm/s, electrochemical rate constant. Default = 1E-2 
kc     = 1E-3;   % [=] 1/s, chemical rate constant. Default = 1E-3
T      = 298.15; % [=] K, temperature. Default = 298.15

%Cdl=1e-1; % F/cm2

% PHYSICAL CONSTANTS %%
F      = 96485;   % [=] C/mol, Faraday's constant
R      = 8.3145;  % [=] J/mol-K, ideal gas constant
f      = F/(R*T); % [=] 1/V, normalized Faraday's constant at room temperature

% SIMULATION VARIABLES %%
%L      = L;    % [=] number of iterations per t_k (pg 790). Default = 500
DM     = 0.45;   % [=] model diffusion coefficient (pg 788). Default = 0.45

% DERIVED CONSTANTS %%
tk  = 2*(etai-etaf)/v1;    % [=] s, characteristic exp. time (pg 790). In this case, total time of fwd and rev scans
Dt  = tk/L;               % [=] s, delta time (Eqn B.1.10, pg 790)
Dx  = sqrt(D*Dt/DM);      % [=] cm, delta x (Eqn B.1.13, pg 791)
j   = ceil(4.2*L^0.5)+5;  % number of boxes (pg 792-793). If L~200, j=65

% REVERSIBILITY PARAMETERS %%
ktk    = kc.*tk;             % dimensionless kinetic parameter (Eqn B.3.7, pg 797)
km     = ktk./L;              % normalized dimensionless kinetic parameter (see bottom of pg 797)
Lambda = k0/(D*f*v1)^0.5;     % dimensionless reversibility parameter (Eqn 6.4.4, pg. 236-239)

% CHEMICAL REVERSIBILITY WARNING %%
if km>0.1
    warning(['k_c*t_k/l equals ' num2str(km) ...
        ', which exceeds the upper limit of 0.1 (see B&F, pg 797)'])
end

% PRE-INITIALIZATION %%
C = C / 1000;           % Convert C from mol/L to mol/cm3
k = 0:L-1;                % time index vector
t = Dt * k;             % time vector
eta1 = etai - v1*t;      % overpotential vector, negative scan
eta2 = etaf + v1*t;      % overpotential vector, positive scan
eta = [eta1(eta1>etaf) eta2(eta2<=etai)]'; % overpotential scan, both directions
Enorm = eta*f;          % normalized overpotential
kf = k0.*exp(  -alpha *n*Enorm); % [=] cm/s, fwd rate constant (pg 799)
kb = k0.*exp((1-alpha)*n*Enorm); % [=] cm/s, rev rate constant (pg 799)

O = C*ones(L+1,j); % [=] mol/cm^3, concentration of O
R = C*zeros(L+1,j);  % [=] mol/cm^3, concentration of R
JO = zeros(1,L+1); % [=] mol/cm^2-s, flux of O at the surface

% START SIMULATION %%
% i1 = time index. i2 = distance index
for i1 = 1:L
    % Update bulk concentrations of O and R
    for i2 = 2:j-1
        O(i1+1,i2) = O(i1,i2) + DM*(O(i1,i2+1)+O(i1,i2-1)-2*O(i1,i2));

        R(i1+1,i2) = R(i1,i2) + DM*(R(i1,i2+1)+R(i1,i2-1)-2*R(i1,i2)) ...
            - km * R(i1,i2);
    end

    % Update flux
    JO(i1+1)   = ( kf(i1+1).*O(i1+1,2) - kb(i1+1).*R(i1+1,2) ) ./ (1 + Dx/D*(kf(i1+1) + kb(i1+1)) );

    % Update surface concentrations
    O(i1+1,1) = O(i1+1,2) - JO(i1+1)*(Dx/D);
    R(i1+1,1) = R(i1+1,2) + JO(i1+1)*(Dx/D) - km*R(i1+1,1);
end

% Calculate current density, Z, from flux of O
Z = -n.*F.*JO * 1000; % [=] A/cm^2 -> mA/cm^2, current density

% charging current
Fs=v1*L;  %sampling freq v(i) is the scan rate
T=1/Fs;

tt=0:1/Fs:1/(v1);    %creating array for time axis
tt_tr=-1/Fs:1/Fs:0.1-1/Fs;
L1=length(tt);

tx=-0.5*sawtooth(2*pi*v1*tt, 1/2); %triangular wave of CV

y=fft(tx);
ys=fftshift(y);
amp=abs(ys);



ly=length(y);
f=(-ly/2:ly/2-1)/ly*Fs;

% setting threshold for phase
th = 1e-6;

% current phasor

nf=max(f)-abs(f)+v1;
Qval = CPE(m1); %pseudo-capacitance value f/cm2
n_cpe = n_exp(m2);  %CPE exponent factor 0< n_cpe <= 1 (1 ideal capacitance) 
xc1= +1./(Qval*((1i*2*pi.*nf).^n_cpe)); %impedance of CPE 1/(Q*(wj)^n)
xc=1./xc1;

m_xc1=abs(xc1);
a_xc1=angle(xc1);

zz=y .* xc;

%inverse fft to obtain current in time domain
ii=ifft(zz, 'symmetric'); %final charging current with CPE

% Sometimes length(eta) = length(Z) + 1. If this is the case, truncate last value
if length(eta) > length(Z)
    eta = eta(1:end-1);
end
loop_index = length(k0_vary)*length(CPE)*length(n_exp)*(i-1)+(m-1)*length(CPE)*length(n_exp)+(m1-1)*length(n_exp)+m2;

data_eta(loop_index,:)=eta; %letdo
data_cur(loop_index,:)=Z;   %letdo
data_cur_noise(loop_index,:)=Z+ii;

% normalization
P1(loop_index,:)=data_cur(loop_index,:);
Q1(loop_index,:)=data_cur_noise(loop_index,:);


P1t(loop_index,:)=P1(loop_index,:)./(sqrt(v1).*C); %./sqrt(v1); .*1.025
Q1t(loop_index,:)=Q1(loop_index,:)./(sqrt(v1).*C); %./sqrt(v1);

P1n(loop_index,:)=(P1t(loop_index,:) -mean(Q1t(loop_index,:)))/std(Q1t(loop_index,:)); %);
Q1n(loop_index,:)=(Q1t(loop_index,:) -mean(Q1t(loop_index,:)))/std(Q1t(loop_index,:)); %);


P1mean(loop_index)=mean(P1(loop_index,:));
Q1mean(loop_index)=mean(Q1(loop_index,:));
P1std(loop_index)=std(P1(loop_index,:));
Q1std(loop_index)=std(Q1(loop_index,:));



app_fre=1/t(end);
Fs=1/Dt;
f= Fs*(0:(L/2))/L;
fre(i,:)=f;

            end
        end
    end
end

%% ML training model
target=zeros(size(P1n));

title("Single-Sided Amplitude Spectrum of X(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")


%-----------neural network
[trainInd,valInd,testInd] = dividerand(set_size,0.8,0.2,0);

trainPredictor=Q1n(trainInd,:);
trainTarget=P1n(trainInd,:);

validatePredictor=Q1n(valInd,:);
validateTarget=P1n(valInd,:);

[numSegments, numFeatures]=size(P1n(1,:));

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


numWeights = 0;
for index = 1:numel(denoiseNetFullyConnected.Layers)
    if isa(denoiseNetFullyConnected.Layers(index),"nnet.cnn.layer.FullyConnectedLayer")
        numWeights = numWeights + numel(denoiseNetFullyConnected.Layers(index).Weights);
    end
end
disp("Number of weights = " + numWeights);
%% save trained deep learning
% CVtest20230608 = denoiseNetFullyConnected;
% save CVtest20230608;
%% prediction
predictor=predict(denoiseNetFullyConnected, Q1n);
%%

temp=predictor.*Q1std(1)+Q1mean(1);
predictor1 = zeros(size(predictor));
[row_no col_no] = size(predictor);

%% data smoothing
predictor1 = smoothdata(predictor,2, "SmoothingFactor", 0.005);

%% check predictor results from deep learning network
idx= 13500;
figure(6)
clf;
plot(data_eta(idx,:), (predictor1(idx,:)*std(Q1(idx,:))+mean(Q1(idx,:))),'k.-'); %*std(P1(idx,:))+mean(P1(idx,:)))
hold on;
plot(data_eta(idx,:), P1(idx,:),'r.-');
hold on;
plot(data_eta(idx,:), Q1(idx,:),'b.-');
xlabel('Overpotential (V)'), ylabel('Current density (mA/cm^2)')
hold off
%str_1 = sprintf('Scan rate: %.3f V/s', scanr);
%title('Scan rate: 0.2 V/s')

P1save = P1(idx,:);
Q1save = Q1(idx,:);
predict1save = (predictor1(idx,:)*std(Q1(idx,:))+mean(Q1(idx,:)));
%%
figure(7);
plot(data_eta(idx,:), P1n(idx,:),'r.');
hold on
plot(over_pot, cur_far_smooth,'k.')
plot(over_pot, cur_norm,'g.');
plot(data_eta(idx,:), Q1n(idx,:),'b.');
hold off
xlabel('Overpotential (V)'), ylabel('Current density (mA/cm^2)')
title('Scan rate: 1 V/s')
%%
figure(8)
plot(data_eta(idx,:), P1n(idx,:), 'k.');
hold on;
plot(data_eta(idx,:), Q1n(idx,:),'r.');
plot(over_pot, cur_far_smooth,'b.');
hold off;

%% 
figure(9)
hold on
plot(cur_real_smooth,'ko');
plot(cur_dens, 'bo')
hold off