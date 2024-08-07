%% variable name
current;
potential;
formal_pot = 0.167;
over_pot_o = potential - formal_pot;
C1 = (1.3e-3)/1000; % mol/L to mol/cm^3

%% resampling CV data to fit deep learning
int_cv_pt = 3201; %1 1601 3201 
fin_cv_pt = 4800; %1600 3200 4800
cur_smooth = smoothdata(current, "SmoothingFactor", 0.0010);
cur_resamp = transpose(cur_smooth(int_cv_pt:fin_cv_pt)); %  current
pot_resamp = transpose(potential(int_cv_pt:fin_cv_pt)); % 
over_pot = pot_resamp - formal_pot;
%% normalization
scanr =  1; %scan rate: 0.02 0.05 0.1 0.2 0.5 1
cur_dens = (cur_resamp.*1000./(pi*(0.08^2)));
cur_dens_mean = cur_dens - mean(cur_dens);
cur_scanr = cur_dens./(sqrt(scanr)*C1); %./sqrt(scanr)
cur_norm = (cur_scanr-mean(cur_scanr))./(std(cur_scanr));

%% predict with model
cur_far = predict(denoiseNetFullyConnected,cur_norm);
cur_far_smooth = smooth(cur_far(1:end), 0.019);

%% plot
figure(12)
plot(over_pot,cur_far(1:end));
%% filter
Fs1 = 1000;

fs = 10;
cur_fil_far = lowpass(cur_far(10:end-10), 0.015, 'Steepness', 0.80);
m_start = mean(cur_far(1:25)); %normally use: cur_far_smooth, testing with cur_far
m_end = mean(cur_far(end-25:end));
cur_far_smooth_pad = [ones(25,1).*m_start; cur_far_smooth; ones(25,1).*m_end];
cur_fil1 = lowpass(cur_far_smooth_pad(1:end), 0.015, 'Steepness', 0.80);
cur_fil = cur_fil1(26:end-25);
cur_fil = smoothdata(cur_fil, "SmoothingFactor", 0.0004);
cur_fil_t = transpose(cur_fil);

%% calculate real current vs non norm predicted current
cur_real_smooth = transpose((cur_far_smooth.*(std(cur_scanr))+mean(cur_scanr)).*sqrt(scanr).*C1- mean(cur_dens));
%cur_real_smooth = transpose((cur_far_smooth.*(std(cur_scanr))+mean(cur_scanr)).*sqrt(scanr).*C1- mean(cur_dens));
cur_real_smooth_fil = transpose((cur_fil.*(std(cur_scanr))+mean(cur_scanr)).*sqrt(scanr).*C1- mean(cur_dens));
cur_real_fil = (cur_fil.*(std(cur_scanr))+mean(cur_scanr)).*sqrt(scanr).*C1- mean(cur_dens);
cur_real = (cur_far(1:end).*(std(cur_scanr))+mean(cur_scanr)).*sqrt(scanr);
cur_real_fil_far = ((cur_fil_far.*(std(cur_scanr))+mean(cur_scanr)).*sqrt(scanr).*C1- mean(cur_dens));
cur_real_smooth_1 = transpose((cur_far_smooth.*(std(cur_far))+mean(cur_far)).*sqrt(scanr));
%% calculate background current from predicted data
nonf_cur_s = cur_dens - cur_real_smooth;
nonf_cur_f = cur_dens - cur_real_smooth_fil;
nonf_cur_fil = cur_dens - cur_real_fil';
nonf_cur = cur_dens - cur_real;
nonf_cur_n = cur_norm - cur_far_smooth';
nonf_cur_smooth = lowpass(nonf_cur_s, 0.0001);


%%
idx2 =6180;
P1mean = (P1(idx2,:) - mean(Q1(idx2,:)))./1000.*(pi*(0.08^2));%- mean(Q1(idx2,:))
Q1mean = (Q1(idx2,:) - mean(Q1(idx2,:)))./1000.*(pi*(0.08^2));%- mean(Q1(idx2,:))
%B1n = Q1mean - P1mean;
B1n = Q1(idx2,:) - P1(idx2,:);

%% real current
cur_abs_sm = (cur_real_fil+mean(cur_dens))./1000.*(pi*(0.08^2));

figure(13)
plot(cur_abs_sm,'r.-'); %over_pot,
hold on
plot(cur_resamp,'k.-') %over_pot,
hold off
%% plot
figure(14)
plot(over_pot,cur_dens_mean,'k-');
hold on
plot(over_pot, cur_real_smooth, 'r-');
plot(over_pot, cur_real_fil, 'b-')
plot(over_pot, nonf_cur_fil,'-');
hold off

xlabel('Overpotential (V)',"FontSize",12)
ylabel('Current density (mA/cm^2)',"FontSize",12)
str = sprintf('Scan rate: %.3f V/s', scanr);
title(str);
%% plot 2
figure(15)

hold on
plot(over_pot, nonf_cur_fil ,'r.-');
hold off
xlabel('Overpotential (V)'), ylabel('Current density (mA/cm^2)')
str = sprintf('Scan rate: %.3f V/s', scanr);
title(str);

%% 
ipa = max(cur_real_smooth);
ipc = min(cur_real_smooth);
sqrt_scnr = sqrt(scanr);

% figure(16)
% plot(sqrt_scnr, ipa, 'ko', sqrt_scnr, ipc, 'ro');
% hold on;

%%
% idx2 = 14460;
figure(17)
plot(over_pot, cur_norm);
hold on 
plot(data_eta(idx2,:),Q1n(idx2,:),'k.');
plot(data_eta(idx2,:),P1n(idx2,:),'r.');
plot(over_pot, cur_far_smooth, 'g.');
hold off

%%
figure(18)

plot(data_eta(idx2,:),P1mean,'r.');
hold on
plot(data_eta(idx2,:),Q1mean,'k.');
plot(over_pot, cur_resamp,'b.');
plot(over_pot, cur_abs_sm,'g.');
hold off
xlabel('Overpotential (V)'), ylabel('Current density (mA/cm^2)')
str = sprintf('Scan rate: %.3f V/s', scanr);
title(str);

%%
cur_abs_fil_far = (cur_real_fil)./1000.*(pi*(0.08^2));

%%
% figure()
% plot(P1(idx2,:),'k','LineWidth',3)
% axis([-50 1650 -0.22 0.18]);

%% test
[ipa1, ipc1] = peakanalyzer(130,830, scanr, over_pot, cur_abs_fil_far);

%% function for peak
function [ipa, ipc] = peakanalyzer(c_idx, a_idx, scanr, pot, current) 

%c_idx = cathode baseline start data point index
%a_idx = anode baseline start data point index
%scanr = scan rate
%pot = potential

%finding cathode and anode peak current and data index position
[cat_min, cat_idx] = min(current);
[anode_max, anode_idx] = max(current);

%selecting data length for linear regression
length_cat = ceil((cat_idx - c_idx)/3);
length_anode = ceil((anode_idx - a_idx)/3);

%linear regression
basecath = polyfit(pot(c_idx:c_idx+length_cat-1),current(c_idx:c_idx+length_cat-1),1); 
baseanode = polyfit(pot(a_idx:a_idx+length_anode-1),current(a_idx:a_idx+length_anode-1),1);

cath_base = polyval(basecath, pot(c_idx:cat_idx));
anode_base = polyval(baseanode, pot(a_idx:anode_idx));

%peak current calculation
ipc = cat_min - cath_base(end);
ipa = anode_max - anode_base(end);

figure()
plot(pot,current);
hold on 
plot(pot(c_idx+1:c_idx+length(cath_base)),cath_base,'k');
plot(pot(a_idx+1:a_idx+length(anode_base)),anode_base,'k');

plot([pot(cat_idx) pot(cat_idx)], [cath_base(end) current(cat_idx)], 'k');
plot([pot(anode_idx) pot(anode_idx)], [anode_base(end) current(anode_idx)], 'k');

plot(pot(c_idx:c_idx+length_cat-1),cath_base(1:length_cat),'r','Linewidth',2);
plot(pot(a_idx:a_idx+length_anode-1),anode_base(1:length_anode),'r','Linewidth',2);
hold off

xlabel('Overpotential (V)'), ylabel('Current (A)')
str = sprintf('Scan rate: %.3f V/s', scanr);
title(str);
text(pot(cat_idx)*1.01, cat_min./2, sprintf('%.3e',ipc), 'HorizontalAlignment','right');
text(pot(anode_idx)*1.01, anode_max./2, sprintf('%.3e',ipa), 'HorizontalAlignment','left');

end