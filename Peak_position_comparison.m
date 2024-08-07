%%
exp_peak_anal = peak_analysis(P1n,predictor1);
%exp_peak_anal = peak_analysis(cur_norm,cur_fil_t);
%exp_peak_anal = peak_analysis(cur_norm,cur_far_smooth');
%%
function errors_anal = peak_analysis(theoretical, predict)
[~, Ip] = max(theoretical, [], 2);
[~, Ipre] = max(predict, [], 2);

[~, ip] = min(theoretical, [], 2);
[~, ipre] = min(predict, [], 2);

po_ano = (Ip-800).*0.001 -0.4;
po_ano_p = (Ipre-800).*0.001 -0.4;

po_cat = 0.4 - ip.*0.001;
po_cat_p = 0.4 - ipre.*0.001;

mape_ano = mape(po_ano_p, po_ano);
mape_cat = mape(po_cat_p, po_cat);

mse_ano = mae(po_ano_p, po_ano);
mse_cat = mae(po_cat_p, po_cat);

errors_anal(1) = mape_ano;
errors_anal(2) = mape_cat;
errors_anal(3) = mse_ano;
errors_anal(4) = mse_cat;
disp(mape_ano)
disp(mape_cat)
disp(mse_ano)
disp(mse_cat)
end


% figure(1)
% plot(po_cat,'k.')
% hold on
% plot(po_cat_p,'r.')
% plot(po_ano_p,'r.')
% plot(po_ano,'k.')
% hold off