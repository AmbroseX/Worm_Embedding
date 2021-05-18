answer=inputdlg('请输入最佳Kstar','The best Kstar');
Kstar=str2num(answer{1});
disp(['Best Kstar is:',num2str(Kstar)]);

y=load(fullfile(savefolder,strcat(wormName,'_','angle_PCA.mat')));
y=y.y;

%% Use SVD to optimize m
%set the value of Kstar as determined from Tpred(K) above
%Maximum embedding dimension corresponds to the maximum number of singular
%vectors to consider.
[ybar, ybar_orig, ybar_orig_idx] = embed_cells(y, Kstar);
[~, tstidx] = ismember(testinds, ybar_orig_idx);
%perform SVD，奇异值分解寻找最佳的m
[U, S, V] = svd(ybar, 'econ');

mmax = 14; %m最大到多少为止
for mi = 1:mmax
    [rmsp_m{mi}, es_m(mi, :), ar_m(mi, :), tpred_m(mi, :), ts_m(mi, :)] = predict_nn(U(:, 1:mi), ybar_orig, Tmax, npred, V(:, 1:mi) * S(1:mi, 1:mi), nb, tstidx);
    disp(mi);
end

save(fullfile(savefolder,strcat(filepath,'_',wormName,'_','angle_PCA_Tpredm.mat')),'tpred_m'); % %存储矩阵
% plot Tpred(m)
figure
hold on
errorbar(1:mmax, tpred_m(:, 1), tpred_m(:, 1) - tpred_m(:, 2), tpred_m(:, 3) - tpred_m(:, 1), 'ko', 'markersize', 8);
xlabel('m')
ylabel('T_{pred}')
hold off

saveas(gcf, fullfile(savefolder,'Tpred(m).jpg'));