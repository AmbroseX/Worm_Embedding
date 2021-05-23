answer=inputdlg(strcat('请输入',filepath,':',wormName,'的最佳Kstar'),'The best Kstar');
Kstar=str2num(answer{1});
disp(['Best Kstar is:',num2str(Kstar)]);


y=load(fullfile(savefolder,strcat(filepath,'_',wormName,'_','angle_PCA_Kmode.mat')));
y=y.y;  % M*Kmode的 PCA  angle_data


%% optimize K
%params
%Kmax, maximum K used in the search. As a first guess, autocorrelation
%time, zero of autocorrelation function,oscillation period or its double.
%Delrange can be used to sample K coarsely, recommended when starting.
maxdel = 80; delrange = 1:1:maxdel;
%number of random points used for prediction. Larger number of points give
%a better estimate, but takes longer.
npred = 7000; %7000
%tau_max. Harder to estimate, set it to several multiples of the
%autocorrelation function. Check and error curves in rmsp_d and verify that
%the error has saturated. If it never saturates then most likely the
%nonstationarity is too strong. Contact tosifahamed@gmail.com in this case
Tmax = 700;

%generate random test indices for the search
[ybar, ~, idx, ~] = embed_cells(y, maxdel); %embed for maximum delay to calculate the valid range of test indices
%center ybar by removing mean
ybar = bsxfun(@minus, ybar, mean(ybar));
nb = find_nb(ybar); %find nearest neighbors and use the same for the entire search
%getting valid test indices
idx = idx(1:end - Tmax, 1);
testinds = idx(randi(length(idx), npred * 2, 1));


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
figure('Name',strcat(filepath,'_',wormName),'NumberTitle','off');
hold on
errorbar(1:mmax, tpred_m(:, 1), tpred_m(:, 1) - tpred_m(:, 2), tpred_m(:, 3) - tpred_m(:, 1), 'ko', 'markersize', 8);
xlabel('m')
ylabel('T_{pred}')
hold off

