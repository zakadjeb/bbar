function [stat_obs, p_perm, stat_perm_dist] = bbar_ersp_permutation(data1, data2, n_perm, test_type, paired)
%--------------------------------------------------------------------------
% Pixel-based permutation test for 2D ERSP data
%
% INPUTS:
%   data1, data2 : [subjects x freq x time] (paired)
%                  or [subjects1 x freq x time], [subjects2 x freq x time]
%   n_perm        : number of permutations (default = 1000)
%   test_type     : 't'  or  'wilcoxon'
%   paired        : true  -> dependent test
%                   false -> independent test
%
% OUTPUTS:
%   stat_obs       : observed t- or z-stat map
%   p_perm         : permutation-based p-value map (FWER corrected)
%   stat_perm_dist : null distribution of max|stat|
%
% Author: Zakaria Djebbara 
%--------------------------------------------------------------------------

if nargin < 3, n_perm = 1000; end
if nargin < 4, test_type = 'wilcoxon'; end
if nargin < 5, paired = true; end

% --- Basic sanity checks ---
assert(ndims(data1) == 3 && ndims(data2) == 3, 'Data must be 3D arrays.');
[n1, nFreq, nTime] = size(data1);
[n2, ~, ~] = size(data2);
if paired, assert(n1 == n2, 'For paired design, sample sizes must match.'); end

fprintf('\nRunning %s-test (%s design) with %d permutations...\n', ...
        test_type, ternary(paired, 'paired', 'independent'), n_perm);

% --- Compute observed statistic map ---
stat_obs = nan(nFreq, nTime);

for f = 1:nFreq
    for t = 1:nTime
        x = squeeze(data1(:, f, t));
        y = squeeze(data2(:, f, t));
        switch lower(test_type)
            case 't'
                if paired
                    [~,~,~,stats] = ttest(x, y);
                else
                    [~,~,~,stats] = ttest2(x, y);
                end
                stat_obs(f, t) = stats.tstat;

            case 'wilcoxon'
                if paired
                    [~,~,stats] = signrank(x, y);
                else
                    [~,~,stats] = ranksum(x, y);
                end
                if isfield(stats, 'zval')
                    stat_obs(f, t) = stats.zval;
                else
                    stat_obs(f, t) = 0;
                end

            otherwise
                error('Unknown test_type: use "t" or "wilcoxon"');
        end
    end
end

% --- Initialize null distribution ---
stat_perm_dist = nan(n_perm,1);
stat_perm = zeros(nFreq, nTime);

% --- Permutation loop ---
for perm_i = 1:n_perm

    if paired
        % Within-subject label flips
        flip_sign = (rand(n1,1) > 0.5);
        perm_data1 = data1;
        perm_data2 = data2;
        for s = 1:n1
            if flip_sign(s)
                tmp = perm_data1(s,:,:);
                perm_data1(s,:,:) = perm_data2(s,:,:);
                perm_data2(s,:,:) = tmp;
            end
        end
    else
        % Shuffle subject labels across groups
        all_data = cat(1, data1, data2);
        perm_idx = randperm(size(all_data,1));
        perm_data1 = all_data(perm_idx(1:n1),:,:);
        perm_data2 = all_data(perm_idx(n1+1:end),:,:);
    end

    % Compute permutation statistic map
    for f = 1:nFreq
        for t = 1:nTime
            x = squeeze(perm_data1(:, f, t));
            y = squeeze(perm_data2(:, f, t));
            switch lower(test_type)
                case 't'
                    if paired
                        [~,~,~,stats] = ttest(x, y);
                    else
                        [~,~,~,stats] = ttest2(x, y);
                    end
                    stat_perm(f, t) = stats.tstat;

                case 'wilcoxon'
                    if paired
                        [~,~,stats] = signrank(x, y);
                    else
                        [~,~,stats] = ranksum(x, y);
                    end
                    if isfield(stats, 'zval')
                        stat_perm(f, t) = stats.zval;
                    else
                        stat_perm(f, t) = 0;
                    end
            end
        end
    end

    % Store max absolute statistic (FWER correction)
    stat_perm_dist(perm_i) = max(abs(stat_perm(:)));

    if mod(perm_i, round(n_perm/10)) == 0
        fprintf('Permutation %d/%d complete\n', perm_i, n_perm);
    end
end

% --- Compute corrected p-values ---
p_perm = nan(nFreq, nTime);
for f = 1:nFreq
    for t = 1:nTime
        p_perm(f,t) = mean(abs(stat_perm_dist) >= abs(stat_obs(f,t)));
    end
end

fprintf('Done.\n');
end

% --- Small helper: ternary operator ---
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
