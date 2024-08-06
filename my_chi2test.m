function [chi2stat, p, df] = my_chi2test(data)
% AK, written with chatGPT 23/06/2024

% Observed data
    observed = data;
    
    % Expected data
    [row_totals, col_totals, total] = marginals(data);
    expected = (row_totals * col_totals) / total;
    
    % Chi-square statistic calculation
    chi2stat = sum((observed(:) - expected(:)).^2 ./ expected(:));
    
    % Degrees of freedom
    [nRows, nCols] = size(data);
    df = (nRows - 1) * (nCols - 1);
    
    % p-value calculation
    p = 1 - chi2cdf(chi2stat, df);
end

function [row_totals, col_totals, total] = marginals(data)
    % Row totals
    row_totals = sum(data, 2);
    
    % Column totals
    col_totals = sum(data, 1);
    
    % Total
    total = sum(row_totals);
end
