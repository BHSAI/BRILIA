%formatDRF will format the D gene reading frame data into either a Mx4 cell
%or 3*Mx2 cell. For the 2-col version, D gene names have ending of '-RF1',
%'-RF2', '-RF3'. For the 4-col version, D gene names do not have '-'RFN',
%and the frequency of each RF of D gene is placed in the respected columns
%2-4 for reading frames 1-3.
%
%  Out = formatDRF(DRF, ColumnFormat)
%
%  INPUT
%    DRF: Mx4 cell matrix of D gene (col1) and frequency of D gene in the
%      reading frame 1-3 (col2-4). 
%    ColumnFormat [2 or 4]: the 2- or 4-column format you want at the end  
%
%  OUTPUT
%    Out: 2-column format of 'DGENE-RFN' [Freq], OR 
%         4-column format of 'DGENE' [RF1_Freq] [RF2-Freq] [RF3Freq]

function Out = formatDRF(DRF, ColumnFormat)
%Determine # of col for output
if nargin == 1 
    if size(DRF, 2) == 4
        ColumnFormat = 2;
    elseif size(DRF, 2) == 2
        ColumnFormat = 4;
    else
        error('%s: The DRF input should be a 4- or 2- column cell array', mfilename);
    end
else
    if ColumnFormat ~= 2 && ColumnFormat ~= 4
        error('%s: The ColumnFormat input must be 2 or 4', mfilename);
    end
end

if size(DRF, 2) == 1 %This is the initial DRF name only, all with frequency 1
    DRF = [DRF num2cell(ones(size(DRF, 1), 1))];
end

if ColumnFormat == size(DRF, 2) %No changes needed
    Out = DRF;
    return
end

if ColumnFormat == 2
    Out = repelem({0}, size(DRF, 1)*3, 2);
    Out(:, 1) = repelem(DRF(:, 1), 3*ones(size(DRF, 1), 1));
    for k = 1:3:size(Out, 1)
        for j = 0:2
            Out{k+j, 1} = sprintf('%s-RF%d', Out{k+j}, j+1);
        end
    end
    Freq = cell2mat(DRF(:, 2:end))';
    Out(:, 2) = num2cell(Freq(:));
else
    if ~all(contains(DRF(:, 1), {'-RF1', '-RF2', '-RF3'}))
        error('%s: 1st column DRF must all have ''-RFN''', mfilename);
    end
    RF = cellfun(@(x) convStr2Num(x(end)), DRF(:, 1));
    GeneName = cellfun(@(x) x(1:end-4), DRF(:, 1), 'unif', false);
    [UnqGeneName, ~, UnqIdx] = unique(GeneName, 'stable');
    Out = repelem({0}, length(UnqGeneName), 4);
    Out(:, 1) = UnqGeneName;
    Out(sub2ind(size(Out), UnqIdx, RF+1)) = DRF(:, 2);
end