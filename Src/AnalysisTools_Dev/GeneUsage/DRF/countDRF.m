%countDRF will get the D gene reading frame frequency data.
%
%  S = countDRF(S, SizeFilter)
%
%  S = countDRF(S, SizeFilter, MinLen)
%
%  INPUT
%    S: structure of VDJdata files
%    SizeFilter: string code for filtering repertoire by clonotype sizes
%      see getGrpIdx.m
%    MinLen: minimum D length to include in data. Longer length helps to
%      remove short-D and/or poorly annoted D genes from data set. Default
%      0.
%
%  OUPUT
%    S: structure of VDJdata files, but with added field SizeFilter_DRF
%      .XC_DRF - 3*Nx2 cell, where column 1 is the D gene name plus -RF1, 
%        -RF2, or -RF3, and column 2 is the frequency.
%
%  See also formatDRF

function S = countDRF(S, SizeFilter, MinLen)
if nargin < 3
    MinLen = 0;
end

for f = 1:length(S)
    GrpData = getGrpIdx(S(f).VDJdata, S(f).VDJheader, SizeFilter);
    GrpIdx  = arrayfun(@(x) x.Idx(1), GrpData);
    Map = getVDJmapper(S(f).VDJheader);
    
    LenD = cell2mat(S(f).VDJdata(GrpIdx, Map. hLength(3)));    
    GrpIdx = GrpIdx(LenD > MinLen);
    
    CDR3S = cell2mat(S(f).VDJdata(GrpIdx, Map.hCDR3(3)));
    DelD5 = cell2mat(S(f).VDJdata(GrpIdx, Map.hDel(2)));
    LenVN = sum(cell2mat(S(f).VDJdata(GrpIdx, Map.hLength(1:2))), 2);

    Chart = [1 3 2;   %Each row is the mod(DDel5,3) +1
             2 1 3;   %Each col is the mod(LenVN-CDR3S+1,3) +1
             3 2 1];
    CharChart = '123';

    C = mod(LenVN-CDR3S+1, 3) + 1;
    R = mod(DelD5, 3) + 1;
    RF = Chart(sub2ind([3 3], R, C));
    RFChar = CharChart(RF);
    %To take advantage of countData.m, need to modify H-D_GeneName to
    %RF1:IGHD0-00*00 format
    DNames = S(f).VDJdata(GrpIdx, Map.hGeneName(2));
    for j = 1:length(DNames)
        if contains(DNames{j}, '|')
            DNames{j} = [RFChar(j) strrep(DNames{j}, '|', ['|' RFChar(j)])];
        else
            DNames{j} = [RFChar(j) DNames{j}];
        end
    end
    
    DNames = renameSameGene(DNames, S(f).Species);
    Data = countData(DNames);
    ReadF = cellfun(@(x) convStr2Num(x(1)), Data(:, 1));
    DGene = cellfun(@(x) x(2:end), Data(:, 1), 'unif', false);
    
    [UnqDGene, ~, UnqIdx] = unique(DGene);
    Freq = zeros(length(UnqDGene), 3);
    for k = 1:length(UnqIdx)
        R = UnqIdx(k);
        C = ReadF(k);
        Freq(R, C) = Freq(R, C) + Data{k, 2};
    end
    
    S(f).([SizeFilter '_DRF']) = formatDRF([UnqDGene num2cell(Freq)]);
end