%countSHMfreq will return the mutation rates of CT and AG nts, in terms of % of
%total CT or AG nts in the sequence. This is done to ensure that mutations
%account for availability of nts that can be mutated, rather than total
%mutation counts.
%
%  S = countSHMfreq(S, SizeFilter)
%
%  INPUT
%
%  OUTPUT
%     S: input S with new data as such '
%       .G2C - germline-to-child mutation count
%       .P2C - parent-to-child mutation count, EXCLUDING 1st sequence,
%              since that's G2C and cannot guarantee P2C.
%         .Header - {'A' 'C' 'G' 'T' '}
%         .Data

function S = countSHMfreq(S, SizeFilter)
for f = 1:length(S)
    G = getGrpIdx(S(f).VDJdata, S(f).VDJheader, SizeFilter);
    Idx = arrayfun(@(x) x.Idx(1), G);
    Map = getVDJmapper(S(f).VDJheader);
    
    %Extract the SHM column index, from Map.hVmut, .hMmut, etc...
    MapCell = struct2cell(Map);
    MutFieldLoc = contains(fieldnames(Map), 'mut');
    MutLoc = vertcat(MapCell{MutFieldLoc});
    MutLoc(MutLoc == 0) = [];

    %Extract the VMDNJ-VJ nt lengths
    LengthFieldLoc = contains(fieldnames(Map), 'Length');
    LengthLoc = vertcat(MapCell{LengthFieldLoc});
    LengthLoc(LengthLoc == 0) = [];
    
    %Computer %SHM of length
    FreqSHM = sum(cell2mat(S(f).VDJdata(Idx, MutLoc)), 2) ./ sum(cell2mat(S(f).VDJdata(Idx, LengthLoc)), 2);
    S(f).([SizeFilter '_SHM_G2C']) = countData(num2cell(FreqSHM));
    
    %PERFORM FOR P2C
    IdxNoFirst = arrayfun(@(x) x.Idx(2:end), G, 'unif', false);
    IdxNoFirst(cellfun('isempty', IdxNoFirst)) = [];
    Idx = horzcat(IdxNoFirst{:});
    %Idx = vertcat(G.Idx);
    FreqSHM = zeros(size(Idx));
    for k = 1:length(Idx)
        MissLoc = ~cmprSeqMEX(S(f).VDJdata{Idx(k), Map.hRefSeq}, S(f).VDJdata{Idx(k), Map.hSeq}, 'n');
        FreqSHM(k) = sum(MissLoc) / length(MissLoc);
    end
    if isempty(FreqSHM)
        S(f).([SizeFilter '_SHM_P2C']) = {0 0};
    else
        S(f).([SizeFilter '_SHM_P2C']) = countData(num2cell(FreqSHM));
    end
end
