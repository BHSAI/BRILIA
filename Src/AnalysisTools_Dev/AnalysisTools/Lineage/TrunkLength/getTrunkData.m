%getTreeSpecs will get the lineage tree summary data, such as tree trunk
%lengths, heights, nodes, and leaves. 
%
%  S = getTrunkData(S)
%
%  INPUT
%    S: structure of BRILIA output files
%
%  OUTPUT
%    S: input with added fields 
%      .TreeTrunks = germline to 1st ancestor hamming %
%      .TreeHeights = ancestor to furthest descendant hamming %
%      .TreeNodes = number of nodes per clonotype
%      .TreeLeaves = number of leaves per clonotype
%
function S = getTreeSpecs(S, varargin)
if ~isfield(S, 'TreeSpecs') 
    S(1).TreeSpecs = struct('Trunks', [], 'Heights', [], 'Nodes', [], 'Leaves', []);
end
for f = 1:length(S)
    Map = getVDJmapper(S(f).VDJheader);
    SeqIdx = nonzeros(vertcat(Map.hSeq, Map.lSeq));
    
    FieldNames = fieldnames(Map);
    MutLoc = endsWith(FieldNames, 'mut');
    MutNames = FieldNames(MutLoc);
    MutIdx = nonzeros(cellfun(@(x) Map.(x), MutNames));
    
    G = getGrpIdx(S(f).VDJdata, S(f).VDJheader, varargin{:});
    
    %Initialize frequency data
    S(f).TreeSpecs.Trunks  = zeros(length(G), 1);
    S(f).TreeSpecs.Heights = zeros(length(G), 1);
    S(f).TreeSpecs.Nodes   = zeros(length(G), 1); 
    S(f).TreeSpecs.Leaves  = zeros(length(G), 1);

    %Compute each groups' tree specs
    for j = 1:length(G)
        AllSHM = sum(cell2mat(S(f).VDJdata(G(j).Idx, MutIdx)), 2);
        SeqLen = sum(cellfun('length', S(f).VDJdata(G(j).Idx(1), SeqIdx)), 2);
        TrunkData(j, 1) = AllSHM(1) / SeqLen;
        TrunkData(j, 2) = (max(AllSHM) - AllSHM(1)) / SeqLen;
    end
    S(f).TrunkData = TrunkData;
end
    