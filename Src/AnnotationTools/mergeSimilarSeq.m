%mergeSimilarSeq will merge sequences with similar sequences, but differing
%by either X nts or Y% difference. This is similar to removeDupSeq.m, but
%instead of looking for identical seq only, this one selectively removes
%sequences that do not change the overall lineage tree drastically.
%
%  VDJdata = mergeSimilarSeq(VDJdata, Map, DiffCutoff)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    Map: struct of indices of VDJdata columns
%    DiffCutoff: integer number of nts or fraction length of sequences for 
%      use in determining cutoff distance at which two sequences will not 
%      be grouped as one. Cutoff = 1 means nts with 1 nt difference will be
%      considered the same ONLY IF the lineage relations remains unchanged.
%
%  OUTPUT
%    VDJdata: main BRILIA data cell
%
function VDJdata = mergeSimilarSeq(varargin)
VDJdata = [];
varargin = cleanCommandLineInput(varargin{:}); %To enable post-processing capabilities

if numel(varargin) >= 3
    [VDJdata, Map, Cutoff] = deal(varargin{1:3});
    VDJdata = mergeSimilarSeq_Calc(VDJdata, Map, Cutoff);
else
    if numel(varargin) == 0 
        FileNames = getBriliaFiles('', true, true);
    elseif numel(varargin) >= 1
        FileNames = getBriliaFiles(varargin{1}, true, true);
    elseif strcmpi(varargin{1}, 'gui')
        FileNames = getBriliaFiles('', true, false);
    else
        FileNames = getBriliaFiles(varargin{1}, true, false);
    end
    if isempty(FileNames); return; end
    if numel(varargin) == 2
        Cutoff = varargin{2};
    else
        Cutoff = input('Merge sequences with differences greater than:\n 0 to 0.99 for seq length, or integer >= 1 for hamming distance: ');
    end
    for f = 1:length(FileNames)
        [VDJdata, VDJheader, ~, ~, Map] = openSeqData(FileNames{f});
        CurNum = size(VDJdata, 1);
        VDJdata = mergeSimilarSeq_Calc(VDJdata, Map, Cutoff);
        NewNum = size(VDJdata, 1);
        if CurNum ~= NewNum
            [Path, ~, Ext, Pre] = parseFileName(FileNames{f});
            Ptrn = ternary(Cutoff >= 1, '%s.Merge%d%s', '%s.Merge%0.2f%s');
            SaveName = fullfile(Path, sprintf(Ptrn, Pre, Cutoff, Ext));
            saveSeqData(SaveName, VDJdata, VDJheader);
            fprintf('%s: Finished with "%s".\n  Reduced %d similar sequences.\n', mfilename, FileNames{f}, CurNum - NewNum);
        else
            fprintf('%s: No sequences to merge in "%s".\n', mfilename, FileNames{f});
        end
    end
    VDJdata = []; %For multiple files, return nothing.
end

function VDJdata = mergeSimilarSeq_Calc(VDJdata, Map, Cutoff)
if isempty(VDJdata) || Cutoff <= 0; return; end %quick return
Map = getVDJmapper(Map); %In case user uses VDJheader instead of Map

IsSpliced = iscell(VDJdata{1});
if ~IsSpliced
    VDJdata = spliceData(VDJdata, Map); %Has to be done to ensure proper group numbers
end

parfor y = 1:length(VDJdata)
    VDJdata{y} = mergeSimilarSeqPerGroup(VDJdata{y}, Map, Cutoff);
end

if ~IsSpliced %Rejoin to ensure cluster-based splicing is preserved when passing this function
    VDJdata = joinData(VDJdata, Map); %Have to join to ensure cluster-based splicing w/ correct GrpNum
end

%CODINGNOTE: Test code for mergSimilarSeqPerGroup
% AncMap = [2 0 2 2; 4 2 2 2 ; 6 4 1 2; 8 4 1 2; 10 8 1 2];
% HSeq = {'ACGTACGTACGT'; 'ACGTCCCTACGT'; 'ACGGACCTCCCT'; 'TTTTACCTACGT'; 'TTCAACCTACGT'}
% 
% VDJdata = [num2cell(AncMap) HSeq];
% Map.SeqNum = 1;
% Map.ParNum = 2;
% Map.Ham = 3;
% Map.Template = 4;
% Map.lSeq = 0;
% Map.hSeq = 5;
% CutLen = 2;

function VDJdata = mergeSimilarSeqPerGroup(VDJdata, Map, Cutoff)
if size(VDJdata, 1) <= 1; return; end

AncMap = zeros(size(VDJdata, 1), 4);
AncMap(:, [1:2, 4]) = cell2mat(VDJdata(:, [Map.SeqNum, Map.ParNum, Map.Template]));
OrigSeqNum = AncMap(:, 1); %Need this later at end to preserve SeqNum!
AncMap = renumberAncMap(AncMap); %So you can use relative indexing

%CODING_NOTE: instead of redoing this comparison, consider saving info in VDJdata after lineage inference
SeqIdx = nonzeros([Map.hSeq; Map.lSeq]);
for j = 1:size(AncMap, 1)
    if AncMap(j, 2) > 0
        Seq1 = horzcat(VDJdata{AncMap(j, 1), SeqIdx}); %The children
        Seq2 = horzcat(VDJdata{AncMap(j, 2), SeqIdx}); %The parent
        AncMap(j, 3) = length(Seq1) - sum(cmprSeqMEX(Seq1, Seq2));
    end
end

if Cutoff >= 1 %Cutoff should be a integer
    CutLen = round(Cutoff);
else %Cutoff is a fraction of the sequence length
    CutLen = round(Cutoff * length(Seq1));
end

%Those that fall under the cutoff will be regrouped
ChildIdx = find(AncMap(:, 3) <= CutLen & AncMap(:, 2) > 0);
for j = 1:length(ChildIdx)
    ParNum = AncMap(ChildIdx(j), 2); %This is the sequence that will remain
    ParIdx = find(AncMap(:, 1) == ParNum);
    
    AncMap(ParIdx, 4) = AncMap(ParIdx, 4) + AncMap(ChildIdx(j), 4);
    AncMap(ChildIdx(j), 4) = 0;
    AncMap(ChildIdx(j), 1) = -ParNum; %use negative to prevent find it
    
    OldChildLoc = AncMap(:, 2) == ChildIdx(j);
    AncMap(OldChildLoc, 2) = ParNum;
end
AncMap(ChildIdx, :) = [];
VDJdata(ChildIdx, :) = [];

%Need to get the real sequence numbers
Idx = find(AncMap(:, 1:2) > 0);
AncMap(Idx) = OrigSeqNum(AncMap(Idx));
VDJdata(:, [Map.SeqNum; Map.ParNum; Map.Template]) = num2cell(AncMap(:, [1:2 4]));
% 
% %     
% %     SamePar = AncMap(ChildIdx(j), 1); %This is the sequence that will be lost
% %     MergeLoc = AncMap(:, 2) == SamePar; %These sequences must have new parents
% %     AncMap(MergeLoc, 2) = ParNum; %Reassigned parent to orphans
% %     AncMap(ParIdx, 4) = AncMap(ParIdx, 4) + AncMap(ChildIdx(j), 4); %Added up the templates
% % end
% % 
% 
% 
% %Those that fall under the cutoff will be regrouped
% Idx = find(AncMap(:, 3) <= CutLen & AncMap(:, 2) > 0);
% for j = 1:length(Idx)
%     MainPar = AncMap(Idx(j), 2); %This is the sequence that will remain
%     MainParIdx = find(AncMap(:, 1) == MainPar);
%     SamePar = AncMap(Idx(j), 1); %This is the sequence that will be lost
%     MergeLoc = AncMap(:, 2) == SamePar; %These sequences must have new parents
%     AncMap(MergeLoc, 2) = MainPar; %Reassigned parent to orphans
%     AncMap(MainParIdx, 4) = AncMap(MainParIdx, 4) + AncMap(Idx(j), 4); %Added up the templates
% end
% 
% AncMap(Idx, :)  = [];
% VDJdata(Idx, :) = [];
% VDJdata(:, [Map.SeqNum; Map.ParNum; Map.Template]) = num2cell(AncMap(:, [1:2 4]));