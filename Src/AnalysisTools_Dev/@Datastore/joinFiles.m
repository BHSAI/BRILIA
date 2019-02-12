function O = joinFiles(O, Option, Level)
if isempty(O.SaveDir)
    warning('%s: Must specify save dir before hand.', mfilename);
    return
end
if nargin < 2
    Option = 'abridged';
end
if  nargin < 3
    Level = 'Clone';
else
    if ~ismember(lower(Level), {'clonotype', 'clone'})
        error('%s: Level "%s" is not valid. Either ''clonotype'' or ''clone''.', mfilename, Level); 
    end
end
SaveFile = fullfile(O.SaveDir, sprintf('ANNOT.%s.%s.BRILIA%s.csv', Option, Level, BRILIA('version')));

%Abridged : removes GeneNum and Align from the data, before joining
%Addon : adds clonal count data to the data before saving
%AddID : adds on 2 columns: group name and also file name
% ValidOptions = {'abridged', 'addon'};
%Option = intersect(lower(Option), ValidOptions);

if startsWith('abridged', Option, 'ignorecase', true)
    Map = O.Map;
    Map.lValign = 0;
    Map.lJalign = 0;
    Map.hValign = 0;
    Map.hDalign = 0;
    Map.hJalign = 0;
    Map.hGeneNum = 0;
    Map.lGeneNum = 0;
    Map.Chain = 0;
    MapIdx = struct2cell(Map);
    MapIdx = nonzeros(vertcat(MapIdx{:}));
    O.SelectedVariableNames = MapIdx;
    VDJheader = O.VDJheader(MapIdx);
else
    VDJheader = O.VDJheader;
end

O.reset;
for f = 1:numel(O.Files)
    f
    Ptr = 'Grp(\w+)-';
    [~, FileName] = parseFileName(O.Files{f});
    Token = regexp(FileName, Ptr, 'tokens');
    if isempty(Token)
        Token = FileName; %use file name default
    else
        Token = Token{1}{1}; %unwrap
    end
       
    [VDJdata, Map] = read(O);
    if strcmpi('clonotype', Level)
        [VDJdata, AddClonotypeHeader, AddClonotypeData] = condenseVDJdata(VDJdata, Map, O.Idx{f});
    else
        AddClonotypeHeader = {};
        AddClonotypeData = {};
    end
    
    N = size(VDJdata, 1);
    TC = cell2mat(VDJdata(:, Map.Template));
    CellFr = TC/sum(TC); %Seq_Freq_Of_Functional_Seq
    
    MutIdx = nonzeros([Map.hVmut, Map.hMmut, Map.hDmut, Map.hNmut, Map.hJmut Map.lVmut, Map.lNmut, Map.lJmut]);
    TotCtSHM = sum(cell2mat(VDJdata(:, MutIdx)), 2); 
    
    TotLen = sum(cell2mat(VDJdata(:, Map.hLength)), 2);
    TotFrSHM = TotCtSHM./TotLen;
    
    Group = repelem(O.GroupName(f), N, 1);
    SampleID = repelem({Token}, N, 1);
    
    AddHeader = horzcat({'Group', 'SampleID', 'CellFraction', 'TotalCountOfSHM', 'TotalFractionOfSHM'}, AddClonotypeHeader{:});
    AddData = horzcat(Group, SampleID, num2cell([CellFr, TotCtSHM, TotFrSHM]), AddClonotypeData);
    
    if f == 1
        writeDlmFile([[AddHeader VDJheader]; [AddData VDJdata]], SaveFile);
    else
        writeDlmFile([AddData VDJdata], SaveFile, 'append');
    end
end

function [VDJdata, AddHeader, AddData] = condenseVDJdata(VDJdata, Map, IdxC)
VDJdata = spliceData(VDJdata, Map, IdxC);
N = numel(VDJdata);
TemplateIdx = Map.Template;
MutIdx = nonzeros([Map.hVmut, Map.hMmut, Map.hDmut, Map.hNmut, Map.hJmut Map.lVmut, Map.lNmut, Map.lJmut]);
TreeTrunk = zeros(N, 1);
TreeHeight = zeros(N, 1);
parfor j = 1:N
    VDJdata{j}{1, TemplateIdx} = sum(cell2mat(VDJdata{j}(:, TemplateIdx)));
    TotCtSHM = sum(cell2mat(VDJdata{j}(:, MutIdx)), 2); 
    TreeTrunk(j) = TotCtSHM(1);
    TreeHeight(j) = max(TotCtSHM);
    VDJdata{j} = VDJdata{j}(1, :);
end
VDJdata = joinData(VDJdata, Map, 'stable');
AddHeader = {'TreeTrunk', 'TreeHeight'};
AddData = num2cell([TreeTrunk TreeHeight]);