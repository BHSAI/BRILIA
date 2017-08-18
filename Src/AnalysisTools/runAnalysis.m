%runAnalysis will perform a series of common analysis to an existing BRILIA
%output data, returning a folder of plots and tables used in analysis.
%
%  INPUT
%    FileName: file name of the BRILIA output file
%    SaveDir: directory to save the file to. If not specified or empty,
%      will save to the same folder as FileName.
%    SaveSubDir: a subdirectory that will ALWAYS be made and used to store
%      the output file. If empty, will use SaveDir.

function varargout = runAnalysis(varargin)
P = inputParser;
addOptional(P,  'InputFile', '', @(x) ischar(x) || iscell(x) || isempty(x));
addParameter(P, 'SaveDir', '', @(x) ischar(x) || isempty(x));
addParameter(P, 'SaveSubDir', 'Analysis', @(x) ischar(x) || isempty(x));
[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return;
end
SaveSubDir = Ps.SaveSubDir;
SaveDir = Ps.SaveDir;
InputFile = Ps.InputFile;

[VDJdata, VDJheader, ~, FilePath] = openSeqData(InputFile);
if isempty(VDJdata)
    return;
end
if isempty(SaveDir) && ~isempty(FilePath)
    SaveDir = FilePath;
end

plotTree(VDJdata, VDJheader, 'SaveDir', SaveDir, 'SaveSubDir', 'Tree', 'SaveAs', 'Tree.png');

MotifData = collectMotifData(VDJdata, VDJheader);
HotMotifData = getHotMotifData(MotifData);

SearchTable = cell(5, 3);
SearchTable(1, :) = {'Title', 'Description', 'FileName'};
j = 2;

try
    SearchTable{j, 3} = printMotifData(MotifData, 'Normalize', 'n', 'SaveDir', SaveDir, 'SaveSubDir', SaveSubDir, 'SaveAs', 'Mutability.csv');
    SearchTable{j, 2} = 'Table of 3-nt motifs'' mutation frequencies from nt X to nt Y.';
    SearchTable{j, 1} = 'Trinucleotide Mutability Data';
    j = j+1;
catch
    error('%s: Could not run printMotifData, not normalized.', mfilename);
end

try
    SearchTable{j, 3} = printMotifData(MotifData, 'Normalize', 'y', 'SaveDir', SaveDir, 'SaveSubDir', SaveSubDir, 'SaveAs', 'NormMutability.csv');
    SearchTable{j, 2} = 'Table of 3-nt motifs'' normalized mutation frequencies from nt X to nt Y.';
    SearchTable{j, 1} = 'Normalized Trinucleotide Mutability Data';
    j = j+1;
catch
    error('%s: Could not run printMotifData, normalize.', mfilename);
end

try
    SearchTable{j, 3} = plotHotMotifBarGraph(HotMotifData, 'SaveDir', SaveDir, 'SaveSubDir', SaveSubDir, 'SaveAs', 'HotMotifBarGraph.png', 'Visible', 'off');
    SearchTable{j, 2} = 'Bar charts of hotspot motifs'' normalized mutation frequencies from nt X to nt Y.';
    SearchTable{j, 1} = 'Hotspot Mutation Frequencies Plot';
    j = j+1;
catch
    error('%s: Could not run plotHotMotifBarGraph.', mfilename);
end

try
    SearchTable{j, 3} = plotHotMotifDendrogram(HotMotifData, 'SaveDir', SaveDir, 'SaveSubDir', SaveSubDir, 'SaveAs', 'HotMotifBarDendrogram.png', 'Visible', 'off');
    SearchTable{j, 2} = 'Dendrogram of similarity in hotspot motifs'' mutation frequencies';
    SearchTable{j, 1} = 'Hotspot Mutation Similarity Plot';
    j = j+1;
catch
    error('%s: Could not run plotHotMotifDendrogram.', mfilename);
end

%--------------------------------------------------------------------------
%Preparing search table for saving

for k = 1:size(SearchTable, 1) 
    [~, SearchTable{k, 3}, ~] = parseFileName(SearchTable{k, 3}, 'ignorefilecheck');
end

SearchTableFileName = prepSaveTarget('SaveDir', SaveDir, 'SaveSubDir', SaveSubDir, 'SaveAs', 'AnalysisSearchTable', 'SaveExt', '.csv');
if exist(SearchTableFileName, 'file') %Open file, replace redundant, and than reseave. EDIT.
    PastSearchTable = readDlmFile(SearchTableFileName);
    SameLoc = findCell(PastSearchTable(2:end, 1), SearchTable(2:end));
    if SameLoc(1) > 0
        PastSearchTable(SameLoc + 1, :) = [];
    end
    SearchTable = cat(1, PastSearchTable, SearchTable(2:end, :));
end
writeDlmFile(SearchTable, SearchTableFileName);

varargout{1} = SearchTableFileName;
