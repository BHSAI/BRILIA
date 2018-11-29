%printMotifData will save the MotifData as a csv table.
%
%  printMotifData(MotifData, prepSaveTargetVar{:}, writeDLMVar{:})
%
%  INPUT
%    prepSaveTargetVar{:} are param-value pairs of prepSaveTarget
%    writeDlmFileVar{:} are param-value pairs of writeDlmFile
%
%  OUTPUT
%    FullSaveName: the save name of the table

function varargout = printMotifData(MotifData, varargin)
P = inputParser;
addParameter(P, 'Normalize', 'y', @(x) ischar(x) && ismember(lower(x), {'y', 'n'}));
addParameter(P, 'SaveSubDir', 'Analysis', @ischar);
addParameter(P, 'DefaultSaveAs', 'MotifData.csv', @ishcar);
addParameter(P, 'Delimiter', ',', @(x) ischar(x) && ismember(lower(x), {',', ';', '\t'}));
[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return;
end
Normalize = Ps.Normalize;
Delimiter = Ps.Delimiter;
SaveSubDir = Ps.SaveSubDir;
DefaultSaveAs = Ps.DefaultSaveAs;

if strcmpi(Normalize, 'y')
    MotifData = normalizeMotifData(MotifData);
end

%Prepare the empty table first
FullTable = cell(length(MotifData.RowHeader)+1, length(MotifData.ColHeader)*4+1);
FullTable(2:end, 1) = MotifData.RowHeader;
ColHeader = repmat(MotifData.ColHeader, 4, 1);
for j = 1:4
    NT = int2nt(j);
    for k = 1:3
        ColHeader{j, k} = [ColHeader{j, k} '-' NT];
    end
end
ColHeader = ColHeader(:)';
FullTable(1, 2:end) = ColHeader;

%Fill in the table
for r = 1:64
    for c = 1:3
        FullTable(r+1, 4*(c-1)+2:4*c+1) = num2cell(MotifData.Data{r, c});
    end
end

FullSaveName = prepSaveTarget(ExpPu{:}, 'SaveSubDir', SaveSubDir, 'DefaultSaveAs', DefaultSaveAs, 'MakeSaveDir', 'y');
writeDlmFile(FullTable, FullSaveName, Delimiter);
varargout{1} = FullSaveName;
