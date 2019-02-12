function printFreqData(S)

Fields = setdiff(fieldnames(S), {'VDJheader', 'VDJdata', 'Species', 'FileName'});

SizeFiltType = cell(length(Fields), 1);
for k = 1:length(SizeFiltType)
    DashLoc = regexp(Fields{k}, '\_', 'start');
    if isempty(DashLoc); continue; end
    SizeFiltType{k} = Fields{k}(1:DashLoc(1)-1);
end
[UnqSizeFiltType, ~, ~] = unique(SizeFiltType, 'stable');

AllTempData = cell(1, length(S));
for f = 1:length(S)
    TempData = cell(length(Fields), 1);
    for t = 1:length(Fields)
        TempData{t} = S(f).(Fields{t});
        for r = 1:size(TempData{t}, 1)
            if isnumeric(TempData{t}{r, 1})
                TempData{t}{r, 1} = sprintf('%s:%d', Fields{t}, TempData{t}{r, 1});
            else
                TempData{t}{r, 1} = sprintf('%s:%s', Fields{t}, TempData{t}{r, 1});
            end
        end
    end
    TempData = vertcat(TempData{:});
    AllTempData{f} = TempData;
end
AllTempData = conformDist(AllTempData{:});

CellHeader = cell(1, length(S)+1); %empty + FileNames
for f = 1:length(S)
    [FilePath, CellHeader{f+1}] = parseFileName(S(f).FileName);
end
FilePath = fileparts(fileparts(FilePath));

for q = 1:length(UnqSizeFiltType)
    FiltIdx = startsWith(AllTempData(:, 1), UnqSizeFiltType{q});
    SaveName = fullfile(FilePath, [UnqSizeFiltType{q} '_Freq.csv']);
    writeDlmFile([CellHeader; AllTempData(FiltIdx, :)], SaveName, ',');
end