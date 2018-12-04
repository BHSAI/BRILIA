%testBRILIA is used to test different function of BRILIA.
%
%  test MouseH
%  test MouseL
%  test MouseHL
%  test HumanH
%  test MacaqueH
%  test All

function test(varargin)

ExampleDir = fullfile(findExeDir, 'Examples');
if isdeployed && ~isdir(ExampleDir)
    [Success, Msg] = copyfile(fullfile(findRoot, 'Examples'), ExampleDir, 'f');
    assert(Success, '%s: Error copying file - %s', mfilename, Msg);
end

[FullDirList, DirList] = dir2(ExampleDir, 'dir');
if nargin == 0
    ChoiceNum = chooseFromList(vertcat(DirList, 'All'), 'Attempt', 1, 'Default', [], 'MultiSelect', 'off');
elseif strcmpi(varargin{1}, 'all')
    ChoiceNum = length(DirList) + 1;
else
    ChoiceNum = find(strcmpi(DirList, varargin{1}));
end
if isempty(ChoiceNum)
    return
elseif ChoiceNum <= length(DirList)
    FullDirList = FullDirList(ChoiceNum);
    DirList = DirList(ChoiceNum);
end

for j = 1:length(FullDirList)
    FileName = dir2(fullfile(FullDirList{j}, '*.*'), 'file');
    DirName = DirList{j};
    CapIdx = find(isstrprop(DirName, 'upper'));
    Species = DirName(1:CapIdx(2)-1);
    Chain = DirName(CapIdx(2):end);
    for f = 1:length(FileName)
        [~, ~, FileExt] = parseFileName(FileName{f});
        if startsWith(FileExt, {'.fa', '.csv'}, 'ignorecase', true)
            OutFile = BRILIA(FileName{f}, 'Species', Species, 'Chain', Chain, 'Resume', 'n');
            runAnalysis(OutFile{1});
        end
    end
end