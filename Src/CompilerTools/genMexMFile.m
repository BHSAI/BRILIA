%genMexMFile will make the .m files for each MEX file so that users can use
%the help function in matlab or showhelp function in BRILIA executable.
%
%  genMexMFile(FileName)
%
%  INPUT
%    FileName: file name or cell array of file names. If empty, will
%      do all *MEX.cpp files within BRILIA src folders
%
%  OUTPUT
%    Text .txt files with the text of each mex file source code above the
%    first #include statement.
function genMexMFile(FileName)
if nargin == 0
    RootDir = findRoot;
    Files = dir(fullfile(RootDir, '**', '*MEX.c*'));
    FileName = fullfile({Files.folder}, {Files.name});
else
    if ischar(FileName)
        FileName = strsplit(FileName, pathsep);
        FileName(cellfun('isempty', FileName)) = [];
    end
    for f = 1:length(FileName)
        [FP, FN] = parseFileName(FileName{f}); %Ensure to get full path
        assert(~isempty(FP), '%s: Could not find the file "%s".', mfilename, FileName{f});
        FileName{f} = fullfile(FP, FN);
    end
end

for f = 1:length(FileName)
    CFileName = FileName{f};
    MFileName = [FileName{f}(1:find(FileName{f} == '.', 1, 'last')) 'm'];
    
    recycle('on')
    if exist(MFileName, 'file')
        delete(MFileName)
    end

    FID1 = fopen(CFileName, 'r');
    assert(FID1 > 0, '%s: Could not read "%s" ]', mfilename, CFileName);
    TXT1 = textscan(FID1, '%s', 'delimiter', '\n', 'whitespace', '');
    TXT1 = TXT1{1};
    fclose(FID1);
    
    EndTxtIdx = find(startsWith(strtrim(TXT1), '#include'), 1) - 1; 
    TXT1 = strrep(TXT1, '/*', '');
    TXT1 = strrep(TXT1, '*/', '');
    TXT1 = strrep(TXT1, '//', '');
    FirstTxtIdx = find(~cellfun('isempty', regexpi(TXT1, '\w+', 'once')), 1);
    
    FID2 = fopen(MFileName, 'w');
    assert(FID2 > 0, '%s: Could not write to "%s"', mfilename, MFileName);
    fprintf(FID2, '%%%s\r\n', TXT1{FirstTxtIdx:EndTxtIdx});
    fclose(FID2);
end
