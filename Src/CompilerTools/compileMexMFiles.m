%makeMexMFiles will make the .m files for each MEX file only.
%
%  compileMexMFiles
%
%  compileMexMFiles(FileNames)
%
%  INPUT
%    FileNames: file name string or cell array of strings. If empty, will
%      do all mex files within BRILIA
function compileMexMFiles(FileNames)
if nargin == 0
    RootDir = findRoot;
    Files = dir(fullfile(RootDir, '**', '*.cpp'));
    FileNames = arrayfun(@(x) fullfile(x.folder, x.name), Files, 'unif', false);
else
    if ischar(FileNames)
        FileNames = {FileNames};
    end
    for f = 1:length(FileNames)
        [FilePath, FileName] = parseFileName(FileNames{f});
        FileNames{f} = fullfile(FilePath, FileName);
    end
end

for f= 1:length(FileNames)

    [FilePath, FileName, FileExt] = parseFileName(FileNames{f});
    
    if ~strcmpi(FileExt, '.cpp')
        warning('%s: Skipping file "%s" since it does not end with .cpp', mfilename, FileNames{f});
    end

    MFileName = fullfile(FilePath, [FileName(1:end-length(FileExt)) '.m']);
    recycle('on')
    if exist(MFileName, 'file')
        delete(MFileName)
    end

    FID1 = fopen(FileName, 'r');
    assert(FID1 > 0, '%s: Could not read "%s" ]', mfilename, FileName);
    FID2 = fopen(MFileName, 'w');
    assert(FID2 > 0, '%s: Could not create "%s"', mfilename, MFileName);

    Bgn = 0;
    End = 0;
    j = 0;
    TextLine = fgetl(FID1);
    while ischar(TextLine)
        j = j + 1;
        if Bgn == 0 && contains(TextLine, {'/*', '//'})
            Bgn = j;
        end
        if Bgn > 0 && contains(TextLine, {'*/', '#include'})
            End = j;
            break
        end
        TextLine = fgetl(FID1);
    end

    fseek(FID1, 0, 'bof');
    wroteFirstLine = false;
    for k = Bgn:End
        TextLine = fgetl(FID1);
        Pos1 = regexpi(TextLine, '\S', 'once');
        if isempty(Pos1) && wroteFirstLine
            fprintf(FID2, '%%\r\n');
            continue
        elseif isempty(Pos1) && ~wroteFirstLine
            continue
        end

        Pos2 = Pos1 + 1;
        if Pos2 > length(TextLine)
            continue 
        end

        if contains(TextLine(Pos1:Pos2), {'//', '/*', '*/'})
            if ~wroteFirstLine && isempty(regexp('/S', TextLine(Pos2+1:end), 'once'))
                continue
            end
            fprintf(FID2, '%%%s\r\n', TextLine(Pos2+1:end));
        else
            fprintf(FID2, '%%%s\r\n', TextLine);
        end
        wroteFirstLine = true;
    end
    fclose(FID1);
    fclose(FID2);
end