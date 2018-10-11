%makeHelpFile will search through all m files and directories, pull out
%the help text, and save it to a RootDir\HelpText\mfilename.txt. This is
%used with showHelp(mfilename) to show the help file in a command line
%environment.
%
%  makeHelpFile
%
function makeHelpFile
RootDir = findRoot;
Files = dir(fullfile(RootDir, '**', '*.m'));
FileName = fullfile({Files.folder}, {Files.name});

TargetDir = fullfile(RootDir, 'HelpText');
[Success, Msg] = mkdir(TargetDir);
assert(Success > 0, '%s: Could not create HelpText folder at "%s".\n  %s', mfilename, TargetDir, Msg)

delete(fullfile(TargetDir, '*.txt'));

for f = 1:length(FileName)
    [~, MFileName] = fileparts(FileName{f});
    TFileName = fullfile(TargetDir, [MFileName '.txt']);
    
    FID1 = fopen(FileName{f}, 'r');
    assert(FID1 > 0, '%s: Could not read "%s" ]', mfilename, FileName{f});
    FID2 = []; %Keep empty until you have something to add
    
    while ~feof(FID1)
        Txt = strtrim(fgetl(FID1));
        if startsWith(Txt, 'function')
            break
        elseif startsWith(Txt, '%')
            if isempty(FID2)
                FID2 = fopen(TFileName, 'w');
                assert(FID2 > 0, '%s: Could not write to "%s" ]', mfilename, TFileName);
            end
            fprintf(FID2, '%s\r\n', Txt);
        end
    end
    
    fclose(FID1);
    if ~isempty(FID2)
        fclose(FID2);
    end
end