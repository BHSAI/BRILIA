%parseFileName will take a filename string and parse it into Path, Name,
%and FileExt. This is better than MATLAB's fileparts, that yields filepaths
%that lack the ending slash, or yields filenames without the extensions.
%
%  [FilePath, FileName, FileExt] = parseFileName(FullName)
%
%  EX: 
%    FullName = 'C:\Users\Desktop\testing.xlsx';
%    [FilePath, FileName, FileExt] = parseFileName(FullName)
%      FilePath =  C:\Users\Desktop\
%      FileName =  testing.xlsx
%      FileExt =  .xlsx
function [FilePath, FileName, FileExt] = parseFileName(FullName)
SlashLoc = regexp(FullName,'\\|\/'); %Location of fwd or bwd slashes
DotLoc =regexp(FullName,'\.'); %Location of dots

if isempty(SlashLoc)
    if ispc
        FullName = [cd '\' FullName];
    else
        FullName = [cd '/' FullName];
    end
    SlashLoc = regexp(FullName,'\\|\/'); %Location of fwd or bwd slashes
    DotLoc =regexp(FullName,'\.'); %Location of dots
end

if isempty(DotLoc)
    FilePath = FullName;
    FileName = [];
    FileExt = [];
else
    %Determine file extension
    if isempty(DotLoc) || (~isempty(DotLoc) && DotLoc(end) < SlashLoc(end))
        FileExt = [];
    else
        FileExt = FullName(DotLoc(end):end);
    end

    %Determine the file name
    if isempty(SlashLoc)
        FileName = FullName;
    else
        FileName = FullName(SlashLoc(end)+1:end);
    end

    %Determine the file path
    if isempty(DotLoc)
        FilePath = [];
    else
        FilePath = FullName(1:SlashLoc(end));
    end
end
