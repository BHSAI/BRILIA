%parseFileName will take a filename string and parse it into Path, Name,
%and FileExt. This is better than MATLAB's fileparts, which yields
%filepaths that lack the ending slash, or yields filenames without the
%extensions.
%
%  [FilePath, FileName, FileExt] = parseFileName(FullName)
%
%  [FilePath, FileName, FileExt] = parseFileName(FullName, 'ignorefilecheck')
%
%  [FilePath, FileName, FileExt, FileNamePre] = parseFileName(FullName, 'ignorefilecheck')
%
%  INPUT
%    FullName: file path and name to be parsed
%    'ignorefilecheck': use this to prevent parseFileName from checking if
%      file exists, ONLY when file path is not provided in the
%      FullName. This is used mainly for parsing a save file name. Note
%      that without this option, FilePath will return [] if file does not
%      exist. If ignorefilecheck, will return current directory as file the
%      path.
%
%  OUTPUT
%    FilePath: File path for the FullName, including end slasshes EX:
%      'C:\Users\User1\PathFolder\'
%    FileName: File name including any extensions, if any 'test.csv'
%    FileExt: File extension only, including the dot. Ex: '.csv'
%    FileNamePre: File name excluding any extensions.
%
%  EX: 
%    F = 'C:\Users\Desktop\testing.xlsx';
%    [FilePath, FileName, FileExt] = parseFileName(F, 'ignorefilecheck')
%      FilePath =  C:\Users\Desktop\
%      FileName =  testing.xlsx
%      FileExt =  .xlsx
function [FilePath, FileName, FileExt, FileNamePre] = parseFileName(FullName, varargin)
%Check for file existence
if any(contains(varargin, 'ignorefilecheck'))
    IgnoreFileCheck = true;
else
    IgnoreFileCheck = false;
end

%Establish defaults
FileName = [];
FileExt = [];
FileNamePre = [];
HasPath = false;

%Switch the incorrect slashes
SlashType  = '\/';
WrongSlash = SlashType(~(SlashType == filesep));
FullName = strrep(FullName, WrongSlash, filesep);

%Case when no file path is given, assume cd is file path
SlashLoc = find(FullName == filesep);
if isempty(SlashLoc) %No path was given
    FilePath = [pwd filesep]; %Assume FilePath is cd for now
    FileName = FullName;
   
else %Path was given
    HasPath = true;
    FilePath = FullName(1:SlashLoc(end));
    if SlashLoc(end) < length(FullName)
        FileName = FullName(SlashLoc(end)+1:end);
    end
end

%Determine ext and filename without ext
if ~isempty(FileName) 
    DotLoc = find(FileName == '.');
    if ~isempty(DotLoc)
        if length(FileName) - DotLoc(end) <= 7 %.xxxx file extension 
            FileExt = FileName(DotLoc(end):end);
            FileNamePre = FileName(1:DotLoc(end)-1);
        end
    end
end

%If no path was given, file doew not exist, and ignorefilecheck option not
%given, return an empty path to show this file does not exist
if ~HasPath && ~exist(fullfile(pwd, FileName), 'file') && ~IgnoreFileCheck
    FilePath = [];
end