%parseFileName will take a filename string and parse it into Path, Name,
%and FileExt. In contrast to MATLAB's fileparts, this will attempt to look
%for the full path if not provided.
%
%  [FilePath, FileName, FileExt, FileNamePre] = parseFileName(FullName)
%
%  [...] = parseFileName(FullName, 'ignorefilecheck')
%
%  INPUT
%    FullName: file path and name to be parsed
%    'ignorefilecheck': use this to prevent parseFileName searching for the
%    full file path ONLY IF the FullName lacks the file path.
%
%  OUTPUT
%    FilePath: File path for the FullName. EX: 'C:\Users\User1\PathFolder\'
%    FileName: File name including any extensions.    Ex: 'test.csv'
%    FileExt: File extension only, including the dot. Ex:     '.csv'
%    FileNamePre: File name excluding any extensions. Ex: 'test'
%
%  EX: 
%    F = 'C:\Users\Desktop\testing.xlsx';
%    [FilePath, FileName, FileExt] = parseFileName(F, 'ignorefilecheck')
%    FilePath =  
%       'C:\Users\Desktop\'
%    FileName =  
%       'testing.xlsx'
%    FileExt  =  
%       '.xlsx'
function [FilePath, FileName, FileExt, FileNamePre] = parseFileName(FullName, varargin)
IgnoreFileCheck = any(strcmpi(varargin, 'ignorefilecheck'));
[FilePath, FileNamePre, FileExt] = fileparts(strtrim(FullName));
FileName = [FileNamePre FileExt];
if isempty(FilePath) && ~IgnoreFileCheck %Want the full file path too
    FilePath = fileparts(which(FileName));
end
if ~isempty(FilePath)
    FilePath = fullfile(FilePath, filesep);
end
% 
% %Establish defaults
% FileName = [];
% FileExt = [];
% FileNamePre = [];
% HasPath = false;
% 
% %Switch the incorrect slashes
% SlashType  = '\/';
% WrongSlash = SlashType(~(SlashType == filesep));
% FullName = strrep(FullName, WrongSlash, filesep);
% 
% %Case when no file path is given, assume cd is file path
% SlashLoc = find(FullName == filesep);
% if isempty(SlashLoc) %No path was given
%     FilePath = [pwd filesep]; %Assume FilePath is cd for now
%     FileName = FullName;
%    
% else %Path was given
%     HasPath = true;
%     FilePath = FullName(1:SlashLoc(end));
%     if SlashLoc(end) < length(FullName)
%         FileName = FullName(SlashLoc(end)+1:end);
%     end
% end
% 
% %Determine ext and filename without ext
% if ~isempty(FileName) 
%     DotLoc = find(FileName == '.');
%     if ~isempty(DotLoc)
%         if length(FileName) - DotLoc(end) <= 7 %.xxxx file extension 
%             FileExt = FileName(DotLoc(end):end);
%             FileNamePre = FileName(1:DotLoc(end)-1);
%         end
%     end
% end
% 
% %If no path was given, file does not exist, and ignorefilecheck option not
% %given, return an empty path to show this file does not exist
% if ~HasPath && ~exist(fullfile(pwd, FileName), 'file') && ~IgnoreFileCheck
%     FilePath = [];
% end