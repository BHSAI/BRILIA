%parseFileName will take a filename string and parse it into Path,
%Name, Ext, and Pre. It is a wrapper for Matlab's fileparts.m
%function, with added check for file or dir feature. 
%
%  [Path, Name, Ext, Pre, Exist] = parseFileName(FullFileName)
%
%  INPUT
%    FullFileName: file path and name to be parsed
%
%  OUTPUT
%    Path: File path for the FullFileName.        EX: 'C:\Users\PathFolder'
%    Name: File name, including extension.        Ex: 'test.csv'
%    Ext:  File extension, including the dot.     Ex: '.csv'
%    Pre:  File name without path and extension.  Ex: 'test'
%    Exist: 1 if the file or dir exist, or 0.     Ex: 0
%
%  EXAMPLE
%    [Path, Name, Ext, Pre, Exist] = parseFileName('C:\UsrFolder\test.csv')
%    Path =  
%       'C:\UsrFolder'
%    Name =  
%       'test.csv'
%    Ext  =  
%       '.csv'
%    Pre  =
%       'test'
%    Exist = 
%       0
%
%    [Path, Name, Ext, Pre, Exist] = parseFileName('c:\UsrFolder\')
%    Path =
%      'c:\UsrFolder'
%    Name =
%      0×0 empty char array
%    Ext =
%      0×0 empty char array
%    Pre =
%      0×0 empty char array
%    Exist =
%      0
%
function [FilePath, FileName, FileExt, FilePre, Exist] = parseFileName(FullFileName, varargin)
if nargin > 1
    error('%s: A 2nd input is no longer supported.', mfilename);
end
if isempty(FullFileName); FullFileName = ''; end
[FilePath, FilePre, FileExt] = fileparts(strtrim(FullFileName));
if isempty(FileExt) %It was a folder, so rejoin
    FilePath = fullfile(FilePath, FilePre);
    FileName = '';
    FilePre = '';
    SearchFor = 'dir';
else
    FileName = [FilePre FileExt];
    SearchFor = 'file';
end
Exist = exist(FullFileName, SearchFor) > 0;