%openMultSeqData will open multiple BRILIA output files into a single
%structure. This is used for comparing repertoire annotations.
%
%  S = openMultSeqData
%
%  S = openMultSeqData(FileNames)
%
%  S = openMultSeqData(..., Header1, Header2, ...)
%
%  INPUT
%    FileNames: char or cell of char of file names to load. Use
%      getBriliaFiles to get a cell of char file names.
%    HeaderN: char of the VDJheader name that you want to extract. Ex:
%      h-seq, h-Vgene, etc. (The "-" and case are ignored in the search)
%
%  OUTPUT
%    S: non-scalar structure of BRILIA output files with 3 fields.
%      S(N).VDJdata: store the VDJdata of the Nth file
%      S(N).VDJheader: store the header names of the Nth file
%      S(N).FileName: full file name of the data file of the Nth file
%
%  See also getBriliaFiles


function S = openMultSeqData(varargin)
if nargin > 0
    if isempty(varargin{1}) 
        FileNames = {};
        varargin(1) = [];
    elseif ischar(varargin{1}) 
        if exist(varargin{1}, 'file') 
            FileNames = {varargin{1}};
            varargin(1) = [];
        else
            FileNames = {};
        end
    elseif iscell(varargin{1})
        FileExist = cellfun(@(x) exist(x, 'file')>0, varargin{1});
        FileNames = varargin{1}(FileExist);
    else
        FileNames = {};
    end
    if isempty(FileNames)
        FileNames = getBriliaFiles;
    end
else
    FileNames = getBriliaFiles;
end

if isempty(FileNames)
    warning('%s: No files were selected or are valid', mfilename);
end

S(1:length(FileNames)) = struct('VDJdata', [], 'VDJheader', [], 'FileName', '');
DelLoc = zeros(length(FileNames), 1, 'logical');
for f = 1:length(FileNames)
    [S(f).VDJdata, S(f).VDJheader, FileName, FilePath] = openSeqData(FileNames{f}, varargin{:});
    if isempty(FileName)
        DelLoc(f) = 1;
        continue
    end
    S(f).FileName = fullfile(FilePath, FileName);
end
S(DelLoc) = [];