%renameVDJheader will each sequence data file VDJheader from one string to
%another. This is used mainly for "cleaning up" data header as new version
%comes out and a more standard header name is adopted.
%
%  renameVDJheader(FileNames, OrigH, NewH)
%
%  renameVDJheader([], OrigH, NewH)
%
%  renameVDJheader('', OrigH, NewH)
%
%  INPUT
%    FileNames: file name string or cell array of strings
%    OrigH: original header string or cell array of strings
%    NewH: new header string or cell array of strings
%
%  NOTE
%    This is case sensitive. Only EXACTLY matching strings are replaced
%    with EXACTLY what is specified.
%
%    If no match are found, no changes will be made.
function renameVDJheader(FileNames, OrigH, NewH)
if ischar(OrigH)
    OrigH = {OrigH};
end
if ischar(NewH)
    NewH = {NewH};
end
if numel(OrigH) ~= numel(NewH)
    error('%s: The number of headers to change must be same.', mfilename);
end

if isempty(FileNames) 
    FileNames = {};
elseif ischar(FileNames) 
    if exist(FileNames, 'file') 
        FileNames = {FileNames};
    else
        FileNames = {};
    end
elseif iscell(FileNames)
    FileExist = cellfun(@(x) exist(x, 'file')>0, FileNames);
    FileNames = FileNames(FileExist);
else
    FileNames = {};
end
if isempty(FileNames)
    FileNames = getBriliaFiles;
end

for f = 1:length(FileNames)
    [VDJheader, Delimiter] = readDlmFile(FileNames{f}, 'LineRange', [1 1]);
    [~, OldIdx, NewIdx] = intersect(VDJheader, OrigH);
    if ~isempty(OldIdx)
        VDJheader(OldIdx) = NewH(NewIdx);
        VDJdata = readDlmFile(FileNames{f});
        VDJdata(1, :) = VDJheader;
        writeDlmFile(VDJdata, FileNames{f}, Delimiter)
    end
end