%prepTempDir will prepare a Temp dir for a particular BRILIA run to
%ensure that the TempDir exists and that no Raw.csv and Final.csv exists.
%This is the preferred method of preparing the Temp dir, just in case the
%user stores files into the Temp folder or has a same name Temp folder, in
%which BRILIA should NOT delete all contents inside the Temp dir and cause
%headaches.
%
%  prepTempDir(TempDir)
%
%  prepTempDir(TempDir, 'delete')
%
%  INPUT
%    TempDir: the path where files will be temporarily stored when BRILIA
%      is running.
%    'delete': will delete the TempDir ONLY if it is empty

function prepTempDir(TempDir, varargin)
SlashType = TempDir(regexp(TempDir, '\\|\/', 'once'));
if isempty(SlashType)
    error('%s: TempDir should have a slash in its name', mfilename);
end
if TempDir(end) ~= SlashType
    TempDir = cat(2, TempDir, SlashType);
end

%Deletes empty Temp dir
if ~isempty(varargin) && ischar(varargin{1}) && strcmpi(varargin{1}, 'delete')
    if ~exist(TempDir, 'dir'); 
        return; 
    end
    FileContents = dir([TempDir '*.*']);
    if length(FileContents) <= 2
        try
            rmdir(TempDir)
        catch
        end
    end
    return;
end

%Make a new Temp dir
if ~exist(TempDir, 'dir')
    try
        mkdir(TempDir);
    catch
        error('%s: Could not make temp dir %s . Check permission.\n', mfilename, TempDir);
    end
    return;
end

%Remove Raw.csv files
RawFileStruct = dir([TempDir '*Raw.csv']);
for j = 1:length(RawFileStruct)
    delete([TempDir RawFileStruct(j).name]);
end

%Remove Final.csv files
RawFileStruct = dir([TempDir '*Final.csv']);
for j = 1:length(RawFileStruct)
    delete([TempDir RawFileStruct(j).name]);
end
