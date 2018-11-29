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
if isempty(regexp(TempDir, filesep, 'once'))
    error('%s: TempDir "%s" should have "%s" in it.', mfilename, TempDir, filesep);
end
if TempDir(end) ~= filesep
    TempDir = [TempDir filesep];
end

%Deletes empty Temp dir
if ~isempty(varargin) && ischar(varargin{1}) && strcmpi(varargin{1}, 'delete')
    if ~exist(TempDir, 'dir')
        return
    end
    FileContents = dir([TempDir '*.*']);
    if length(FileContents) <= 2
        [Success, Msg] = rmdir(TempDir);
        assert(Success == 1, '%s: Error removing temp dir "%s".\n  %s\n', mfilename, TempDir, Msg);
    end
    return
end

%Make a new Temp dir
if ~exist(TempDir, 'dir')
    [Success, Msg] = mkdir(TempDir);
    assert(Success == 1, '%s: Error making temp dir "%s".\n  %s\n', mfilename, TempDir, Msg);
    return
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
