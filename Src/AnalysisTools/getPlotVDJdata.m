%getPlotVDJdata will parse varargin and determine if VDJdata and
%VDJheader are provided, or file name is provided, or nothing is provided.
%If some info leading to VDJdata is provided, will return the varargin
%without that information.
%
%  [VDJdata, VDJheader, FileName, FilePath, varargin] =
%                                              getPlotVDJdata(varargin{:});
%  INPUT
%    varargin: 1xM cell of inputs for other function
% 
%  OUTPUT
%    VDJdata: main BRILIA output data
%    VDJheader: header cell for VDJdata
%    FileName: file name of the brilia csv file, or empty if varargin
%      specified VDJdata, VDJheader 
%    FilePath: file path of the BRILIA csv file, or empty of varargin
%      specified VDJdata, VDJheader
%    varargin: same as varargin, but with any filename or VDJdata,
%      VDJheader variables removed.

function [VDJdata, VDJheader, FileName, FilePath, varargin] = getPlotVDJdata(varargin)
FilePath = '';
FileName = '';

%See if user gave VDJdata
if isempty(varargin) %User specified nothing
    [VDJdata, VDJheader, FileName, FilePath] = openSeqData;
else
    if isempty(varargin{1}) %User specified nothing
        [VDJdata, VDJheader, FileName, FilePath]  = openSeqData;
        varargin(1) = [];
    elseif ischar(varargin{1}) %Maybe user specified file name
        try %User attempted to specify file name
            [VDJdata, VDJheader, FileName, FilePath] = openSeqData(varargin{1});
            varargin(1) = [];
        catch %Specified file name did not work or is not a file name.
            [VDJdata, VDJheader, FileName, FilePath] = openSeqData; %Don't delete from varargin.
        end
    elseif nargin >= 2 && iscell(varargin{1}) && (isstruct(varargin{2}) || (iscell(varargin{2}) && size(varargin{1}, 2) == size(varargin{2}, 2))) %User specified VDJdata nd VDJheader
        VDJdata = varargin{1};
        VDJheader = varargin{2};
        varargin(1:2) = [];
    else %If all fails, just ask user to get VDJdata
        [VDJdata, VDJheader, FileName, FilePath] = openSeqData;
    end        
end
