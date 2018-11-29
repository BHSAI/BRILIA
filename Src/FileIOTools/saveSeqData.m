%saveSeqData will save VDJdata and VDJheader to an output delimited file.
%
%  saveSeqData(FullSaveName, VDJdata, VDJheader)
%
%  saveSeqData(FullSaveName, VDJdata, VDJheader, 'Delimiter', Delimiter)
%
%  saveSeqData(..., 'append')
%
%  INPUT
%    VDJdata: BRILIA data table
%    VDJheader: BRILIA header table for VDJdata
%    Delimiter: can be ', ' or ';' or '\t'
%    'append': will append VDJdata after existing file, or create new one.
%      Note that append can be specified anywhere.
%
%  NOTE
%    If the delimiter character exists in a string field, it will replace
%    the delimiter character with '_' to prevent confusion in opening
%    delimited files later.

function saveSeqData(FullSaveName, VDJdata, VDJheader, varargin)
%Parse the input for append mode and delimiter
AppendThis = 'n';
Delimiter = ',';
for j = 1:length(varargin)
    if ischar(varargin{j}) 
        if ismember(lower(varargin{j}), {'-append', 'append'})
            AppendThis = 'y';
        elseif strcmpi(varargin{j}, 'Delimiter') && length(varargin) > j+1
            if ismember(lower(varargin{j+1}), {';' ',' '\t'})
                Delimiter = varargin{j+1};
            end
        end
    end
end

%Write the data
if AppendThis == 'y'
    if exist(FullSaveName, 'file') == 0 %First entry, so include header
        writeDlmFile([VDJheader; VDJdata], FullSaveName, Delimiter, 'append');    
    else %Do not include header when appending.
        writeDlmFile(VDJdata, FullSaveName, Delimiter, 'append');
    end
else 
    writeDlmFile([VDJheader; VDJdata], FullSaveName, Delimiter);
end
