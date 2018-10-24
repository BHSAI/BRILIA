%getBlankDataTable will create an empty VDJdata for BRILIA depending on
%what type of annotation is being done (e.g., heavy chain, light chain,
%both).
%
%  [VDJdata, VDJheader] = getBlankDataTable(N, Option)
%
%  INPUT
%    N: Number of rows to have in the table (equal to number of sequences)
%    Option ['H', 'L', 'HL']: Specifies which data table to return 
%
%  OUTPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%
function [VDJdata, VDJheader] = getBlankDataTable(N, Option)
if nargin < 2
    Option = 'HL';
end

%Obtain the VDJheader info for output format
HeaderData = readDlmFile('DataHeaderInfo.csv', 'Delimiter', ','); 
VDJheader = HeaderData(2:end, 1)';

%Determine which columns to use for the IG chain
CommonIdx = find(startsWith(VDJheader, '#COMMON', 'ignorecase', true), 1);
HeavyIdx  = find(startsWith(VDJheader, '#HEAVY', 'ignorecase', true), 1);
LightIdx  = find(startsWith(VDJheader, '#LIGHT', 'ignorecase', true), 1);
MiscIdx   = find(startsWith(VDJheader, '#MISC', 'ignorecase', true), 1);

%Generate the blank data table
switch upper(Option)
    case 'H'
        DataCol = [CommonIdx+1:HeavyIdx-1  HeavyIdx+1:LightIdx-1];
    case 'L'
        DataCol = [CommonIdx+1:HeavyIdx-1  LightIdx+1:MiscIdx-1];
    otherwise
        DataCol = [CommonIdx+1:HeavyIdx-1  HeavyIdx+1:LightIdx-1  LightIdx+1:MiscIdx-1];
end
VDJheader = VDJheader(DataCol);
VDJdata = cell(N, length(VDJheader));
