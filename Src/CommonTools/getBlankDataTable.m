%getBlankDataTable will create an empty VDJdata for BRILIA depending on
%what type of annotation is being done (e.g., heavy chain, light chain,
%both.
%
%  [VDJdata, VDJheader] = getBlankDataTable(N, Option)
%
%  INPUT
%    N: Number of rows to have in the table (equal to number of sequences)
%    Option ['H', 'L', 'HL']: Specifies which data table to return. 
%
%  OUTPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%
%  See also convertInput2VDJdata

function [VDJdata, VDJheader] = getBlankDataTable(N, Option)
%Obtain the VDJheader info for output format
HeaderData = readDlmFile('DataHeaderInfo.csv', 'Delimiter', ','); 
VDJheader = HeaderData(2:end, 1)';

%Determine which columns to use for the IG chain
CommonLoc = findCell(VDJheader, '#COMMON', 'MatchWord', 'Any');
HeavyLoc  = findCell(VDJheader, '#HEAVY', 'MatchWord', 'Any');
LightLoc  = findCell(VDJheader, '#LIGHT', 'MatchWord', 'Any');
MiscLoc   = findCell(VDJheader, '#MISC', 'MatchWord', 'Any');
Option = upper(Option);
switch Option
    case 'H'
        DataCol = [CommonLoc+1:HeavyLoc-1 HeavyLoc+1:LightLoc-1];
    case 'L'
        DataCol = [CommonLoc+1:HeavyLoc-1 LightLoc+1:MiscLoc-1];
    otherwise
        DataCol = [CommonLoc+1:HeavyLoc-1 HeavyLoc+1:LightLoc-1 LightLoc+1:MiscLoc-1];
end

%Select only relevant columns
VDJheader = VDJheader(DataCol);
VDJdata = cell(N, length(VDJheader));
