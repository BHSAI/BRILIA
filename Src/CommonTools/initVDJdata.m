%initVDJdata will initialize a new VDJ data structure with empty fields.
%
%  S = initVDJdata(N)
%
%  S = initVDJdata(N, Chain)
%
%  INPUT
%    N: number of entries
%    Chain ['H', 'L', or 'HL']: heavy and/or light chain data fields to use
%                               Defaults to 'H' if empty.
%
%  OUTPUT
%    S: a nonscalar structure with all the fields needed for BRILIA

function S = initVDJdata(N, varargin)
if N < 1
    error('%s: Input must be an integer > 0', mfilename);
else
    N = round(N);
end
if nargin == 1
    Chain = 'H';
elseif ischar(varargin{1}) && ismember(upper(varargin{1}), {'H', 'L', 'HL'})
    Chain = upper(varargin{1});
else
    error('%s: Chain must be H, L or HL.', mfilename);
end

%Common fields
Fields{1} = {
    'SeqName'   , '';
    'SeqNum'    , '';
    'GrpNum'    , '';
    'Chain'     , '';
    'Copies'    , '';
    'Children'  , ''};

%Heavy chain fields
Fields{2} = {
    'hSeq'      , '';
    'hRefSeq'   , '';
    'hTrimSeq5' , '';
    'hTrimSeq3' , '';
    'hFunction' , '';
    'hCDR1'     , '';
    'hCDR2'     , '';
    'hCDR3'     , '';
    'hVgene'    , '';
    'hDgene'    , '';
    'hJgene'    , '';
    'hValign'   , '';
    'hDalign'   , '';
    'hJalign'   , '';

    'hCDR1_Bgn' , '';
    'hCDR1_End' , '';
    'hCDR1_Len' , '';		
    'hCDR2_Bgn' , '';
    'hCDR2_End' , '';
    'hCDR2_Len' , '';
    'hCDR3_Bgn' , '';
    'hCDR3_End' , '';
    'hCDR3_Len' , '';
    
    'hLen_V'    , '';
    'hLen_M'    , '';
    'hLen_D'    , '';
    'hLen_N'    , '';
    'hLen_J'    , '';
    'hSHM_V'    , '';
    'hSHM_M'    , '';
    'hSHM_D'    , '';
    'hSHM_N'    , '';
    'hSHM_J'    , '';
    'hDel_V3'   , '';
    'hDel_D5'   , '';
    'hDel_D3'   , '';
    'hDel_J5'   , '';
    'hMap_V'    , '';
    'hMap_D'    , '';
    'hMap_J'    , '';
    'hScore_V'  , '';
    'hScore_D'  , '';
    'hScore_J'  , '';

    'hCgene'    , '';
    'hC_Bgn'    , '';
    'hC_End'    , ''};			

%Light chain fields
Fields{3} = {
    'lSeq'      , ''; 
    'lRefSeq'   , '';
    'lTrimSeq5' , '';
    'lTrimSeq3' , '';
    'lFunction' , '';
    'lCDR1'     , '';
    'lCDR2'     , '';
    'lCDR3'     , '';
    'lVgene'    , '';
    'lJgene'    , '';
    'lValign'   , '';
    'lJalign'   , '';

    'lCDR1_Bgn' , '';
    'lCDR1_End' , '';
    'lCDR1_Len' , '';
    'lCDR2_Bgn' , '';
    'lCDR2_End' , '';
    'lCDR2_Len' , '';
    'lCDR3_Bgn' , '';
    'lCDR3_End' , '';
    'lCDR3_Len' , '';
    
    'lLen_V'    , '';
    'lLen_N'    , '';
    'lLen_J'    , '';
    'lSHM_V'    , '';
    'lSHM_N'    , '';
    'lSHM_J'    , '';
    'lDel_V3'   , '';
    'lDel_J5'   , '';
    'lMap_V'    , '';
    'lMap_J'    , '';
    'lScore_V'  , '';
    'lScore_J'  , ''};

switch upper(Chain)
    case 'HL'
        Idx = [1 2 3];
    case 'H'
        Idx = [1 2];
    case 'L'
        Idx = [1 3];
end
Fields = cat(1, Fields{Idx})';
S(1:N, 1) = struct(Fields{:});