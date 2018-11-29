%getLightHeaderVar is used to return the column number associated with the
%data type of VDJdata via a short naming system and ending with "Loc". It
%is NOT an enumeration class because newer version of VDJdata can have
%different headers, hence once must alway compute the new column locations.
%
%  L = getLightHeaderVar(VDJheader)
%
%  [L, N] = getLightHeaderVar(VDJheader)
%
%  [L, N, C] = getLightHeaderVar(VDJheader)

%  INPUT
%    VDJheader: main BRILIA header cell
%
%  OUTPUT
%    L: standardized header name variable
%      L.SeqNameLoc
%      L.SeqNumLoc
%      L.GrpNumLoc
%      L.TemplateLoc
%      L.ChildCountLoc
%      L.SeqLoc
%      L.RefSeqLoc
%      L.FunctLoc
%      L.CDR1Loc
%      L.CDR2Loc
%      L.CDR3Loc
%      L.VJmapLoc
%      L.VJnameLoc
%      L.LengthLoc
%      L.DelLoc
%      L.VmutLoc
%      L.NmutLoc
%      L.JmutLoc
%    N: Matrix of data column indices holding numerical values
%    C: Matrix of data column indices holding character values
%
%  See also getHeavyHeaderVar, getAllHeaderVar
function [L, varargout] = getLightHeaderVar(VDJheader)
warning('%s: This is deprecated. Use getVDJmapper.', mfilename);

%Extract common fields
L.SeqNameLoc    = findCell(VDJheader, 'SeqName', 'MatchCase', 'Any');
L.SeqNumLoc     = findCell(VDJheader, 'SeqNum', 'MatchCase', 'Any');
L.GrpNumLoc     = findCell(VDJheader, 'GroupNum', 'MatchCase', 'Any');
L.TemplateLoc   = findCell(VDJheader, {'TemplateCount', 'count (templates)'}, 'MatchCase', 'Any');
L.ChildCountLoc = findCell(VDJheader, 'TreeChildCount', 'MatchCase', 'Any');

%Extract light chain fields
L.SeqLoc        = findCell(VDJheader, 'L-Seq', 'MatchCase', 'Any');
L.RefSeqLoc     = findCell(VDJheader, 'L-RefSeq', 'MatchCase', 'Any');
L.OverSeq5Loc   = findCell(VDJheader, 'L-OverhangSeq5', 'MatchCase', 'Any');
L.OverSeq3Loc   = findCell(VDJheader, 'L-OverhangSeq3', 'MatchCase', 'Any');
L.FunctLoc      = findCell(VDJheader, 'L-Functional', 'MatchCase', 'Any');

L.CDR1Loc       = findCell(VDJheader, {'L-CDR1_AminoAcid', 'L-CDR1_Length', 'L-CDR1_Start', 'L-CDR1_End'}, 'MatchCase', 'Any');
L.CDR2Loc       = findCell(VDJheader, {'L-CDR2_AminoAcid', 'L-CDR2_Length', 'L-CDR2_Start', 'L-CDR2_End'}, 'MatchCase', 'Any');
L.CDR3Loc       = findCell(VDJheader, {'L-CDR3_AminoAcid', 'L-CDR3_Length', 'L-CDR3_Start', 'L-CDR3_End'}, 'MatchCase', 'Any');

L.GeneNumLoc    = findCell(VDJheader, {'L-V_MapNum', 'L-J_MapNum'}, 'MatchCase', 'Any');
L.GeneNameLoc   = findCell(VDJheader, {'L-V_GeneName', 'L-J_GeneName'}, 'MatchCase', 'Any');

L.LengthLoc     = findCell(VDJheader, {'L-Length_V', 'L-Length_N', 'L-Length_J'}, 'MatchCase', 'Any');
L.DelLoc        = findCell(VDJheader, {'L-V_Deletion3', 'L-J_Deletion5'}, 'MatchCase', 'Any');
L.VmutLoc       = findCell(VDJheader, 'L-SHM_V', 'MatchCase', 'Any');
L.NmutLoc       = findCell(VDJheader, 'L-SHM_N', 'MatchCase', 'Any');
L.JmutLoc       = findCell(VDJheader, 'L-SHM_J', 'MatchCase', 'Any');

L.ValignLoc     = findCell(VDJheader, 'L-V_Align');
L.JalignLoc     = findCell(VDJheader, 'L-J_Align');

if nargout >= 2
    %Specify the data type
    N.SeqNameLoc    = 0;
    N.SeqNumLoc     = 1;
    N.GrpNumLoc     = 1;
    N.TemplateLoc   = 1;
    N.ChildCountLoc = 1;
    N.SeqLoc        = 0;
    N.RefSeqLoc     = 0;
    N.OverSeq5Loc   = 0;
    N.OverSeq3Loc   = 0;
    N.FunctLoc      = 0;
    N.CDR1Loc       = [0 1 1 1];
    N.CDR2Loc       = [0 1 1 1];
    N.CDR3Loc       = [0 1 1 1];
    N.GeneNumLoc    = [1 1];
    N.GeneNameLoc   = [0 0];
    N.LengthLoc     = [1 1 1];
    N.DelLoc        = [1 1];
    N.VmutLoc       = 1;
    N.NmutLoc       = 1;
    N.JmutLoc       = 1;
    N.ValignLoc     = 0;
    N.JalignLoc     = 0;
       
    FieldNames = fieldnames(N);
    LocType = zeros(50, 2); %Preallocate 50 indices, for now. Col 1 is data index, Col 2 is the data type (0 = char, 1 = num);
    k = 1; %LocType counter
    for j = 1:length(FieldNames)
        FieldValues = L.(FieldNames{j});
        FieldTypes  = N.(FieldNames{j});
        LocType(k:k+length(FieldValues)-1, 1) = FieldValues;
        LocType(k:k+length(FieldTypes)-1, 2)  = FieldTypes;
        k = k + length(FieldValues);
    end
    LocType(k:end, :) = [];
    LocType(LocType(:,1) == 0, :) = [];
    varargout{1} = LocType(LocType(:, 2) == 1, 1); %Numerical indices
    if nargout >= 3
        varargout{2} = LocType(LocType(:, 2) == 0, 1); %Character indices
    end
end
