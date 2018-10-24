%getHeavyHeaderVar is used to return the column number associated with the
%data type of VDJdata via a short naming system and ending with "Loc". It
%is NOT an enumeration class because newer version of VDJdata can have
%different headers, hence once must alway compute the new column locations.
%
%  H = getHeavyHeaderVar(VDJheader)
%
%  [H, N] = getHeavyHeaderVar(VDJheader)
%
%  [H, N, C] = getHeavyHeaderVar(VDJheader)
%
%  INPUT
%    VDJheader: main BRILIA header cell
%
%  OUTPUT
%    H: standardized header name variable
%      H.SeqNameLoc
%      H.SeqNumLoc
%      H.GrpNumLoc
%      H.TemplateLoc
%      H.ChildCountLoc
%      H.SeqLoc
%      H.RefSeqLoc
%      H.FunctLoc
%      H.CDR1Loc
%      H.CDR2Loc
%      H.CDR3Loc
%      H.VDJmapLoc
%      H.VDJnameLoc
%      H.LengthLoc
%      H.DelLoc
%      H.VmutLoc
%      H.MmutLoc
%      H.DmutLoc
%      H.NmutLoc
%      H.JmutLoc
%    N: Matrix of data column indices holding numerical values
%    C: Matrix of data column indices holding character values
%
%  See also getLightHeaderVar, getAllHeaderVar

function [H, varargout]= getHeavyHeaderVar(VDJheader)
warning('%s: This is deprecated. Use getVDJmapper.', mfilename);

%Extract common fields
H.SeqNameLoc    = findCell(VDJheader, 'SeqName', 'MatchCase', 'Any');
H.SeqNumLoc     = findCell(VDJheader, 'SeqNum', 'MatchCase', 'Any');
H.GrpNumLoc     = findCell(VDJheader, 'GroupNum', 'MatchCase', 'Any');
H.TemplateLoc   = findCell(VDJheader, {'TemplateCount', 'count (templates)'}, 'MatchCase', 'Any');
H.ChildCountLoc = findCell(VDJheader, 'TreeChildCount', 'MatchCase', 'Any');

%Extract heavy chain fields
H.SeqLoc        = findCell(VDJheader, 'H-Seq', 'MatchCase', 'Any');
H.RefSeqLoc     = findCell(VDJheader, 'H-RefSeq', 'MatchCase', 'Any');
H.OverSeq5Loc   = findCell(VDJheader, 'H-OverhangSeq5', 'MatchCase', 'Any');
H.OverSeq3Loc   = findCell(VDJheader, 'H-OverhangSeq3', 'MatchCase', 'Any');
H.FunctLoc      = findCell(VDJheader, 'H-Functional', 'MatchCase', 'Any');

H.CDR1Loc       = findCell(VDJheader, {'H-CDR1_AminoAcid', 'H-CDR1_Length', 'H-CDR1_Start', 'H-CDR1_End'}, 'MatchCase', 'Any');
H.CDR2Loc       = findCell(VDJheader, {'H-CDR2_AminoAcid', 'H-CDR2_Length', 'H-CDR2_Start', 'H-CDR2_End'}, 'MatchCase', 'Any');
H.CDR3Loc       = findCell(VDJheader, {'H-CDR3_AminoAcid', 'H-CDR3_Length', 'H-CDR3_Start', 'H-CDR3_End'}, 'MatchCase', 'Any');

H.GeneNumLoc    = findCell(VDJheader, {'H-V_MapNum', 'H-D_MapNum', 'H-J_MapNum'}, 'MatchCase', 'Any');
H.GeneNameLoc   = findCell(VDJheader, {'H-V_GeneName', 'H-D_GeneName', 'H-J_GeneName'}, 'MatchCase', 'Any');

H.LengthLoc     = findCell(VDJheader, {'H-Length_V', 'H-Length_Nvd', 'H-Length_D', 'H-Length_Ndj', 'H-Length_J'}, 'MatchCase', 'Any');
H.DelLoc        = findCell(VDJheader, {'H-V_Deletion3', 'H-D_Deletion5', 'H-D_Deletion3', 'H-J_Deletion5'}, 'MatchCase', 'Any');
H.VmutLoc       = findCell(VDJheader, 'H-SHM_V', 'MatchCase', 'Any');
H.MmutLoc       = findCell(VDJheader, 'H-SHM_Nvd', 'MatchCase', 'Any');
H.DmutLoc       = findCell(VDJheader, 'H-SHM_D', 'MatchCase', 'Any');
H.NmutLoc       = findCell(VDJheader, 'H-SHM_Ndj', 'MatchCase', 'Any');
H.JmutLoc       = findCell(VDJheader, 'H-SHM_J', 'MatchCase', 'Any');

H.ValignLoc     = findCell(VDJheader, 'H-V_Align');
H.DalignLoc     = findCell(VDJheader, 'H-D_Align');
H.JalignLoc     = findCell(VDJheader, 'H-J_Align');

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
    N.GeneNumLoc    = [1 1 1];
    N.GeneNameLoc   = [0 0 0];
    N.LengthLoc     = [1 1 1 1 1];
    N.DelLoc        = [1 1 1 1];
    N.VmutLoc       = 1;
    N.MmutLoc       = 1;
    N.DmutLoc       = 1;
    N.NmutLoc       = 1;
    N.JmutLoc       = 1;
    N.ValignLoc     = 0;
    N.DalignLoc     = 0;
    N.JalignLoc     = 0;
    
    FieldNames = fieldnames(N);
    LocType = zeros(50, 2); %Preallocate 50 indices, for now. Col 1 is data index, Col 2 is the data type (0 = char, 1 = num);
    k = 1; %LocType counter
    for j = 1:length(FieldNames)
        FieldValues = H.(FieldNames{j});
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
