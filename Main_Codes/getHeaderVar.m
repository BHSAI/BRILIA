%H = getHeaderVar(VDJheader); is used to return the column number associated with the data
%type of VDJdata via a short naming system and ending with "Loc". It is NOT
%an enumeration class because newer version of VDJdata can have different
%headers, hence once must alway compute the new column locations.
%
%  H = H = getHeaderVar(VDJheader); 
%
%  H = H = getHeaderVar(VDJheader);(VDJheader)
%
%  INPUT
%    VDJheader: 1xM cell headers for the VDJdata cell matrix.
%
%  OUTPUT
%    H: standardized header name variable
%      H.SeqLoc
%      H.SeqNumLoc
%      H.SeqNameLoc
%      H.TemplateLoc
%      H.CDR3Loc
%      H.GrpNumLoc
%      H.RefSeqLoc
%      H.FamLoc
%      H.LengthLoc
%      H.DelLoc
%      H.MatchLoc
%      H.SHMLoc
%      H.AlignLoc
%      H.FormClassLoc
%      H.ChildCountLoc
%      H.FunctLoc
%      H.VmutLoc
%      H.MmutLoc
%      H.DmutLoc
%      H.JmutLoc
%      H.LeftTrimLoc
%      H.RightTrimLoc
%      H.MiscLoc
function H = getHeaderVar(VDJheader)
%Extract relevant column locations
H.SeqLoc      = findCell(VDJheader,{'nucleotide','Seq'});
H.SeqNumLoc   = findCell(VDJheader,'SeqNum');
H.SeqNameLoc  = findCell(VDJheader,'SeqName');

H.TemplateLoc = findCell(VDJheader,{'count (templates)','TemplateCount'});
H.CDR3Loc     = findCell(VDJheader,{'CDR3_AminoAcid','CDR3_Length','CDR3_Start','CDR3_End'});

H.GrpNumLoc   = findCell(VDJheader,'GroupNum');
H.RefSeqLoc   = findCell(VDJheader,'RefSeq');
H.FamNumLoc   = findCell(VDJheader,{'vMapNum','dMapNum','jMapNum','V_MapNum','D_MapNum','J_MapNum'});
H.FamLoc      = findCell(VDJheader,{'vMaxResolved','dMaxResolved','jMaxResolved','V_GeneName','D_GeneName','J_GeneName'});

H.LengthLoc   = findCell(VDJheader,{'vAlignLength' 'n2AlignLength' 'dAlignLength' 'n1AlignLength' 'jAlignLength','Length_V','Length_Nvd','Length_D','Length_Ndj','Length_J'});
H.DelLoc      = findCell(VDJheader,{'vDeletion','d5Deletion','d3Deletion','jDeletion','V_Deletion3','D_Deletion5','D_Deletion3','J_Deletion5'});

H.MatchLoc    = findCell(VDJheader,{'vMatchCount' 'dMatchCount', 'jMatchCount'});
H.SHMLoc      = findCell(VDJheader,{'SHM_V' 'SHM_D', 'SHM_J'});    

H.AlignLoc    = findCell(VDJheader,{'vAlignment','dAlignment','jAlignment','V_Alignment','D_Alignment','J_Alignment'});
H.FormClassLoc = findCell(VDJheader,{'FormattedSeq'; 'Classifier'});
H.ChildCountLoc = findCell(VDJheader,'TreeChildCount');
    
H.FunctLoc    = findCell(VDJheader,{'Functional'});

H.VmutLoc     = findCell(VDJheader,{'vMutCt_Germline','SHM_V'});
H.MmutLoc     = findCell(VDJheader,{'mMutCt_Germline','SHM_Nvd'});
H.DmutLoc     = findCell(VDJheader,{'dMutCt_Germline','SHM_D'});
H.NmutLoc     = findCell(VDJheader,{'nMutCt_Germline','SHM_Ndj'});
H.JmutLoc     = findCell(VDJheader,{'jMutCt_Germline','SHM_J'});

H.LeftTrimLoc     = findCell(VDJheader,'TrimmedSeqLeft');
H.RightTrimLoc    = findCell(VDJheader,'TrimmedSeqRight');
H.MiscLoc     = findCell(VDJheader,{'Misc'});
