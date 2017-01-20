%getHeaderVar is used to return the column number associated with the data
%type of VDJdata. It is NOT an enumeration class because newer version of
%VDJdata can have different headers, hence once must alway compute the new
%column locations.
%
%  getHeaderVar;

%Extract relevant column locations
SeqLoc      = findCell(NewHeader,{'nucleotide','Seq'});
SeqNumLoc   = findCell(NewHeader,'SeqNum');
SeqNameLoc  = findCell(NewHeader,'SeqName');

TemplateLoc = findCell(NewHeader,{'count (templates)','TemplateCount'});
CDR3Loc     = findCell(NewHeader,{'CDR3_AminoAcid','CDR3_Length','CDR3_Start','CDR3_End'});
%CDR3Loc     = findCell(NewHeader,{'aminoAcid', 'cdr3Length','CDR3_AminoAcid','CDR3_Length','CDR3_Start','CDR3_End'});

GrpNumLoc   = findCell(NewHeader,'GroupNum');
RefSeqLoc   = findCell(NewHeader,'RefSeq');
FamNumLoc   = findCell(NewHeader,{'vMapNum','dMapNum','jMapNum','V_MapNum','D_MapNum','J_MapNum'});
FamLoc      = findCell(NewHeader,{'vMaxResolved','dMaxResolved','jMaxResolved','V_GeneName','D_GeneName','J_GeneName'});

LengthLoc   = findCell(NewHeader,{'vAlignLength' 'n2AlignLength' 'dAlignLength' 'n1AlignLength' 'jAlignLength','Length_V','Length_Nvd','Length_D','Length_Ndj','Length_J'});
DelLoc      = findCell(NewHeader,{'vDeletion','d5Deletion','d3Deletion','jDeletion','V_Deletion3','D_Deletion5','D_Deletion3','J_Deletion5'});

MatchLoc    = findCell(NewHeader,{'vMatchCount' 'dMatchCount', 'jMatchCount'});
SHMLoc      = findCell(NewHeader,{'SHM_V' 'SHM_D', 'SHM_J'});    

AlignLoc    = findCell(NewHeader,{'vAlignment','dAlignment','jAlignment','V_Alignment','D_Alignment','J_Alignment'});
FormClassLoc    = findCell(NewHeader,{'FormattedSeq'; 'Classifier'});
ChildCountLoc   = findCell(NewHeader,'TreeChildCount');
    
FunctLoc    = findCell(NewHeader,{'Functional'});

VmutLoc     = findCell(NewHeader,{'vMutCt_Germline','SHM_V'});
MmutLoc     = findCell(NewHeader,{'mMutCt_Germline','SHM_Nvd'});
DmutLoc     = findCell(NewHeader,{'dMutCt_Germline','SHM_D'});
NmutLoc     = findCell(NewHeader,{'nMutCt_Germline','SHM_Ndj'});
JmutLoc     = findCell(NewHeader,{'jMutCt_Germline','SHM_J'});

LeftTrimLoc     = findCell(NewHeader,'TrimmedSeqLeft');
RightTrimLoc    = findCell(NewHeader,'TrimmedSeqRight');

MiscLoc     = findCell(NewHeader,{'Misc'});