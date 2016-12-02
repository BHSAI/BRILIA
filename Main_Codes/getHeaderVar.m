%Extract relevant column locations
SeqLoc = findHeader(NewHeader,{'nucleotide','Seq'});
if ~isempty(SeqLoc) && sum(SeqLoc)>0
    SeqLoc(SeqLoc==0) = []; SeqLoc = SeqLoc(1); 
end
SeqNumLoc = findHeader(NewHeader,'SeqNum');

SeqNameLoc = findHeader(NewHeader,'SeqName');

GrpNumLoc = findHeader(NewHeader,'GroupNum');

FamNumLoc = findHeader(NewHeader,{'vMapNum','dMapNum','jMapNum','V_MapNum','D_MapNum','J_MapNum'});
    FamNumLoc(FamNumLoc==0) = [];

FamLoc = findHeader(NewHeader,{'vMaxResolved','dMaxResolved','jMaxResolved','V_GeneName','D_GeneName','J_GeneName'});
    FamLoc(FamLoc==0) = []; 

LengthLoc = findHeader(NewHeader,{'vAlignLength' 'n2AlignLength' 'dAlignLength' 'n1AlignLength' 'jAlignLength','Length_V','Length_Nvd','Length_D','Length_Ndj','Length_J'});
    LengthLoc(LengthLoc==0) = [];

DelLoc = findHeader(NewHeader,{'vDeletion','d5Deletion','d3Deletion','jDeletion','V_Deletion3','D_Deletion5','D_Deletion3','J_Deletion5'});
    DelLoc(DelLoc==0) = [];

MatchLoc = findHeader(NewHeader,{'vMatchCount' 'dMatchCount', 'jMatchCount'});

SHMLoc = findHeader(NewHeader,{'SHM_V' 'SHM_D', 'SHM_J'});    

AlignLoc = findHeader(NewHeader,{'vAlignment','dAlignment','jAlignment','V_Alignment','D_Alignment','J_Alignment'});
    AlignLoc(AlignLoc==0) = [];

FormClassLoc = findHeader(NewHeader,{'FormattedSeq'; 'Classifier'});

RefSeqLoc = findHeader(NewHeader,'RefSeq');

CDR3Loc = findHeader(NewHeader,{'aminoAcid', 'cdr3Length','CDR3_AminoAcid','CDR3_Length','CDR3_Start','CDR3_End'});
    if sum(CDR3Loc(1:2)) == 0
        CDR3Loc = CDR3Loc(3:end);
    end
    
TemplateLoc = findHeader(NewHeader,{'count (templates)','TemplateCount'});
    TemplateLoc(TemplateLoc==0) = [];
    
ChildCountLoc = findHeader(NewHeader,'TreeChildCount');
    
FunctLoc = findHeader(NewHeader,{'Functional'});

VmutLoc = findHeader(NewHeader,{'vMutCt_Germline'});
MmutLoc = findHeader(NewHeader,{'mMutCt_Germline'});
DmutLoc = findHeader(NewHeader,{'dMutCt_Germline'});
NmutLoc = findHeader(NewHeader,{'nMutCt_Germline'});
JmutLoc = findHeader(NewHeader,{'jMutCt_Germline'});
