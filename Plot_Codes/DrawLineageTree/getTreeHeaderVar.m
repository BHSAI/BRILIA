%getTreeHeaderVar is used to return the column number associated with the
%data type of TreeData. It is NOT an enumeration class because newer
%version of BRILIA can change the header locations.
%
%  TH = getTreeHeaderVar(TreeHeader);
function TH = getTreeHeaderVar(TreeHeader)
%Extract relevant column locations
TH.ChildSeqNumLoc  = findCell(TreeHeader,'ChildSeqNum'); 
TH.ParentSeqNumLoc = findCell(TreeHeader,'ParentSeqNum');
TH.TemplateLoc  = findCell(TreeHeader,'TemplateCount');
TH.SHMdistLoc   = findCell(TreeHeader,'SHMdist');
TH.HAMdistLoc   = findCell(TreeHeader,'HAMdist');
TH.CDR3seqLoc   = findCell(TreeHeader,'CDR3seq');
TH.GrpNameLoc   = findCell(TreeHeader,'GrpName');
TH.GrpNumLoc    = findCell(TreeHeader,'GrpNum');
TH.SeqNumLoc    = findCell(TreeHeader,'ChildSeqNum'); %Same as TH.ChildSeqNumLoc
