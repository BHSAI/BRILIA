%getTreeHeaderVar is used to return the column number associated with the
%data type of TreeData. It is NOT an enumeration class because newer
%version of BRILIA can change the header locations.
%
%  getTreeHeaderVar;

%Extract relevant column locations
ChildSeqNumLoc  = findCell(TreeHeader,'ChildSeqNum'); 
ParentSeqNumLoc = findCell(TreeHeader,'ParentSeqNum');
TemplateLoc  = findCell(TreeHeader,'TemplateCount');
SHMdistLoc   = findCell(TreeHeader,'SHMdist');
HAMdistLoc   = findCell(TreeHeader,'HAMdist');
CDR3seqLoc   = findCell(TreeHeader,'CDR3seq');
GrpNameLoc   = findCell(TreeHeader,'GrpName');
GrpNumLoc    = findCell(TreeHeader,'GrpNum');
SeqNumLoc    = findCell(TreeHeader,'ChildSeqNum'); %Same as ChildSeqNumLoc
