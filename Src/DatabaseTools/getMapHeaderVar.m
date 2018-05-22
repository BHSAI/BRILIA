%getMapHeaderVar is used to return the column number associated with the
%data type of Xmap via a short naming system and ending with "Loc". It is
%NOT an enumeration class because newer version can have different headers,
%hence once must alway compute the new column locations.
%
%  M = getMapHeaderVar(MapHeader); 
%
%  INPUT
%    MapHeader: 1xM cell headers for the germline sequences.
%
%  OUTPUT
%    M: structure containing standardized header name variable
%      M.Seq
%      M.Gene
%      M.Funct
%      M.Strain
%      M.EntryNum
%      M.GapInfo
%      M.CDR1s
%      M.CDR1e
%      M.CDR2s
%      M.CDR2e
%      M.Anchor
function M = getMapHeaderVar(MapHeader)
%Extract relevant column locations
M.Seq      = findCell(MapHeader,'Seq');
M.Gene     = findCell(MapHeader,'GeneName');
M.Funct    = findCell(MapHeader,'Function');
M.Strain   = findCell(MapHeader,'Strain');
M.EntryNum = findCell(MapHeader,'EntryNum');
M.GapInfo  = findCell(MapHeader,'GapInfo');
M.CDR1s    = findCell(MapHeader,'CDR1start');
M.CDR1e    = findCell(MapHeader,'CDR1end');
M.CDR2s    = findCell(MapHeader,'CDR2start');
M.CDR2e    = findCell(MapHeader,'CDR2end');
M.Anchor   = findCell(MapHeader,'AnchorDist');
