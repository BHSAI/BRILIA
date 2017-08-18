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
%      M.SeqLoc
%      M.GeneLoc
%      M.FunctLoc
%      M.StrainLoc
%      M.EntryNumLoc
%      M.GapInfoLoc
%      M.CDR1sLoc
%      M.CDR1eLoc
%      M.CDR2sLoc
%      M.CDR2eLoc
%      M.AnchorLoc
function M = getMapHeaderVar(MapHeader)
%Extract relevant column locations
M.SeqLoc      = findCell(MapHeader,'Seq');
M.GeneLoc     = findCell(MapHeader,'GeneName');
M.FunctLoc    = findCell(MapHeader,'Function');
M.StrainLoc   = findCell(MapHeader,'Strain');
M.EntryNumLoc = findCell(MapHeader,'EntryNum');
M.GapInfoLoc  = findCell(MapHeader,'GapInfo');
M.CDR1sLoc    = findCell(MapHeader,'CDR1start');
M.CDR1eLoc    = findCell(MapHeader,'CDR1end');
M.CDR2sLoc    = findCell(MapHeader,'CDR2start');
M.CDR2eLoc    = findCell(MapHeader,'CDR2end');
M.AnchorLoc   = findCell(MapHeader,'AnchorDist');
