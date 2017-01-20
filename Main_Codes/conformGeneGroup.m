%conformGeneGroup will force all sequences in the same groups to have the
%same Gene family and VMDNJ segment lengths as the first sequence in the
%cluster, which is suppose to be close to the germline sequence.
%
%  VDJdata = conformGeneGroup(VDJdata,NewHeader)
%
function VDJdata = conformGeneGroup(VDJdata,NewHeader,varargin)
%Extract the VDJ database
if length(varargin) >= 3
    [Vmap,Dmap,Jmap] = deal(varargin{1:3});
else
    [Vmap,Dmap,Jmap] = getCurrentDatabase;
end
getHeaderVar;

%Perform 1st group therapy: Conforming VMDNJ lengths to group
GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
UnqGrpNum = unique(GrpNum);
for y = 1:length(UnqGrpNum)
    try
        %Determine if group correction is possible
        IdxLoc = find(UnqGrpNum(y) == GrpNum);
        if length(IdxLoc) <= 1; continue; end %Nothing to correct, single entry.

        %Standardize the other fields as well (1st one is the root)
        Tdata = VDJdata(IdxLoc,:);    
        Tdata(:,[LengthLoc(:); FamLoc(:); FamNumLoc(:); DelLoc(:)]) = repmat(Tdata(1,[LengthLoc(:); FamLoc(:); FamNumLoc(:); DelLoc(:)]),size(Tdata,1),1);

        VDJdata(IdxLoc,:) = Tdata;  
    catch
        WarningMsg = sprintf('Warning at %s, sequence group # %d',mfilename,UnqGrpNum(y));
        disp(WarningMsg);
    end
end

%Update all data
VDJdata = buildRefSeq(VDJdata,NewHeader,'germline','first'); %must do first seq of each cluster
VDJdata = updateVDJdata(VDJdata,NewHeader,varargin);

