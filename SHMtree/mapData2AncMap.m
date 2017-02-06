%mapData2AncMap will take Tdata (subset of VDJdata for same clonotype) and
%AncMap, and then map the RefSeq to ParentSeq;
%
%  Tdata = mapData2AncMap(Tdata,VDJheader,AncMap) %Will build refseq by
%  default using current database;
%  Tdata = mapData2AncMap(Tdata,VDJheader,AncMap,'SkipRef') %Will NOT build
%  refseq, since it assumes you already have it.

function Tdata = mapData2AncMap(Tdata,VDJheader,AncMap,varargin)
H = getHeaderVar(VDJheader);

%Does the sequence for [AncMap(x,2) == 0] already have the germline RefSeq?
SkipBuildRef = 0; %Assume no, so build RefSeq.
if ~isempty(varargin)
    if strcmpi(varargin{1},'SkipRef');
        SkipBuildRef = 1; %Skip finding germline seq
    end
end

%Make sure AncMap is relatively mapped by Tdata position. 
%EX: AncMap(1,1) --> Tdata(1,:)
AncMap = renumberAncMap(AncMap);

%For each child seq, you want the RefSeq to be the Seq of the parent.
for j = 1:size(Tdata,1)
    if AncMap(j,2) == 0;
        if SkipBuildRef == 0
            Tdata(AncMap(j,1),:) = buildRefSeq(Tdata(AncMap(j,1),:),VDJheader);
        end
    else
        Tdata{AncMap(j,1),H.RefSeqLoc} = Tdata{AncMap(j,2),H.SeqLoc};
    end
end
% 
% %Make sure the 1st one is the 0.
% AncMap = sortrows(AncMap,[2 1]);
% Tdata = Tdata(AncMap(:,1),:);
