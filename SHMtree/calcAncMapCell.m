%calcAncMapCell will take VDJdata that has already the parent=child
%location made up, and return the AncMapCell variable required to make the
%phylogenetic trees.
%
%  AncMapCell = calcAncMapCell(VDJdata,NewHeader)

function AncMapCell = calcAncMapCell(VDJdata,NewHeader,varargin)

OutputDistMode = 'shmham';
if ~isempty(varargin)
    OutputDistMode = varargin{1};
end
getHeaderVar;

GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
UnqGrpNum = unique(GrpNum);
AncMapCell = cell(length(UnqGrpNum),2);
for y = 1:length(UnqGrpNum)
    IdxLoc = UnqGrpNum(y) == GrpNum;
    Tdata = VDJdata(IdxLoc,:);
    
    AncMapT = zeros(size(Tdata,1),4);
    AncMapT(:,1) = 1:size(Tdata,1);
    AncMapT(:,4) = cell2mat(Tdata(:,TemplateLoc));
    
    RefSeq = char(Tdata(:,RefSeqLoc));
    SamSeq = char(Tdata(:,SeqLoc));
    
    for j = 1:size(Tdata,1)
        for k = 1:size(Tdata,1)
            if sum(RefSeq(j,:) ~= SamSeq(k,:)) == 0 || j == 1%Complete match required
                AncMapT(j,2) = k;
                if strcmpi(OutputDistMode,'ham')
                    AncMapT(j,3) = sum(RefSeq(j,:) ~= SamSeq(j,:));
                elseif strcmpi(OutputDistMode,'shmham')
                    AncMapT(j,3) =  calcSHMHAMdist(RefSeq(j,:),SamSeq(j,:));
                end
                break
            end
        end
    end
    
    AncMapT(1,2) = 0;
    AncMapCell{y,1} = AncMapT;
    AncMapCell{y,2} = Tdata(AncMapT(:,1),CDR3Loc(1));    
end