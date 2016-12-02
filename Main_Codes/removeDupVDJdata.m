%removeDupVDJdata will look for duplicate nts, and then remove them. This
%could arise after fixIndel, where 2 unique sequences converge to be the
%same.

function VDJdata = removeDupVDJdata(VDJdata,NewHeader)

getHeaderVar;

SamSeq = VDJdata(:,SeqLoc);

%Check for duplicate nts (rare, but possible due to indel issues)
[UnqNT,~,UnqNTidx] = unique(SamSeq);

DelIdx = zeros(size(VDJdata,1),1) > 1;
if max(UnqNTidx) ~= length(SamSeq) %Duplicates found
    for q = 1:length(UnqNT)
        SubIdx = find(UnqNTidx == q);
        if length(SubIdx) > 1
            TempCt = cell2mat(VDJdata(SubIdx,TemplateLoc));
            MaxLoc = find(TempCt == max(TempCt));
            if isempty(MaxLoc)
                continue
            end
            MaxLoc = MaxLoc(1);
            MaxIdx = SubIdx(MaxLoc);
            TotTemp = sum(TempCt);
            VDJdata{MaxIdx,TemplateLoc} = TotTemp;
            SubIdx(MaxLoc) = [];
            DelIdx(SubIdx) = 1;           
        end
    end
end

VDJdata(DelIdx,:) = [];
        