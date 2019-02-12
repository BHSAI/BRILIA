function ChargeChanges = calcChargeChanges(VDJdata, VDJheader)

Map = getVDJmapper(VDJheader);
G = getGrpIdx(VDJdata, VDJheader);

CDRs = cell(1, 3*length(Map.Chain));
for c = 1:length(Map.Chain)
    for k = 1:3
        CDRs{(c-1)+k} = sprintf('%sCDR%d', lower(Map.Chain(c)), k);
    end
end
TempStruct = [CDRs; repmat({{cell(length(G), 2)}}, 1, length(CDRs))];
ChargeChanges = struct(TempStruct{:});

for k = 1:length(CDRs)
    for y = 1:length(G)
        
        if length(G(y).Idx) == 1; continue; end

        if min(VDJdata{G(y).Idx(1), Map.(CDRs{k})(3:4)}) <= 0 %Skip if no CDR annotated
            continue
        end
            
        CDRbgn = VDJdata{G(y).Idx(1), Map.(CDRs{k})(3)};
        CDRend = VDJdata{G(y).Idx(1), Map.(CDRs{k})(4)};
        DeltaCharge = zeros(length(G(y).Idx), 1);
        for j = 1:length(G(y).Idx)
            RefSeq = VDJdata{G(y).Idx(j), Map.hRefSeq};
            RefCDR = nt2aa(RefSeq(CDRbgn:CDRend), 'acgtonly', false);
            RefProp = convAA2PropMEX(RefCDR);%;, ReducedLetter);

            CDR = VDJdata{G(y).Idx(j), Map.(CDRs{k})};
            Prop = convAA2PropMEX(CDR);%, ReducedLetter);
            if strcmp(RefProp, Prop); continue; end % Same CDR, no change. skip.
            
            [~, DeltaCharge(j)] = convProp2Charge(RefProp, Prop);
        end
        ChargeChanges.(CDRs{k})(y, :) = {length(RefCDR) DeltaCharge(DeltaCharge ~= 0)};        
    end
end




