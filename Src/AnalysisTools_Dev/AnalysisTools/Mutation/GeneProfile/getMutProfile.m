%getMutProfile gets the mutation profile per unique V, D, J gene per
%repertoire or repertoires.


function MutProfile = getMutProfile(varargin)
FileNames = openFileDialog('*.csv', 'Select the BRILIA files', 'multiselect', 'on');
DB = getGeneDatabase('mouse');

ProfV  = getProfileStruct(DB, 'V');
ProfD  = getProfileStruct(DB, 'D');
ProfJ  = getProfileStruct(DB, 'J');
TotalV = ProfV;
TotalD = ProfD;
TotalJ = ProfJ;

Vfield = fieldnames(ProfV);

for f = 1:length(FileNames)
    f
    [VDJdata, VDJheader] = openSeqData(FileNames{f});
    Map = getVDJmapper(VDJheader);
    GrpNum = cell2mat(VDJdata(:, Map.GrpNum));
    UnqGrpNum = unique(GrpNum);
    
    for y = 1:length(UnqGrpNum)
        Idx = find(UnqGrpNum(y) == GrpNum);
        VMDNJ = VDJdata{Idx(1), Map.hLength};
        RefSeq = VDJdata{Idx(1), Map.hRefSeq};
        CDR3bgn = VDJdata{Idx(1), Map.hCDR3(3)};
        ReadFrame = mod(CDR3bgn-1, 3) + 1;
        GeneNames = getCleanGeneName(VDJdata(Idx(1), Map.hGeneName));
        GeneNums = find(contains(Vfield, GeneNames{1}));
        if isempty(GeneNums); continue; end
        %GeneNums = cellfun(@(x) x(1), VDJdata(Idx(1), Map.hGeneNum));        
       
        ShmCt = zeros(1, length(RefSeq));
        for j = 1:length(Idx)
            Seq = VDJdata{Idx(1), Map.hSeq};
            MutInfo = getMutInfo(RefSeq, Seq, 'Frame', ReadFrame);
            MutIdx = [MutInfo.Idx];
            NonSynIdx = MutIdx([MutInfo.Type] == 'N');
            ShmCt(NonSynIdx) = ShmCt(NonSynIdx) + 1;
            SynIdx = MutIdx([MutInfo.Type] == 'S');
            ShmCt(SynIdx) = ShmCt(SynIdx) - 1;
        end
        
        Vdel = VDJdata{Idx(1), Map.hDel(1)};
        Vunused = length(DB.Vmap{GeneNums(1), 1}) - VMDNJ(1) - Vdel;
        CurVshm = ProfV.(GeneNames{1});
        CurVtotal = TotalV.(GeneNames{1});
        ProfV.(GeneNames{1})(Vunused+(1:VMDNJ(1))) = CurVshm(Vunused+(1:VMDNJ(1))) + ShmCt(1:VMDNJ(1));
        TotalV.(GeneNames{1})(Vunused+(1:VMDNJ(1))) = CurVtotal(Vunused+(1:VMDNJ(1))) + length(Idx);
    end
end



Fields = fieldnames(ProfV);
for c = 1:length(Fields)
    Ytotal = TotalV.(Fields{c});
    if max(Ytotal) < 20; continue; end

    Ydata = ProfV.(Fields{c});
    if c == 1        
        Bx = bar(Ydata);
        Tx = title(DB.Vmap{c, 2});
        resizeFigure(gca, 5, 2);
        xlim([200, 330]);
        ylim([-100, 100])
        resizeSubplots(gca);
    else
        Bx.YData = Ydata;
        Tx.String = DB.Vmap{c, 2};
    end
    
    savePlot(gcf, 'SaveAs', ['GeneNum_' num2str(c) '.png']);
end