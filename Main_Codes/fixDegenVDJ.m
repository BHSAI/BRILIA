%fixDegenVDJ will remove pseudogene (P) and open reading frame (ORF)
%annotations if a functional gene annotations are also suggested by BRILIA.
%
%  VDJdata = fixPseudoVDJ(VDJdata,NewHeader)
%
%  VDJdata = fixPseudoVDJ(VDJdata,NewHeader,Vmap,Dmap,Jmap)
%
%  VDJdata = fixPseudoVDJ(VDJdata,NewHeader,Vmap,Dmap,Jmap,VDJ)
%
%  INPUT
%    VDJ ['V','D','J'] or any combo: Specifies which pseudogene to remove.
%    Default is VDJ = 'VDJ' to do all V, D, J genes. 
%
%  See also BRILIA
function VDJdata = fixDegenVDJ(VDJdata,NewHeader,varargin)
%Extract the V database
if length(varargin) >= 3 
    [Vmap,Dmap,Jmap] = deal(varargin{1:3});
else
    [Vmap,Dmap,Jmap] = getCurrentDatabase;
end       
getHeaderVar;

%Header for the Vmap;
MapFunctLoc = 7;
MapStdNameLoc = 3;

%Determine which gene to perform this on
VDJ = 'VDJ'; %default
if length(varargin) == 4
    VDJ = varargin{4};
end

UpdateIdx = zeros(size(VDJdata,1),1,'logical');
for x = 1:length(VDJ)
    switch upper(VDJ(x))
        case 'V'
            FamLocX = FamLoc(1);
            FamNumLocX = FamNumLoc(1);
            Xmap = Vmap;
        case 'D'
            FamLocX = FamLoc(2);
            FamNumLocX = FamNumLoc(2);
            Xmap = Dmap;
        case 'J'
            FamLocX = FamLoc(3);
            FamNumLocX = FamNumLoc(3);
            Xmap = Jmap;
    end

    %Identify location of pseudo genes
    PseudoLoc = zeros(size(Xmap,1),1,'logical');
    for v = 1:size(Xmap,1)
        if ~isempty(regexpi(Xmap{v,MapFunctLoc},'P|ORF'))
            PseudoLoc(v) = 1;
        end
    end

    %Iteratively find groups
    GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
    UnqGrpNum = unique(GrpNum);
    for y = 1:length(UnqGrpNum)
        try
            %Find the sequences per grp, and the Xnum
            IdxLoc = find(UnqGrpNum(y) == GrpNum);   
            Xnum = VDJdata{IdxLoc(1),FamNumLocX};
            if length(Xnum) <= 1; continue; end

            %Look for any pseudogene
            DelThis = zeros(length(Xnum),1,'logical');
            for k = 1:length(Xnum)
                if PseudoLoc(Xnum) == 1
                    DelThis(k) = 1;
                end
            end
            HavePseudo = max(DelThis);
            HowManyPseudo = sum(DelThis);

            if HavePseudo %Figure out what to do
                if HowManyPseudo >= 1 && (HowManyPseudo < length(Xnum)) %Must have at least 1 instance, but not all instances
                    %Need to adjust the name
                    KeepThis = ~DelThis;
                    TempCell = cell(sum(KeepThis),6); %Making temp cell so that reduceFamily works
                    TempCell(:,1) = num2cell(Xnum(KeepThis)); %Get the new Vnums w/o pseudogenes
                    TempCell(:,2) = Xmap(Xnum(KeepThis),MapStdNameLoc); %Save the gene names
                    TempResults = reduceFamily(TempCell); %Reduce the names into single line.

                    %Update VDJdata family number and name
                    VDJdata(IdxLoc,[FamNumLocX FamLocX]) = repmat(TempResults(1,1:2),length(IdxLoc),1);
                    UpdateIdx(IdxLoc) = 1;
                end
            end
        catch
            WarningMsg = sprintf('Warning at %s, sequence group # %d',mfilename,UnqGrpNum(y));
            disp(WarningMsg);
        end
    end
end

%Update those that have changed
VDJdata(UpdateIdx,:) = buildRefSeq(VDJdata(UpdateIdx,:),NewHeader,'germline','first'); %must do first seq of all cluster
VDJdata(UpdateIdx,:) = updateVDJdata(VDJdata(UpdateIdx,:),NewHeader,varargin);