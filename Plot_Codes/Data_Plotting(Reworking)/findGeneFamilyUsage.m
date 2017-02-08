%findGeneFamilyUsage will extract the family number for V,D,J  usage.
%Unresolved/empty = 0, and inverse gene matches are negative values.

%Option = 1, individual gene usage data
%Option = 2, group by gene family

function [VDJcombo, varargout] = findGeneFamilyUsage(VDJdata,VDJheader,Option,varargin)
%Extract the VDJ database
if length(varargin) == 3
    Vmap = varargin{1};
    Dmap = varargin{2};
    Jmap = varargin{3};
else
    [Vmap, Dmap, Jmap] = getCurrentDatabase;
end

H = getHeaderVar(VDJheader);

if Option == 1
    %Make the 3D coordinate matrix
    VDJcombo = zeros(size(Vmap,1),size(Dmap,1),size(Jmap,1));

    %Get all the gene usage. Degenerate matches contribute fractions. 
    GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
    UnqGrpNum = unique(GrpNum);
    for y = 1:length(UnqGrpNum)
        IdxLoc = find(UnqGrpNum(y) == GrpNum);
        IdxLoc = IdxLoc(1); %You only need 1 of the cluter data.
        Vnum = VDJdata{IdxLoc,H.FamNumLoc(1)};
        Dnum = VDJdata{IdxLoc,H.FamNumLoc(2)};
        Jnum = VDJdata{IdxLoc,H.FamNumLoc(3)};
        if isempty(Vnum) || isempty(Dnum) || isempty(Jnum)
            continue
        elseif Vnum(1)*Dnum(1)*Jnum(1) == 0
            continue
        elseif isnan(Vnum) 
            continue
        elseif isnan(Dnum)
            continue
        elseif isnan(Jnum)
            continue
        end

        Weight = 1/(length(Vnum)*length(Dnum)*length(Jnum));

        for v = 1:length(Vnum)
            for d = 1:length(Dnum)
                for j = 1:length(Jnum)
                    VDJcombo(Vnum(v),Dnum(d),Jnum(j)) = VDJcombo(Vnum(v),Dnum(d),Jnum(j)) + Weight;
                end
            end
        end
    end
    if nargout > 1
        AxisLabels{1} = Vmap(:,2);
        AxisLabels{2} = Dmap(:,2);
        AxisLabels{3} = Jmap(:,2);
        varargout{1} = AxisLabels;
    end

elseif Option == 2
    %If Option = 2, need to regroup by family
    %Create the regroup map

    [VunqFam, ~, Vidx] = unique(Vmap(:,4));
    [DunqFam, ~, Didx] = unique(Dmap(:,4));
    [JunqFam, ~, Jidx] = unique(Jmap(:,4));

    %Make the 3D coordinate matrix
    VDJcombo = zeros(length(VunqFam),length(DunqFam),length(JunqFam));

    %Get all the gene usage. Degenerate matches contribute fractions. 
    GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
    UnqGrpNum = unique(GrpNum);
    for y = 1:length(UnqGrpNum)
        IdxLoc = find(UnqGrpNum(y) == GrpNum);
        IdxLoc = IdxLoc(1); %You only need 1 of the cluter data.
        Vnum = VDJdata{IdxLoc,H.FamNumLoc(1)};
        Dnum = VDJdata{IdxLoc,H.FamNumLoc(2)};
        Jnum = VDJdata{IdxLoc,H.FamNumLoc(3)};

        if isempty(Vnum) || isempty(Dnum) || isempty(Jnum)
            continue
        end
                
        Weight = 1/(length(Vnum)*length(Dnum)*length(Jnum));

        for v = 1:length(Vnum)
            for d = 1:length(Dnum)
                for j = 1:length(Jnum)
                    VDJcombo(Vidx(Vnum(v)),Didx(Dnum(d)),Jidx(Jnum(j))) = VDJcombo(Vidx(Vnum(v)),Didx(Dnum(d)),Jidx(Jnum(j))) + Weight;
                end
            end
        end
    end
    if nargout > 1
        AxisLabels{1} = VunqFam;
        AxisLabels{2} = DunqFam;
        AxisLabels{3} = JunqFam;
        varargout{1} = AxisLabels;
    end
    
else
    %If Option = 3, need to regroup by family and number, but not allele
    [VunqFam, ~, Vidx] = unique(Vmap(:,5));
    [DunqFam, ~, Didx] = unique(Dmap(:,5));
    [JunqFam, ~, Jidx] = unique(Jmap(:,5));

    %Make the 3D coordinate matrix
    VDJcombo = zeros(length(VunqFam),length(DunqFam),length(JunqFam));

    %Get all the gene usage. Degenerate matches contribute fractions. 
    GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
    UnqGrpNum = unique(GrpNum);
    for y = 1:length(UnqGrpNum)
        IdxLoc = find(UnqGrpNum(y) == GrpNum);
        IdxLoc = IdxLoc(1); %You only need 1 of the cluter data.
        Vnum = VDJdata{IdxLoc,H.FamNumLoc(1)};
        Dnum = VDJdata{IdxLoc,H.FamNumLoc(2)};
        Jnum = VDJdata{IdxLoc,H.FamNumLoc(3)};

        if isempty(Vnum) || isempty(Dnum) || isempty(Jnum)
            continue
        end
        
        Weight = 1/(length(Vnum)*length(Dnum)*length(Jnum));

        for v = 1:length(Vnum)
            for d = 1:length(Dnum)
                for j = 1:length(Jnum)
                    VDJcombo(Vidx(Vnum(v)),Didx(Dnum(d)),Jidx(Jnum(j))) = VDJcombo(Vidx(Vnum(v)),Didx(Dnum(d)),Jidx(Jnum(j))) + Weight;
                end
            end
        end
    end
    if nargout > 1
        AxisLabels{1} = VunqFam;
        AxisLabels{2} = DunqFam;
        AxisLabels{3} = JunqFam;
        varargout{1} = AxisLabels;
    end
end
