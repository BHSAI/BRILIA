%findGeneFamilyUsage will extract the family number for V,D,J  usage.
%Unresolved/empty = 0, and inverse gene matches are negative values.

%Option = 1, individual gene usage data
%Option = 2, group by gene family

function [Dcombo, varargout] = findGeneFamilyUsage(VDJdata,VDJheader,Option,varargin)
%Extract the VDJ database
if length(varargin) == 3
    Vmap = varargin{1};
    Dmap = varargin{2};
    Jmap = varargin{3};
else
    [Vmap, Dmap, Jmap] = getCurrentDatabase;
    [Vmap, Dmap, Jmap] = filterRefGene(Vmap,Dmap,Jmap); %Selelect host strain
end

H = getHeaderVar(VDJheader);

if Option == 1
    %Make the 3D coordinate matrix
    Dcombo = zeros(size(Dmap,1),3);

    %Get all the gene usage. Degenerate matches contribute fractions. 
    GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
    UnqGrpNum = unique(GrpNum);
    for y = 1:length(UnqGrpNum)
        IdxLoc = find(UnqGrpNum(y) == GrpNum);
        IdxLoc = IdxLoc(1); %You only need 1 of the cluter data.
        Dnum = VDJdata{IdxLoc,H.FamNumLoc(2)};
        if isempty(Dnum)
            continue
        end
        
        Weight = 1/length(Dnum);
        
        %Determine the reading frame of the D gene
        [VMDNJ] = cell2mat(VDJdata(IdxLoc,H.LengthLoc));
        D5del = VDJdata{IdxLoc,H.DelLoc(2)};
        Dstart = VMDNJ(1) + VMDNJ(2) + 1 - D5del;
        CDR3start = VDJdata{IdxLoc,H.CDR3Loc(3)};
        RFloc = mod(Dstart - CDR3start + 1,3);
        if RFloc == 0; RFloc = 3; end

        for d = 1:length(Dnum)
            Dcombo(Dnum(d),RFloc) = Dcombo(Dnum(d),RFloc) + Weight;
        end
    end
    
%     %Remove the empties
%     DelThis = sum(Dcombo,2) == 0;
%     Dcombo(DelThis,:) = [];
%     Dmap(DelThis,:) = [];
        
    if nargout > 1
        varargout{1} = Dmap(:,2);
    end

elseif Option == 2
    %If Option = 2, need to regroup by family
    %Create the regroup map

    [DunqFam, ~, Didx] = unique(Dmap(:,4));

    %Make the 3D coordinate matrix
    Dcombo = zeros(length(DunqFam),3);

    %Get all the gene usage. Degenerate matches contribute fractions. 
    GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
    UnqGrpNum = unique(GrpNum);
    for y = 1:length(UnqGrpNum)
        IdxLoc = find(UnqGrpNum(y) == GrpNum);
        IdxLoc = IdxLoc(1); %You only need 1 of the cluter data.
        Dnum = VDJdata{IdxLoc,H.FamNumLoc(2)};

        if isempty(Dnum)
            continue
        end
                
        Weight = 1/length(Dnum);
        
        %Determine the reading frame of the D gene
        [VMDNJ] = cell2mat(VDJdata(IdxLoc,H.LengthLoc));
        D5del = VDJdata{IdxLoc,H.DelLoc(2)};
        Dstart = VMDNJ(1) + VMDNJ(2) + 1 - D5del;
        CDR3start = VDJdata{IdxLoc,H.CDR3Loc(3)};
        RFloc = mod(Dstart - CDR3start + 1,3);
        if RFloc == 0; RFloc = 3; end
        
        for d = 1:length(Dnum)
            Dcombo(Didx(Dnum(d)),RFloc) =  Dcombo(Didx(Dnum(d)),RFloc) + Weight;
        end
    end
    
%     %Remove the empties
%     DelThis = sum(Dcombo,2) == 0;
%     Dcombo(DelThis,:) = [];
%     DunqFam(DelThis,:) = [];

    if nargout > 1
        varargout{1} = DunqFam;
    end
end
%     
% else
%     %If Option = 3, need to regroup by family and number, but not allele
%     [VunqFam, ~, Vidx] = unique(Vmap(:,5));
%     [DunqFam, ~, Didx] = unique(Dmap(:,5));
%     [JunqFam, ~, Jidx] = unique(Jmap(:,5));
% 
%     %Make the 3D coordinate matrix
%     Dcombo = zeros(length(VunqFam),length(DunqFam),length(JunqFam));
% 
%     %Get all the gene usage. Degenerate matches contribute fractions. 
%     GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
%     UnqGrpNum = unique(GrpNum);
%     for y = 1:length(UnqGrpNum)
%         IdxLoc = find(UnqGrpNum(y) == GrpNum);
%         IdxLoc = IdxLoc(1); %You only need 1 of the cluter data.
%         Vnum = VDJdata{IdxLoc,H.FamNumLoc(1)};
%         Dnum = VDJdata{IdxLoc,H.FamNumLoc(2)};
%         Jnum = VDJdata{IdxLoc,H.FamNumLoc(3)};
% 
%         if isempty(Vnum) || isempty(Dnum) || isempty(Jnum)
%             continue
%         end
%         
%         Weight = 1/(length(Vnum)*length(Dnum)*length(Jnum));
% 
%         for v = 1:length(Vnum)
%             for d = 1:length(Dnum)
%                 for j = 1:length(Jnum)
%                     Dcombo(Vidx(Vnum(v)),Didx(Dnum(d)),Jidx(Jnum(j))) = Dcombo(Vidx(Vnum(v)),Didx(Dnum(d)),Jidx(Jnum(j))) + Weight;
%                 end
%             end
%         end
%     end
%     if nargout > 1
%         AxisLabels{1} = VunqFam;
%         AxisLabels{2} = DunqFam;
%         AxisLabels{3} = JunqFam;
%         varargout{1} = AxisLabels;
%     end
% end
