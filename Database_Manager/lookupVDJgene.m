%lookupVDJgene will take a Max Resolved VDJ gene name and attempt to find
%the matchine sequence information from the
%IMGT_Mouse_VDJgene_AllStrain.mat database file.

function [GeneSeq, varargout] = lookupVDJgene(GeneName, varargin)
%Extract the VDJ database
if length(varargin) == 3
    Vmap = varargin{1};
    Dmap = varargin{2};
    Jmap = varargin{3};
else
    [Vmap, Dmap, Jmap] = getCurrentDatabase;
end

%Determine the gene class
GeneClass = regexp(GeneName,'IGH[VDJ]','match');

%Make sure to account for  unknown names.
if isempty(GeneClass)
    GeneSeq = '';
    if nargout == 2;
        varargout{1} = 0;
    end
    return
end

switch GeneClass{1}
    case 'IGHV'
        Xmap = Vmap;
    case 'IGHD'
        Xmap = Dmap;
    case 'IGHJ'
        Xmap = Jmap;
end

%Look for the allele info, if any
StarLoc = regexp(GeneName,'*','once');
DashLoc = regexp(GeneName,'-','once');
if isempty(StarLoc) == 0 %Allele info provided, lookup full name column
    LookUpCol = [2 3]; %Check both my and IMGT formats
elseif isempty(DashLoc) == 0 %Up to gene name provided. (Should be common)
    LookUpCol = 5;
else %Up to family name provided
    LookUpCol = 4;
end

%Look for the match locations
MatchLoc = zeros(size(Xmap,1),length(LookUpCol));
for k = 1:length(LookUpCol)
    for j = 1:size(Xmap,1)
        if strcmpi(Xmap{j,LookUpCol(k)},GeneName)
            MatchLoc(j,k) = 1;
        end
    end
end
MatchLoc = sum(MatchLoc,2)>0;

if max(MatchLoc) == 0 %In case there is no match
    GeneSeq = '';
    if nargout == 2;
        varargout{1} = 0;
    end
else     %Output the selected sequence
    GeneSeq = Xmap(MatchLoc,1);
    if nargout == 2 %Also outputs the index location in the gene map
        varargout{1} = find(MatchLoc == 1);
    end    
end

