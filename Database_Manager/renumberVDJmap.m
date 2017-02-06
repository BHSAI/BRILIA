%renumberVDJmap will take renumber all V_MapNum, D_MapNum, and J_MapNum
%based on the full gene name. Use this when changing/updating databases.
%
%  VDJdata = renumberVDJmap(VDJdata,VDJheader)
%
%  VDJdata = renumberVDJmap(VDJdata,VDJheader,Vmap,Dmap,Jmap)

function VDJdata = renumberVDJmap(VDJdata,VDJheader,varargin)
%Extract the VDJ database
if length(varargin) == 3
    [Vmap,Dmap,Jmap] = deal(varargin{:});
else
    [Vmap, Dmap, Jmap] = getCurrentDatabase;
end
H = getHeaderVar(VDJheader);

%Renumber the V, D, J reference gene map number
for j = 1:size(VDJdata,1)
    for q = 1:3
        %Extract the gene name(s)
        Xfam = VDJdata{j,H.FamLoc(q)};
        if sum(isempty(Xfam))>0 || sum(isnan(Xfam)) > 0; 
            continue
        end
        
        %Remove the general name (ex: "IGHV01:")
        ColonLoc = regexp(Xfam,':');
        if ~isempty(ColonLoc)
            Xfam = Xfam(ColonLoc+1:end);
            Xfam = strrep(Xfam,' ','');
        end
        
        %Get the individual gene names, and reassign the ref number
        XfamAll = regexp(Xfam,'\|','split');
        XmapNum = zeros(1,length(XfamAll));
        for v = 1:length(XfamAll)
            [~, XmapNum(v)] = lookupVDJgene(XfamAll{v},Vmap,Dmap,Jmap);
        end
        
        VDJdata{j,H.FamNumLoc(q)} = XmapNum;
    end
end
