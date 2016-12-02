%renumberVDJmap will take VDJdata, and renumber all the vMapNum, dMapNum,
%jMapNum, based on the full gene name. Use this when changing/updating
%databases.

function VDJdata = renumberVDJmap(VDJdata,NewHeader,varargin)
%Extract the VDJ database
if length(varargin) == 3
    Vmap = varargin{1};
    Dmap = varargin{2};
    Jmap = varargin{3};
else
    [Vmap, Dmap, Jmap] = getCurrentDatabase;
end

getHeaderVar;

for j = 1:size(VDJdata,1)
    for q = 1:3
        Xfam = VDJdata{j,FamLoc(q)};
        if sum(isempty(Xfam))>0 || sum(isnan(Xfam)) > 0; 
            continue
        end
        
        ColonLoc = regexp(Xfam,':');
        if ~isempty(ColonLoc)
            Xfam = Xfam(ColonLoc+1:end);
            Xfam = strrep(Xfam,' ','');
        end
        
        XfamAll = regexp(Xfam,'\|','split');
        XmapNum = zeros(1,length(XfamAll));
        for v = 1:length(XfamAll)
            [~, XmapNum(v)] = lookupVDJgene(XfamAll{v},Vmap,Dmap,Jmap);
        end
        VDJdata{j,FamNumLoc(q)} = XmapNum;
    end
    if mod(j,1000) == 0
        [j/size(VDJdata,1)]
    end
end