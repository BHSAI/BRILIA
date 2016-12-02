function VDJdata = updateVDJdata(VDJdata,NewHeader,varargin)
%Extract the VDJ database
if length(varargin) == 3
    Vmap = varargin{1};
    Dmap = varargin{2};
    Jmap = varargin{3};
else
    [Vmap, Dmap, Jmap] = getCurrentDatabase;
end

%Fill in the details now
VDJdata = buildRefSeq(VDJdata,NewHeader,'germline','first'); %must do singles, since group therapy not done.
VDJdata = buildVDJalignment(VDJdata,NewHeader,Vmap,Dmap,Jmap); %Alignment Info
VDJdata = makeClassifier(VDJdata,NewHeader); %Classifier + FormattedSeq
VDJdata = appendMutCt(VDJdata,NewHeader); %SHM infor on the VMDNJ segments
VDJdata = findCDR3(VDJdata,NewHeader); %Get the CDR3 seq and info 
VDJdata = labelNonprodVDJ(VDJdata,NewHeader);