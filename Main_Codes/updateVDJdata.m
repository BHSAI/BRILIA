%updateVDJdata will perform simple updates to VDJ data, EXCEPT for building
%the reference sequences as this requires some user input.
%
%  VDJdata = updateVDJdata(VDJdata,NewHeader)
%
%  VDJdata = updateVDJdata(VDJdata,NewHeader,Vmap,Dmap,Jmap)

function VDJdata = updateVDJdata(VDJdata,NewHeader,varargin)

VDJdata = buildVDJalignment(VDJdata,NewHeader,varargin); %Alignment Info
VDJdata = makeClassifier(VDJdata,NewHeader); %Classifier + FormattedSeq
VDJdata = findCDR3(VDJdata,NewHeader,varargin); %Get the CDR3 seq and info 
VDJdata = countSHM(VDJdata,NewHeader); %SHM info on the VMDNJ segments
VDJdata = labelNonprodVDJ(VDJdata,NewHeader,varargin); %Double check what is functional or not.
