%updateVDJdata will perform simple updates to VDJ data, EXCEPT for building
%the reference sequences as this requires some user input.
%
%  VDJdata = updateVDJdata(VDJdata,VDJheader)
%
%  VDJdata = updateVDJdata(VDJdata,VDJheader,Vmap,Dmap,Jmap)

function VDJdata = updateVDJdata(VDJdata,VDJheader,varargin)

VDJdata = buildVDJalignment(VDJdata,VDJheader,varargin); %Alignment Info
VDJdata = makeClassifier(VDJdata,VDJheader); %Classifier + FormattedSeq
VDJdata = findCDR3(VDJdata,VDJheader,varargin); %Get the CDR3 seq and info 
VDJdata = countSHM(VDJdata,VDJheader); %SHM info on the VMDNJ segments
VDJdata = labelNonprodVDJ(VDJdata,VDJheader,varargin); %Double check what is functional or not.
