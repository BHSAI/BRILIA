%updateVDJdata will perform simple updates to VDJ data EXCEPT for building
%the reference sequences as this requires some user input.
%
%  VDJdata = updateVDJdata(VDJdata,VDJheader)
%
%  VDJdata = updateVDJdata(VDJdata,VDJheader,DB)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%    DB: Gene database structure (getCurrentDatabase.m)
%
%  OUTPUT
%    VDJdata: modified VDJdata germline-child SHM counts, CDR3 start and
%      end anchors, and sequence functionality data.

function VDJdata = updateVDJdata(VDJdata, Map, DB)
VDJdata = countSHM(VDJdata, Map); %SHM info on the VMDNJ segments
VDJdata = findCDR(VDJdata, Map, DB, 3); %Get the CDR3 seq and info, using anchor method

