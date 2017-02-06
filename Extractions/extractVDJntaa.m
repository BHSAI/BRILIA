%extractVDJntaa will extract the full-length nt and aa of the VDJ segments,
%using the reference V and J segments + sample's variable NDN junction.
%This script it used mainly to for protein-digest evaluation of the VDJ
%segment.

function [AllRefSeq, AllRefAA] = extractVDJntaa()

%Select multiple files
[FileNames, FilePath] = uigetfile('*.xlsx; *.mat','Select the VDJdata files','multiselect','on');

%Convert to cells for consistency
if ischar(FileNames)
    FileNames = {FileNames};
end

AllRefSeq = {};
AllRefAA = {};
for j = 1:length(FileNames)
    FileName = FileNames{j};
    FullName = [FilePath FileName];
    [VDJdata, VDJheader, ~, ~] = openSeqData(FullName);
    
    %Determine the RefSeq
    RefSeq = buildRefSeq(VDJdata,VDJheader,'single','VJ','full');
    RefAA = translateNT(RefSeq);
    
    AllRefSeq = cat(1,AllRefSeq,RefSeq);
    AllRefAA = cat(1,AllRefAA,RefAA);
end
