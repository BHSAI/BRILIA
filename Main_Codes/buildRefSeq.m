%buildRefSeq will take VDJdata and rebuild the reference sequence, but with
%several varations in methods. 
%
%  RefSeq = buildRefSeq(VDJdata,NewHeader,LengthOpt,OptionA,OptionB)
%  
%    LengthOpt = 'full'   will return full length references
%    LengthOpt = 'same'   will return same-length references as input seq
%
%    InnerOpt = 'germline' will return same-length sequences, subsituting
%      VDJ with germline genes. N1 and N2 are either from Seq1 or each seq,
%      depending on what the GroupOpt is set at.
%    InnerOpt = 'preserve' = wil match germline V + J ends with preserved,
%    raw data sequence (useful for protein digest studies, required all
%    sequence base the VDJ region).
%
%    GroupOpt = 'single' Find RefSeq for each sequence
%    GroupOpt = 'group'  Find RefSeq for entire group
%    GroupOpt = 'first'  Find RefSeq for just 1st seq in a group (in case
%    you used the buildTree first, which will fix RefSeq to ParentSeq).

function [VDJdata, varargout] = buildRefSeq(VDJdata,NewHeader,varargin)
%Extract the full VDJ database
[Vmap, Dmap, Jmap] = getCurrentDatabase; %Always use the full database for this
getHeaderVar;

%Setup the options
LengthOpt = 'same';
InnerOpt = 'germline';
GroupOpt  = 'first';
for k = 1:length(varargin)
    switch lower(varargin{k})
        case 'full'
            LengthOpt = 'full';
        case 'same'
            LengthOpt = 'same';
        case 'germline'
            InnerOpt = 'germline';
        case 'preserve'
            InnerOpt = 'preserve';
        case 'single'
            GroupOpt = 'single';
        case 'group'
            GroupOpt = 'group';
        case 'first'
            GroupOpt = 'first';
    end
end

%Change group number depending on if you'll process all sequence, or just 1
%per group.
if strcmpi(GroupOpt,'single')
    GrpNum = [1:size(VDJdata,1)]';
    UnqGrpNum = GrpNum;
else    
    GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
    UnqGrpNum = unique(GrpNum);
end

NewRange = zeros(size(VDJdata,1),2);
for y = 1:length(UnqGrpNum)
    %Determine if group correction is possible
    IdxLoc = find(UnqGrpNum(y) == GrpNum);
    Tdata = VDJdata(IdxLoc,:);

    %Extract the other informations about ref sequences
    DelCt = unique(cell2mat(VDJdata(IdxLoc,DelLoc)),'rows');    
    VmapNum = VDJdata{IdxLoc(1),FamNumLoc(1)}(1);
    DmapNum = VDJdata{IdxLoc(1),FamNumLoc(2)}(1);
    JmapNum = VDJdata{IdxLoc(1),FamNumLoc(3)}(1);
    
    %Assemble the full length germline VDJ first
    Seq = Tdata{1,SeqLoc};
    VMDNJ = unique(cell2mat(Tdata(:,LengthLoc)),'rows');
    if size(VMDNJ,1) > 1
        error('This data has not be conformed to unity per group. Perform group treatment first');
    end
       
    VrefSeq = Vmap{VmapNum(1),1}(1:end-DelCt(1));
    DrefSeq = Dmap{DmapNum(1),1}(1+DelCt(2):end-DelCt(3));
    JrefSeq = Jmap{JmapNum(1),1}(1+DelCt(end):end);
    if VMDNJ(2) ~= 0
        N2Seq = Seq(VMDNJ(1)+1:sum(VMDNJ(1:2)));
    else 
        N2Seq = '';
    end
    if VMDNJ(4) ~= 0
        N1Seq = Seq(sum(VMDNJ(1:3))+1:sum(VMDNJ(1:4)));
    else 
        N1Seq = '';
    end
    FullRefSeq = sprintf('%s%s%s%s%s',VrefSeq,N2Seq,DrefSeq,N1Seq,JrefSeq);
    S1 = length(VrefSeq) - VMDNJ(1) + 1; 
    S2 = S1 + length(Seq) - 1;
    
    %Switch the inner components if needed
    if strcmpi(InnerOpt,'preserve')
        FullRefSeq(S1:S2) = Seq;
    end
    
    %Trim the seq if needed
    if strcmpi(LengthOpt,'same')
        FullRefSeq = FullRefSeq(S1:S2);
        S1 = 1;
        S2 = length(FullRefSeq);
    end
    
    %Fill in RefSeq as needed
    switch GroupOpt
        case 'first'
            VDJdata(IdxLoc(1),RefSeqLoc) = {FullRefSeq};
            VDJdata(IdxLoc(2:end),RefSeqLoc) = Tdata(2:end,RefSeqLoc);
        case 'group'
            VDJdata(IdxLoc,RefSeqLoc) = repmat({FullRefSeq},length(IdxLoc),1);
        case 'single'
            VDJdata(IdxLoc(1),RefSeqLoc) = {FullRefSeq};
    end
    NewRange(IdxLoc,1) = S1;
    NewRange(IdxLoc,2) = S2;
end

if nargout >= 2
    varargout{1} = NewRange;
end