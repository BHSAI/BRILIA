%buildRefSeq will fill in the predicted germline reference sequence for
%each sequence cluster.
%
%  VDJdata = buildRefSeq(VDJdata, VDJheader)
%
%  [VDJdata, NewRange] = buildRefSeq(VDJdata, VDJheader)
%
%  VDJdata = buildRefSeq(VDJdata, VDJheader, LengthOpt, InnerOpt, GroupOpt)
%
%  INPUT OPTIONS, which can be in any order in inputs
%    LengthOpt: 
%      'full'      returns full length references
%      'same'      returns same-length references as input seq
%    InnerOpt: 
%      'germline'  returns same-length sequences, subsituting VDJ with
%                  germline genes. Nvd and Ndj are either from Seq1 or each
%                  seq, depending on what the GroupOpt is set at.
%      'preserve'  returns germline V + J ends with preserved, raw data
%                  sequence (useful for protein digest studies, required
%                  all sequence base the VDJ region).
%    GroupOpt: 
%      'single'    finds one RefSeq for each sequence
%      'group'     finds and applies one RefSeq for entire group
%      'first'     finds RefSeq for only the 1st seq in a group, in 
%                  case you used the buildTree first, which will fix RefSeq
%                  to ParentSeq.
%
%  OUTPUT
%    NewRange: the start and end nt sequence positions in Seq that the
%      reference sequence is extracted from.
function [VDJdata, varargout] = buildRefSeq(VDJdata,VDJheader,varargin)
%Extract the full VDJ database
[Vmap, Dmap, Jmap] = getCurrentDatabase; %Always use the full database for this
H = getHeaderVar(VDJheader);

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

%Change group number temporarily if you'll process all sequence
if strcmpi(GroupOpt,'single')
    GrpNum = [1:size(VDJdata,1)]';
    UnqGrpNum = GrpNum;
else    
    GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
    UnqGrpNum = unique(GrpNum);
end

NewRange = zeros(size(VDJdata,1),2);
for y = 1:length(UnqGrpNum)
    %Determine if group correction is possible
    IdxLoc = find(UnqGrpNum(y) == GrpNum);
    Tdata = VDJdata(IdxLoc,:);
    
    try
        %Extract the other informations about ref sequences
        DelCt = unique(cell2mat(VDJdata(IdxLoc,H.DelLoc)),'rows');    
        VmapNum = VDJdata{IdxLoc(1),H.FamNumLoc(1)}(1);
        DmapNum = VDJdata{IdxLoc(1),H.FamNumLoc(2)}(1);
        JmapNum = VDJdata{IdxLoc(1),H.FamNumLoc(3)}(1);

        %Assemble the full length germline VDJ first
        Seq = Tdata{1,H.SeqLoc};
        VMDNJ = unique(cell2mat(Tdata(:,H.LengthLoc)),'rows');
        if size(VMDNJ,1) > 1
            error('This data has not be conformed to unity per group. Perform group treatment first');
        end
        
        %Determine VMDNJ reference sequence, if possible
        VrefSeq = Vmap{VmapNum(1),1}(1:end-DelCt(1));
        DrefSeq = Dmap{DmapNum(1),1}(1+DelCt(2):end-DelCt(3));
        JrefSeq = Jmap{JmapNum(1),1}(1+DelCt(end):end);
        if VMDNJ(2) ~= 0
            Mseq = Seq(VMDNJ(1)+1:sum(VMDNJ(1:2)));
        else 
            Mseq = '';
        end
        if VMDNJ(4) ~= 0
            Nseq = Seq(sum(VMDNJ(1:3))+1:sum(VMDNJ(1:4)));
        else 
            Nseq = '';
        end
        
        %Assemble the full length sequence, left side
        ExtraLeft = VMDNJ(1) - length(VrefSeq);
        if ExtraLeft > 0
            VrefSeq = [repmat('X',1,ExtraLeft) VrefSeq];
        elseif ExtraLeft < 0
            VrefSeq(1:abs(ExtraLeft)) = [];
        end
        
        %Assemble the full length sequence, right side
        ExtraRight = VMDNJ(5) - length(JrefSeq);
        if ExtraRight > 0 
            JrefSeq = [JrefSeq repmat('X',1,ExtraRight)];
        elseif ExtraRight < 0
            JrefSeq(end-(abs(ExtraRight))+1:end) = [];
        end
        
        %Assemble the full sequence and identify start to end locations
        FullRefSeq = sprintf('%s%s%s%s%s',VrefSeq,Mseq,DrefSeq,Nseq,JrefSeq);
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
                VDJdata(IdxLoc(1),H.RefSeqLoc) = {FullRefSeq};
                VDJdata(IdxLoc(2:end),H.RefSeqLoc) = Tdata(2:end,H.RefSeqLoc);
            case 'group'
                VDJdata(IdxLoc,H.RefSeqLoc) = repmat({FullRefSeq},length(IdxLoc),1);
            case 'single'
                VDJdata(IdxLoc(1),H.RefSeqLoc) = {FullRefSeq};
        end
        NewRange(IdxLoc,1) = S1;
        NewRange(IdxLoc,2) = S2;
    catch
        ErrorMsg = sprintf('Errored at %s, sequence # %d',mfilename,y);
        disp(ErrorMsg);
        VDJdata(IdxLoc,H.MiscLoc) = repmat({ErrorMsg},length(IdxLoc),1);
    end
end

if nargout >= 2
    varargout{1} = NewRange;
end
