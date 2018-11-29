%buildRefSeq will fill in the predicted germline reference sequence for
%each sequence cluster.
%
%  [VDJdata, NewRange] = buildRefSeq(VDJdata, VDJheader, DB)
%
%  VDJdata = buildRefSeq(VDJdata, VDJheader, DB, Chain, LengthOpt, InnerOpt, GroupOpt)
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
function [VDJdata, varargout] = buildRefSeq(VDJdata, Map, DB, varargin)
%Setup the options
LengthOpt = 'same';
InnerOpt  = 'germline';
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
        case 'h'
            Map.Chain = 'H';
        case 'l'
            Map.Chain = 'L';
        case 'hl'
            Map.Chain = 'HL';
    end
end

% Extract the full VDJ database
for k = 1:length(Map.Chain)
    if Map.Chain(k) == 'H'
        GrpNumLoc  = Map.GrpNum;
        DelLoc  = Map.hDel;
        GeneNumLoc  = Map.hGeneNum; 
        SeqLoc  = Map.hSeq;
        RefSeqLoc  = Map.hRefSeq;
        LengthLoc  = Map.hLength;
        GeneNameLoc  = Map.hGeneName;
    else
        GrpNumLoc  = Map.GrpNum;
        DelLoc  = Map.lDel;
        GeneNumLoc  = Map.lGeneNum; 
        SeqLoc  = Map.lSeq;
        RefSeqLoc  = Map.lRefSeq;
        LengthLoc  = Map.lLength;
        GeneNameLoc  = Map.lGeneName;
    end

    %Change group number temporarily if you'll process all sequence
    if strcmpi(GroupOpt, 'single')
        GrpNum = [1:size(VDJdata, 1)]';
        UnqGrpNum = GrpNum;
    else    
        GrpNum = cell2mat(VDJdata(:, GrpNumLoc));
        UnqGrpNum = unique(GrpNum);
    end

    NewRangeH = zeros(size(VDJdata, 1), 2);
    NewRangeL = zeros(size(VDJdata, 1), 2);
    for y = 1:length(UnqGrpNum)
        %Determine if group correction is possible
        IdxLoc = find(UnqGrpNum(y) == GrpNum);
        Tdata = VDJdata(IdxLoc, :);

        if Map.Chain(k) == 'H' %Heavy chain only
            %Extract the other informations about ref sequences
            DelCt = unique(cell2mat(VDJdata(IdxLoc, DelLoc)), 'rows');    
            Vnum = VDJdata{IdxLoc(1), GeneNumLoc(1)}(1);
            Dnum = VDJdata{IdxLoc(1), GeneNumLoc(2)}(1);
            Jnum = VDJdata{IdxLoc(1), GeneNumLoc(3)}(1);

            %Make sure nothing is empty
            if isempty(DelCt); continue; end
            if isempty(Vnum); continue; end
            if isempty(Dnum); continue; end
            if isempty(Jnum); continue; end
            Vnum = Vnum(1);
            Dnum = Dnum(1);
            Jnum = Jnum(1);
            
            %Assemble the full length germline VDJ first
            Seq = Tdata{1, SeqLoc};
            VMDNJ = unique(cell2mat(Tdata(:, LengthLoc)), 'rows');
            if size(VMDNJ, 1) > 1
                warning('%s: Group #%d annotation is not uniform.', mfilename, UnqGrpNum(y));
                continue
            end

            %Determine VMDNJ reference sequence, if possible
            try
                VrefSeq = DB.Vmap{Vnum, 1}(1:end-DelCt(1));
                DrefSeq = DB.Dmap{Dnum, 1}(1+DelCt(2):end-DelCt(3));
                JrefSeq = DB.Jmap{Jnum, 1}(1+DelCt(end):end);
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
            catch ME
                disp(ME)
                error('%s: failed to build germline seq.', mfilename)
            end

            %Assemble the full length sequence, left side
            ExtraLeft = VMDNJ(1) - length(VrefSeq);
            if ExtraLeft > 0
                VrefSeq = sprintf('%s%s', repmat('N', 1, ExtraLeft), VrefSeq);
            elseif ExtraLeft < 0
                VrefSeq(1:abs(ExtraLeft)) = [];
            end

            %Assemble the full length sequence, right side
            ExtraRight = VMDNJ(5) - length(JrefSeq);
            if ExtraRight > 0 
                JrefSeq = sprintf('%s%s', JrefSeq, repmat('N', 1, ExtraRight));
            elseif ExtraRight < 0
                JrefSeq(end-(abs(ExtraRight))+1:end) = [];
            end

            %Assemble the full sequence and identify start to end locations
            FullRefSeq = sprintf('%s%s%s%s%s', VrefSeq, Mseq, DrefSeq, Nseq, JrefSeq);
            S1 = length(VrefSeq) - VMDNJ(1) + 1; 
            S2 = S1 + length(Seq) - 1;

            %Switch the inner components if needed
            if strcmpi(InnerOpt, 'preserve')
                FullRefSeq(S1:S2) = Seq;
            end

            %Trim the seq if needed
            if strcmpi(LengthOpt, 'same')
                FullRefSeq = FullRefSeq(S1:S2);
                S1 = 1;
                S2 = length(FullRefSeq);
            end

            %Fill in RefSeq as needed
            switch GroupOpt
                case 'first'
                    VDJdata(IdxLoc(1), RefSeqLoc) = {FullRefSeq};
                    VDJdata(IdxLoc(2:end), RefSeqLoc) = Tdata(2:end, RefSeqLoc);
                case 'group'
                    VDJdata(IdxLoc, RefSeqLoc) = repmat({FullRefSeq}, length(IdxLoc), 1);
                case 'single'
                    VDJdata(IdxLoc(1), RefSeqLoc) = {FullRefSeq};
            end
            NewRangeH(IdxLoc, 1) = S1;
            NewRangeH(IdxLoc, 2) = S2;
        end

        if Map.Chain(k) == 'L' %Light chain only
            %Extract the other informations about ref sequences
            DelCt = unique(cell2mat(VDJdata(IdxLoc, DelLoc)), 'rows');    
            Vname = VDJdata{IdxLoc(1), GeneNameLoc};
            Vnum = VDJdata{IdxLoc(1), GeneNumLoc(1)};
            Jnum = VDJdata{IdxLoc(1), GeneNumLoc(end)};

            %Make sure nothing is empty
            if isempty(DelCt); continue; end
            if isempty(Vname); continue; end
            if isempty(Vnum); continue; end
            if isempty(Jnum); continue; end
            Vnum = Vnum(1);
            Jnum = Jnum(1);

            %Assemble the full length germline VDJ first
            Seq = Tdata{1, SeqLoc};
            VNJ = unique(cell2mat(Tdata(:, LengthLoc)), 'rows');
            if size(VNJ, 1) > 1
                warning('buildRefSeq: This data has not be conformed to unity per group. Perform group treatment first');
                continue;
            end

            %Determine which DB to use
            StrPat = 'IG[LK]';
            LocusLoc = regexp(Vname, StrPat, 'once');
            Locus = Vname(LocusLoc+2);

            %Determine VNJ reference sequence, if possible
            if Locus == 'K'
                VrefSeq = DB.Vkmap{Vnum, 1}(1:end-DelCt(1));
                JrefSeq = DB.Jkmap{Jnum, 1}(1+DelCt(end):end);
            else
                VrefSeq = DB.Vlmap{Vnum, 1}(1:end-DelCt(1));
                JrefSeq = DB.Jlmap{Jnum, 1}(1+DelCt(end):end);
            end
            if VNJ(2) ~= 0
                Nseq = Seq(VNJ(1)+1:sum(VNJ(1:2)));
            else 
                Nseq = '';
            end

            %Assemble the full length sequence, left side
            ExtraLeft = VNJ(1) - length(VrefSeq);
            if ExtraLeft > 0
                VrefSeq = sprintf('%s%s', repmat('N', 1, ExtraLeft), VrefSeq);
            elseif ExtraLeft < 0
                VrefSeq(1:abs(ExtraLeft)) = [];
            end

            %Assemble the full length sequence, right side
            ExtraRight = VNJ(3) - length(JrefSeq);
            if ExtraRight > 0 
                JrefSeq = sprintf('%s%s', JrefSeq, repmat('N', 1, ExtraRight));
            elseif ExtraRight < 0
                JrefSeq(end-(abs(ExtraRight))+1:end) = [];
            end

            %Assemble the full sequence and identify start to end locations
            FullRefSeq = sprintf('%s%s%s', VrefSeq, Nseq, JrefSeq);
            S1 = length(VrefSeq) - VNJ(1) + 1; 
            S2 = S1 + length(Seq) - 1;

            %Switch the inner components if needed
            if strcmpi(InnerOpt, 'preserve')
                FullRefSeq(S1:S2) = Seq;
            end

            %Trim the seq if needed
            if strcmpi(LengthOpt, 'same')
                FullRefSeq = FullRefSeq(S1:S2);
                S1 = 1;
                S2 = length(FullRefSeq);
            end

            %Fill in RefSeq as needed
            switch GroupOpt
                case 'first'
                    VDJdata(IdxLoc(1), RefSeqLoc) = {FullRefSeq};
                    VDJdata(IdxLoc(2:end), RefSeqLoc) = Tdata(2:end, RefSeqLoc);
                case 'group'
                    VDJdata(IdxLoc, RefSeqLoc) = repmat({FullRefSeq}, length(IdxLoc), 1);
                case 'single'
                    VDJdata(IdxLoc(1), RefSeqLoc) = {FullRefSeq};
            end
            NewRangeL(IdxLoc, 1) = S1;
            NewRangeL(IdxLoc, 2) = S2;
        end    
    end
end

if nargout >= 2
    varargout{1} = NewRangeH;
    if nargout >= 3
        varargout{2} = NewRangeL;
    end
end
