%generateVDJlib will generated a VDJ repertoire 
%
%  [VDJdata, NewHeader] = generateVDJrep(RepSize,GrpSize,SHMrate,TDTon)
%     RepSize = number of unique clones
%     GrpSize = number of Seq per clone
%     SHMrate = number of SHM mutation per 100 bp
%     TDTon = enable/disable TDT insertions. 1 = enable, 0 = disable.
%function [VDJdata, NewHeader] = generateVDJlib(RepSize,GrpSize,SHMrate,TDTon)

function [VDJdata, NewHeader] = generateVDJlib()

%Ask users to setup the parameters
Species = input('What is the species? 1 = mouse (Default), 2 = human');
if isempty(Species)
    Species = 1;
end

SetLen = input('What is the length of the nt sequence? Default 125');
if isempty(SetLen)
    SetLen = 125;
end

RepSize = input('How many in-frame CDR3 sequences do you want? Default 1000');
if isempty(RepSize)
    RepSize = 1000;
end

TDTon = input('Do you want TDT to be added? 1 = yes (default), 0 = no');
if isempty(TDTon)
    TDTon = 1;
end
if TDTon == 0
    InsCap = 0; %no tdt
else
    InsCap = 13; %Maximum N length
end

if Species == 1 %Mouse
    %These were collected from the C57BL6 mice data set from A Collins, 2015,
    %processed with MAVRIC Dev 9.
    MeanVdel = 1;
    MeanD5del = 4.5;
    MeanD3del = 3.4;
    MeanJdel = 3.9;
    MeanMlen = 3.8;
    MeanNlen = 2.9;
    DelCap = 13; %Maximum deletion
    
    %Make sure to choose the right database
    [Vmap,Dmap,Jmap] = getCurrentDatabase('change');
    [Vmap,Dmap,Jmap] = filterRefGene(Vmap,Dmap,Jmap);

elseif Species == 2 %Human
    %Souto-Carneiro, M.M., et al., Characterization of the Human Ig Heavy Chain Antigen Binding Complementarity Determining Region 3 Using a Newly Developed Software Algorithm, JOINSOLVER. The Journal of Immunology, 2004. 172(11): p. 6790-6802.
    MeanVdel = 2.1;
    MeanD5del = 4.5;
    MeanD3del = 5.0;
    MeanJdel = 6.8;
    MeanMlen = 7;
    MeanNlen = 7;
    DelCap = 13; %Maximum deletion

    %Make sure to choose the right database
    [Vmap,Dmap,Jmap] = getCurrentDatabase('change');
    [Vmap,Dmap,Jmap] = filterRefGene(Vmap,Dmap,Jmap,0);
end

%TDT insertion bias, based on mouse data. 
ACGTprob = [0.25 0.08 0.60 0.07]; %Prob of A, C, G, T, respectively.
FLIPprobVD = 0.25; %Probability of flipping Nvd side.
FLIPprobDJ = 0.45; %Probability of flipping Ndj side.

%Create the valid map idx 
Vidx = 1:size(Vmap,1);
for j = 1:length(Vidx)
    if isempty(Vmap{j,1})
        Vidx(j) = 0;
    end
end
Vidx(Vidx==0) = [];

%Create the valid map idx 
Didx = 1:size(Dmap,1);
for j = 1:length(Didx)
    if isempty(Dmap{j,1})
        Didx(j) = 0;
    end
end
Didx(Didx==0) = [];

%Create the valid map idx 
Jidx = 1:size(Jmap,1);
for j = 1:length(Jidx)
    if isempty(Jmap{j,1})
        Jidx(j) = 0;
    end
end
Jidx(Jidx==0) = [];

%Rearranging columns tsv files as csv files for MatLab
[~, ~, StandardData] = xlsread('Headers_BRILIA.xlsx');
NewHeaderLoc = findHeader(StandardData(1,:),'VDJdata');
NewHeader = StandardData(2:end,NewHeaderLoc);
for j = 1:length(NewHeader)
    if isnan(NewHeader{j}); break; end
end
NewHeader(j:end) = [];
NewHeader = NewHeader';

getHeaderVar;

VDJdata = cell(RepSize,length(NewHeader));

%--------------------------------------------------------------------------
q = 1;
while q <= RepSize
    %Select a V, D, J gene
    Vsel = ceil(rand(1)*length(Vidx));
    Dsel = ceil(rand(1)*length(Didx));
    Jsel = ceil(rand(1)*length(Jidx));
    Vnum = Vidx(Vsel);
    Dnum = Didx(Dsel);
    Jnum = Jidx(Jsel);

    %Determine NT deletions and N1/N2 lengths
    Vdel = round(-MeanVdel*log(rand(1)));
    D5del = round(-MeanD5del*log(rand(1)));
    D3del = round(-MeanD3del*log(rand(1)));
    Jdel = round(-MeanJdel*log(rand(1)));
    Mlen = round(-MeanMlen*log(rand(1)));
    Nlen = round(-MeanNlen*log(rand(1)));
    
    %Make sure to cap the deletion and insertion lengths
    if Vdel > DelCap;  Vdel = DelCap;  end
    if D5del > DelCap; D5del = DelCap; end
    if D3del > DelCap; D3del = DelCap; end
    if Jdel > DelCap;  Jdel = DelCap;  end
    if Mlen > InsCap;  Mlen = InsCap;  end
    if Nlen > InsCap;  Nlen = InsCap;  end

    %Makesure the deletion makes sense
    Wloc = Jmap{Jnum,end} + 2; %Want to full codon
    if Jdel >= Wloc; continue; end
    Jseq = Jmap{Jnum,1}(Jdel+1:Wloc);

    Dseq = Dmap{Dnum,1};
    if (D5del + D3del) >= length(Dseq); continue; end
    Dseq = Dseq(D5del+1:end-D3del);

    Mseq = generateNregion(Mlen,ACGTprob,FLIPprobVD);
    Nseq = generateNregion(Nlen,ACGTprob,FLIPprobDJ);

    Cloc = Vmap{Vnum,end};
    if Vdel >= Cloc; continue; end
    Vlen = SetLen - sum(length([Mseq Dseq Nseq Jseq]));
    Vseq = Vmap{Vnum,1}(end-Vdel-Vlen+1:end-Vdel);    
    VMDNJ = [length(Vseq) length(Mseq) length(Dseq) length(Nseq) length(Jseq)];
    
    %Assemble the NT sequence, in a similar fashion to the real 125 bp data
    NTseq = sprintf('%s%s%s%s%s',Vseq,Mseq,Dseq,Nseq,Jseq);
    
    %Check to see if there is a productive CDR3
    VallowedDel = Vmap{Vnum,end};
    CDR3nt = NTseq(VMDNJ(1)+Vdel-VallowedDel+1:end);
    if mod(length(CDR3nt),3) ~= 0
        continue
    end
    
    CDR3aa = nt2aa(CDR3nt,'ACGTOnly','false');
    if ~isempty(regexp(CDR3aa,'*','once'))
        continue
    end  

    %Save informations
    VDJdata(q,[SeqNumLoc GrpNumLoc SeqNameLoc]) = repmat({q},1,3);
    VDJdata{q,SeqLoc} = NTseq;
    VDJdata{q,FamNumLoc(1)} = Vnum;
    VDJdata{q,FamLoc(1)} = Vmap{Vnum,3};
    VDJdata{q,FamNumLoc(2)} = Dnum;
    VDJdata{q,FamLoc(2)} = Dmap{Dnum,3};
    VDJdata{q,FamNumLoc(3)} = Jnum;
    VDJdata{q,FamLoc(3)} = Jmap{Jnum,3};
    VDJdata(q,LengthLoc) = num2cell(VMDNJ);
    VDJdata(q,DelLoc) = num2cell([Vdel D5del D3del Jdel]);
    VDJdata(q,CDR3Loc) = {CDR3aa length(CDR3aa) VMDNJ(1)+Vdel-Cloc+1 sum(VMDNJ(1:4))-Jdel+Wloc};
    
    VDJdata(q,:) = labelNonprodVDJ(VDJdata(q,:),NewHeader);
    
    %Final check for nonprodVDJ
    if strcmpi(VDJdata(q,FunctLoc),'n')
        continue
    else
        q = q+1;
    end
end

%Fill in the details now
VDJdata = buildVDJalignment(VDJdata,NewHeader,Vmap,Dmap,Jmap); %Alignment Info
VDJdata = appendMutCt(VDJdata,NewHeader); %SHM infor on the VMDNJ segments
VDJdata = buildRefSeq(VDJdata,NewHeader,'germline','single'); %must do singles, since group therapy not done.
VDJdata = makeClassifier(VDJdata,NewHeader); %Classifier + FormattedSeq
VDJdata = findCDR3(VDJdata,NewHeader); %Get the CDR3 seq and info 

%Also save the parameters used
Settings = {'MeanVdel' MeanVdel; 'MeanD5del' MeanD5del; 'MeanD3del' MeanD3del; 'MeanJdel' MeanJdel; 'MeanMlen' MeanMlen; 'MeanNlen' MeanNlen; 'TDT_ACGTprob' mat2str(ACGTprob); 'TDT_FLIPprobVD' FLIPprobVD; 'TDT_FLIPprobDJ' FLIPprobDJ};

%Setup the file name
[SaveNamePre, SavePath] = uiputfile('*.*','Save Library as');
SaveName = sprintf('%s_RepSize%0.0d',SaveNamePre,RepSize);

%Restructure the alignment data
VDJdata = reformatAlignment(VDJdata,1);
if ispc
    xlswrite([SavePath SaveName '.xlsx'],cat(1,NewHeader,VDJdata));
    xlswrite([SavePath SaveName '.xlsx'],Settings,'Settings');    
else
    writeDlmFile(cat(1,NewHeader,VDJdata),[SavePath SaveName '.csv'],'\t');
    writeDlmFile(Settings,[SavePath SaveName 'Settings.csv'],'\t');
end