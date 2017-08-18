%generateNregion will use create an N region segment based on the
%composition probability of nucleotides in N, and complement addition
%probability of N.
%
%  Nregion = generationNregion(Nlength)
%  Nregion = generationNregion(Nlength,ACGTprob,FLIPprob)
%
%  Default ACGTprob = [0.25 0.08 0.60 0.07] for A,C,G,T
%  Default FLIPprob = 0.33, meaning 33% chance to take complement seq. 

function Nregion = generateNregion(Nlength,varargin)
%Just in case Nlength is 0, quick return.
if Nlength == 0
    Nregion = '';
    return
end

%Prepare input for the TDTmatrix
if isempty(varargin)
    ACGTprob = [0.25 0.08 0.60 0.07]; %Prob of A, C, G, T, respectively.
    FLIPprob = 0.33; %Probability the N region adds to complment strand.
elseif length(varargin) == 2
    ACGTprob = varargin{1};
    FLIPprob = varargin{2};
end

%Generate the N region
RandN = rand(1,Nlength);
NTint = zeros(1,Nlength);
Pcomp = cumsum(ACGTprob);

%Assign the first NT
for j = 1:Nlength
    NTloc = find(RandN(j) <= Pcomp);
    NTint(j) = NTloc(1);
end
Nregion = int2nt(NTint);

if rand(1) <= FLIPprob;
    Nregion = seqcomplement(Nregion);
end
