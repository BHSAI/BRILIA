%loadCodonSynMap will return a map storing a 4x3 codon synonymous mutations
%map. This codonSRMap stores for each codon, a 4x3 logical array for when
%the nth (column) nt in a codon mutations to A, C, G, or T (row) and leads
%to no mutations (1) or mutations (0);
%
%  codonSRMap = loadCodonSynMap();
%
%  OUTPUT
%    codonSRMap: map object storing the codon as keys, and a 4x3 logical
%      array of silent/replacement mutation to the codon. 1 means
%      mutation leads to replacements.
%
%  EXAMPLE
%    codonSRMap = loadCodonSynMap;
%    codonSRMap('ACT') =
%      0     1     0   %ACT -> AAT changes AA
%      1     0     0   %ACT -> CCT changes AA
%      1     1     0   %ACT -> GGT changes AA
%      1     1     0   %ACT -> TTT changes AA
function codonSRMap = loadCodonSynMap()

%Generate 64 codons and it's silent/replacement mutation matrix SRM 
IntCodonList = cell(64, 1);
j = 1;
for x = 1:4
    for y = 1:4
        for z = 1:4
            IntCodonList{j} = [x y z];
            j = j + 1;
        end
    end
end

%Store the key-value pair in a java map object
codonSRMap = containers.Map; %Java object so using Java naming convention lowercase 
for k = 1:length(IntCodonList)
    IntCodon = IntCodonList{k};
    IntAA = nt2aa(IntCodon, 'alternative', false);
    MutAA = zeros(4, 3); %Stores the new AA
    for n = 1:3
        for x = 1:4
            NewIntCodon = IntCodon;
            NewIntCodon(n) = x;
            MutAA(x, n) = nt2aa(NewIntCodon, 'alternative', false);
        end
    end
    
    %Add to codon map
    Codon = int2nt(IntCodon);
    codonSRMap(Codon) = ~(MutAA == IntAA);
end
