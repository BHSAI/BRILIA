%loadRelMutTable will return the relative mutation matrix from nt X (row)
%to nt Y (col).
%
%REF - Mouse L Chain) A Cui, et al., A Model of Somatic Hypermutation
%Targeting in Mice Based on High-Throughput Ig Sequencing Data. Journal of
%Immunology. 2016. 197: 3566-3547.
%
%REF - Mouse H Chain) D Lee, et al., BRILIA: Integrated Tool for High-
%Throughput Annotation and Lineage Tree Assembly of B-Cell Repertoires.
%Frontiers in Immunology. 2017.
%
%  relMutMat = loadRelMutTable();
%
%  INPUT
%    Species ['mouse' or 'human']: species
%    Chain ['H' or 'L']: Ig chain
%
%  OUTPUT
%    relMutMat: 4x4 matrix showing the relative probability of mutating nt
%      X (row) to nt Y (col). sum(MutMat, 2) = ones(4,1)

function MutMat = loadRelMutTable(Species, Chain)
if strcmpi(Species, 'mouse')
    if strcmpi(Chain, 'L')
        %Mouse L Chain, A CUI
        MutMat = [0     0.16  0.54  0.30  ;
                  0.13  0     0.09  0.78  ;
                  0.77  0.14  0     0.09  ;
                  0.25  0.57  0.18  0     ];

    elseif strcmpi(Chain, 'H')
        %Mouse H Chain, D LEE
        % Original values caluclation
        % MutMat = [0     9288  30650 12086; 
        %           7320  0     7874  23588;
        %           22093 8200  0     5886 ;
        %           5798  19160 7492  0    ;]
        % MutMat = MutMat ./ repmat(sum(MutMat, 2), 1, 4);
        % Precalculated for speed
        MutMat = [0         0.1785    0.5892    0.2323;
                  0.1887    0         0.2030    0.6083;
                  0.6107    0.2267    0         0.1626;
                  0.1787    0.5904    0.2309    0     ];
    end
elseif strcmpi(Species, 'human')
    if strcmpi(Chain, 'L')
        %Human L Chain, A Cui
        MutMat = [0     0.26  0.50  0.24  ;
                  0.25  0     0.30  0.45  ;
                  0.53  0.31  0     0.16  ;
                  0.20  0.57  0.23  0     ];

    elseif strcmpi(Chain, 'H')
        %Human H Chain, A Cui
        MutMat = [0     0.28  0.49  0.23  ;
                  0.21  0     0.35  0.44  ;
                  0.49  0.35  0     0.16  ;
                  0.19  0.44  0.37  0     ];
    end
end
