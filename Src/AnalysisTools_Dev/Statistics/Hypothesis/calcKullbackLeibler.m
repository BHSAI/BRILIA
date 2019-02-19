%calcKullbackLeibler will calculate the Kullback-Leibler Divergence
%for comparing two discrete propbability distibutions.
%
%  [Dpq, Dqp] = calcKullbackLeibler(P, Q)
%
%  INPUTS
%    P: Mx1 matrix of frequencies of dataset 1
%    Q: Mx1 matrix of frequencies of dataset 2
%
%  OUTPUTS
%    Dpq : D(p,q) divergence for predicting Q from P 
%    Dqp : D(q,p) divergence for predicting P from Q
%
%  NOTE
%    Uses the natural log.
%
%  Dpq = sum(P.*ln(P./Q));
%  Dpq = sum(Q.*ln(Q./P));

function [Dpq, Dqp] = calcKullbackLeibler(P, Q)

%Normalize distribution first
P = P/sum(P);
Q = Q/sum(Q);

Dpq = sum(P.*log(P./Q));
Dqp = sum(Q.*log(Q./P));
