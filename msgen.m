%%%%%%%%%%%%%%%%%%%%%%%%%% msgen.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function is used to generate m sequence.
%
% date : 2025.4.2  GuRX
%
% mseq = msgen(pair, L)
%
% ***********************************************
% pair : generator polynomial pair of m sequence
% L    : order of m sequence
% mseq : generated m sequence
% ***********************************************

function mseq = msgen(pair, L)

rs    = [zeros(1, L-1) 1];            % initial state
N     = 2^L - 1;                      % length of m sequence
mseq  = zeros(size(pair,2), N);       % initialization of m sequence
bpair = dec2bin(base2dec(num2str(pair.'), 8))-'0';

for i = 1:size(bpair, 1)
    mSeqGen = comm.PNSequence( ...
        'Polynomial', bpair(i, :), ...
        'InitialConditions', rs, ...
        'SamplesPerFrame', N);
    mseq(i, :) = mSeqGen();
end