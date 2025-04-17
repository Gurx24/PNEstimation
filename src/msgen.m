%%%%%%%%%%%%%%%%%%%%%%%%%% msgen.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function is used to generate m sequence.
% Multiple m sequences can be generated based on 
% multiple generator polynomials.
%
% date : 2025.4.2  GuRX
%
% mseq = msgen(pair, L)
%
% *********************************************************
% polys : octal generator polynomial vector of m sequence
% L     : order of m sequence
% mseq  : generated m sequence
% *********************************************************

function mseq = msgen(polys, L, inist)

switch nargin
    case 2
        state = [zeros(1, L-1) 1];      % default initial state
    case 3
        state = inist;                  % specified initial state
end

N     = 2^L - 1;                      % length of m sequence
mseq  = zeros(size(polys,2), N);      % initialization of m sequence
bpair = dec2bin(base2dec(num2str(polys.'), 8))-'0';

for i = 1:size(bpair, 1)
    mSeqGen = comm.PNSequence( ...
        'Polynomial', bpair(i, :), ...
        'InitialConditions', state, ...
        'SamplesPerFrame', N);
    mseq(i, :) = mSeqGen();
end