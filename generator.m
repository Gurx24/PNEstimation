%%%%%%%%%%%%%%%%%%%%%%%%%% generator %%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function generates an m-sequence or a Gold sequence 
% based on the given polynomial and initial state.
%
% date : 2025.4.2  GuRX
%
% seq = generator(poly, inistate)
%
% ***********************************************
% poly     : generator polynomial
% inistate : initial state
% seq      : generated sequence
% ***********************************************

function seq = generator(poly, inistate)

lfsr = inistate;
L    = length(lfsr);
len  = 2^L - 1;
seq  = zeros(len, 1);

for i = 1:len
    seq(i)   = lfsr(end);
    feedback = mod(sum(poly(2:end).*lfsr), 2);
    lfsr     = [feedback, lfsr(1:end-1)];
end

%************************ end of file *************************