%%%%%%%%%%%%%%%%%%%%% mStateEst.m %%%%%%%%%%%%%%%%%%%%%%%%
%
% This program is used for blind estimation of 
% m sequence initial state.
%
% date: 2025.3.25  GuRX
%
% [state] = mStateEst(ems, mpoly)
%
% ******************************************************
% ems   : estimated m sequence
% mpoly : estimated generator polynomial of m sequence
% state : estimated initial state of m sequence
% ******************************************************

function [state] = mStateEst(ems, mpoly)

N = length(ems);                % length of sequence
n = log2(N+1);                  % order of LFSR
L = n;                          % L=n : m sequence
                                % L=2n: Gold sequence
%%% method 1 : Walsh transform
X   = zeros(N-L, L);
c   = flip(mpoly(1:L));

for j = 1:L
    X(1, j) = c(L+1-j);
end

for i = 2:N-L
    if i <= L
        for j = 1:i-1
            temp = 0;
            for k = 1:i-1
                temp = temp + c(k) * X(i-k, j);
            end
            X(i, j) = mod(temp, 2);
        end
        for j = i:L
            temp = 0;
            for k = 1:i-1
                temp = temp + c(k) * X(i-k, j);
            end
            X(i, j) = mod(c(L+i-j)+temp, 2);
        end
    else
        for j = 1:L
            temp = 0;
            for k = 1:L
                temp = temp + c(k) * X(i-k, j);
            end
            X(i, j) = mod(temp, 2);
        end
    end
end

A = [X, flip(ems(L+1:N))];
V = zeros(N-L, 2^(L+1));
x = 2.^ (0:L).';
v = A * x;

for ii = 1 : N-L
    V(ii, v(ii)+1) = 1;
end

cv      = sum(V);                  % coefficient vector
H       = hadamard(2^(L+1));       % Hadamard matrix
cwh     = cv * H;
cwh1    = cwh(2:end);
[~,idx] = max(cwh1);
state   = dec2bin(idx) - '0';
state   = flip(state(2:end));

% figure;
% xx = 0:length(cwh)-1;
% stem(xx, cwh);
% title('Walsh Spectrum of state estimation');
% xlabel('Length of Walsh spectrum');
% ylabel('Amplitude');

% %%% method 2 : traverse state
% a = 2 * ems - 1;
% states = dec2bin(0:2^L-1) - '0';   % states
% vstates = states(2:end, :);
% 
% for i = 1:size(vstates, 1)
%     mSeqGen = comm.PNSequence( ...
%     'Polynomial', mpoly, ...                               
%     'InitialConditions', vstates(i, :), ...
%     'SamplesPerFrame', N);
%     b = 2*mSeqGen() - 1;
%     R = zeros(1, N);
%     for t = 1 : N
%         R(t) = 1/N * a.' * circshift(b, t-1);
%     end
%     [~, idx] = max(R);
%     if idx == 1
%         state = vstates(i, :);
%         break;
%     end
% end

% %%% method 3 : directly reconstruct m sequence
% a  = 2 * ems - 1;
% rs = [0 0 0 0 0 0 0 0 1];          % random state
% mSeqGen = comm.PNSequence( ...
%     'Polynomial', mpoly, ...
%     'InitialConditions', rs, ...
%     'SamplesPerFrame', N);
% b = 2*mSeqGen() - 1;
% R = zeros(1, N);
% for t = 1 : N
%     R(t) = 1/N * a.' * circshift(b, t-1);
% end
% [~, idx] = max(R);
% b        = circshift(b, idx-1);
% figure;
% x = 0 : length(R)-1;
% stem(x, R);
% title('Autocorrelation Function');
% xlabel('$\tau$', 'Interpreter', 'latex');
% ylabel('Amplitude of normalization')
% xlim([0 N-1]);

%*********************** end of file ************************