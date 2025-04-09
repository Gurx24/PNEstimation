%%%%%%%%%%%%%%%%%%%%% GoldStateEst.m %%%%%%%%%%%%%%%%%%%%%%%%
%
% This program is used for blind estimation of 
% Gold sequence initial state.
%
% date: 2025.3.25  GuRX
%
% [state] = GoldStateEst(egs)
%
% ******************************************
% egs   : estimated Gold sequence
% state : initial state of Gold sequence
% ******************************************

function [state] = GoldStateEst(egs, gpoly)

N = length(egs);              % length of sequence
n = log2(N+1);                % order of LFSR
L = 2*n;                      % L=n : m sequence
                              % L=2n: Gold sequence
X = zeros(N-L, L);
c = flip(gpoly(1:L));

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
            X(i, j) = mod(c(L+i-j) + temp, 2);
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

A = [X, flip(egs(L+1:N))];
x = zeros(L+1, 1);
v = zeros(N-L, 1);
V = zeros(N-L, 2^(L+1));

for i = 1 : L+1
    x(i) = 2^(i-1);
end

for ii = 1 : N-L
    v(ii)          = A(ii, :) * x;
    V(ii, v(ii)+1) = 1;
end

c        = sum(V);       % coefficient vector
l        = 9;            % index of segment
ls       = 2^l;          % length of segment
V1       = reshape(c, [ls, 2^(L+1-l)]).';
H1       = hadamard(2^l);
cwh1     = sum(abs(V1 * H1));
cwh11    = cwh1(2:end);
[~,idx1] = max(cwh11);
C_rows   = dec2bin(idx1, l) - '0';
H2       = hadamard(2^(L+1-l));
cwh2     = sum(abs(V1.' * H2));
cwh22    = cwh2(2:end);
[~,idx2] = max(cwh22);
C_cols   = dec2bin(idx2, L+1-l) - '0';
state    = [C_cols, C_rows];
state    = flip(state(2:end));

% figure;
% subplot(1, 2, 1);
% stem(0:length(cwh1)-1, cwh1);
% title('Walsh Spectrum of Row');
% xlabel('Length of Walsh spectrum');
% ylabel('Amplitude');
% 
% subplot(1, 2, 2);
% stem(0:length(cwh2)-1, cwh2);
% title('Walsh Spectrum of Column');
% xlabel('Length of Walsh spectrum');
% ylabel('Amplitude');