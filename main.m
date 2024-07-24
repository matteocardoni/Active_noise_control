clear all

%% Initial settings
M = 64;    % subbands
D = 16;    % decimation factor
[signal, fs] = audioread('rand.wav');
% Convert signal into a row vector
signal = signal';
buf_shift = 100;
x_len = 2*M;

l_SAF = 64; % Length of each subband adaptive filter
N = l_SAF * M;  % Length of the adaptive filter W
mu = 0.0001;
epsilon = 0.0001;
% Read filters values form .dat file
p_file = fopen('p.dat', 'r');
p = fread(p_file);  % Primary path
s_file = fopen('s.dat', 'r');
s = fread(s_file);  % Secondary path
% Adaptive filter w:
w = zeros(M, 1);

%% Creation of the DFT matrix:
F = zeros(M,M);
for m = 1:M
    for n = 1:M
        F(m,n) = exp(-1j*2*pi*(m-1)*(n-1)/M); % there is k-1 and i-1
        % because matlab starts from 0
    end
end

%% Create the decimated DFT matrix
DF = zeros(M, M/D);
for i = 1:M
    % Implementation of decimation:
    % h_0 = h(M*n), h_1 = h(M*n + 1), ..., h_(M-1) = h(M*n + M - 1)
    % (but in matlab we have to start from 1)
    test = F(i, i:D:length(F(i, :)));
    % Start point from the decimatin, comes back to 1 when it reaches D.
    decimation_start = i - ceil((i - D)/D)*D ;
%     % In order to start from numbers greater than D and than come back to
%     % the first samples, circular shift (left) is implemented after the
%     % decimation
%     DF(i, :) = circshift(F(i, decimation_start:D:length(F(i, :))), - ceil((i - D)/D));
      DF(i, :) = F(i, decimation_start:D:length(F(i, :)));
end

%% Buffer
buf = zeros(1, x_len + D);  % Buffer to store M samples in
iter = 1;   % iteration

%% Store a number of samples
% Store M samples in buf, until there are M still available
buf(1, 1:x_len) = signal(1, 1:x_len);
% Instance the matrix of for the w_k filters (in time)
w_k = zeros(M, x_len);
% Instance the matrix of for the W vectors, to stack to obtain the
% w filter (in frequency):
W_k = zeros(M, N);
% Iterate until the sound signal ends:
for iter = 2:floor(length(signal)/M)
    x = buf(1, 1:x_len);
    % Filter
    d = filter(p, 1, x);
    y = filter(w, 1, x);
    y_filtered = filter(s, 1, y);
    x_estimated = filter(s, 1, x);
    e = d - y_filtered;
    % Instance the matrix of x_k (in time)
    x_k = zeros(M, x_len);
    % Instance the matrix of e_k (in time)
    e_k = zeros(M, x_len);
    % Filter x (already filtered for S') for the UDFT filter bank
    for i = 1:M
        x_k(i, :) = filter(DF(i, :), 1, x_estimated);
    end
    % Filter e for the UDFT filter bank
    for i = 1:M
        e_k(i, :) = filter(DF(i, :), 1, e);
    end
    % Update the w^SAF_k coefficients 
    for i = 1:M
        % 1) IN THIS FORMULA I HAVE USED THE ELEMENT-WISE MULTIPLICATION
        %       BETWEEN x AND e
        % 2) I DON'T HAVE TO FILTER x_k AND e_k WITH w^SAF_k, RIGHT?
        w_k(i, :) = w_k(i, :) + mu * conj(x_k(i, :)).*e_k(i, :)/(epsilon + norm(x_k(i, :))^2);
        % Perform DFT at N points of the w^SAF_k filters
        W_k(i, :) = fft(w_k(i,:), N);
    end
    % Do the weights stacking
    
    % Shift the signal in the buffer
    buf = circshift(buf, buf_shift);
    % Add new samples from signal (and check if the signal length is not
    % exceeded.
    if x_len + buf_shift*iter > length(signal)
        break
    else
        buf(1, 1:buf_shift) = signal(1 ,x_len + buf_shift*(iter-1)+1: x_len + buf_shift*iter);
    end
end