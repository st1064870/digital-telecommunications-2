function [smt] = M_PSK_Transmitter(inpt, M, encoding)

R = log2(M);

% We have 4 samples (T_sample) per carrier period (Tc) and
% 10 carrier periods (Tc) per symbol period (T_sym).

% T_scale = 4 * 10^-6;

T_sample = 1;
T_c = 4;
f_c = 1/T_c;
T_sym = 40;

g_t = sqrt(2 / T_sym);

% Constellation
constellation = zeros(M, 2);
for m = 0:M-1
    constellation(m+1,:) = [cos(2 * pi * m/M) sin(2 * pi * m/M)];
end

% -----------------------------Mapper--------------------------------------

blocks_inpt = reshape(inpt.', R, []).';
if strcmp(encoding, 'bin')
    sm = bi2de(blocks_inpt, 'left-msb');
else
    sm = bi2de(blocks_inpt, 'left-msb');
    sm = bin2gray(sm, 'psk', M);
end

tmp = repelem(sm, T_sym);
figure('Position', [10 10 900 600])
subplot(4,1,1);
plot(1:T_sample:T_sample*T_sym*10, tmp(1:T_sample*T_sym*10));
title("Subplot 1: Original signal")


% ----------------------------Modulation-----------------------------------

smt = zeros(length(sm) * T_sym, 1);
t = (0:T_sample:(T_sym - T_sample))';

for i = 1:length(sm)
    for j = 1:length(t)
        smt((i - 1) * T_sym + j) = ...
            g_t * cos(2 * pi * (sm(i)/M)) * cos(2 * pi * f_c * t(j)) + ...
            g_t * sin(2 * pi * (sm(i)/M)) * sin(2 * pi * f_c * t(j));
    end
end

subplot(4,1,2);
plot(1:T_sample:T_sample*T_sym*10, smt(1:T_sample*T_sym*10));
title("Subplot 2: Transmitter's output");


end

