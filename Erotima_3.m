L_b = randsrc(24^4, 1, [0 1]);

smt_4 = M_PSK_Transmitter(L_b, 4, 'bin');
smt_8 = M_PSK_Transmitter(L_b, 8, 'bin');

T_sample = 0.1 * 10^-6;
f_sample = 1/T_sample;
n = 2048;
f = (0:n-1)*(f_sample/n);

% figure;
% 
% x_t = smt_4;
% [Pxx,F] = periodogram(x_t,[],length(x_t),f_sample);
% plot(F,10*log10(Pxx)); hold on;
% 
% x_t = smt_8;
% [Pxx,F] = periodogram(x_t,[],length(x_t),f_sample);
% plot(F,10*log10(Pxx));

figure;

x_t = smt_4;
x_t = reshape(x_t, n, []);
P1 = mean((abs(fft(x_t)).^2), 2);

x_t = smt_8;
x_t = reshape(x_t, n, []);
P2 = mean((abs(fft(x_t)).^2), 2);

semilogy(f,P1); hold on;
semilogy(f,P2);

% plot(f,P1); hold on;
% plot(f,P2);

grid on;
title(sprintf('Power Spectrum for M-PSK'));
xlabel('f (Hz)');
ylabel('P (Amplitude^2)');
legend('4-PSK', '8-PSK');