L_b = randsrc(24^4, 1, [0 1]);

% BINARY (UNIPOLAR) ENCODING

M = [8];
SNR = 0:2:16;
[X, Y] = meshgrid(M, SNR);
Z = zeros(size(X));

for i=1:size(Z,1)
    for j=1:size(Z,2)
        [~, BER, ~] = M_PSK(L_b, X(i,j), Y(i,j), 'bin', 1);
        Z(i,j)=BER;
    end
end

Z(isinf(Z)) = 10^-3;

figure;
a = axes;
s = mesh(X,Y,Z);
title('BER for M-PAM with binary encoding');
set(a,'ZScale','log')
s.FaceColor = 'flat';
xlabel('M(bits)');
ylabel('SNR(dB)');
zlabel('BER');

% GRAY ENCODING

M = [4 8];
SNR = 0:2:16;

[X, Y] = meshgrid(M, SNR);
Z = zeros(size(X));

for i=1:size(Z,1)
    for j=1:size(Z,2)
        m = X(i,j);
        snr = Y(i,j);
        [~, BER, ~] = M_PSK(L_b, m, snr, 'gray', 0);
        Z(i,j)=BER;

    end
end

Z(isinf(Z)) = 10^-3;

figure;
a = axes;
s = mesh(X,Y,Z);
set(a,'ZScale','log')
title('BER for M-PSK with gray encoding');
s.FaceColor = 'flat';
xlabel('M(bits)');
ylabel('SNR(dB)');
zlabel('BER');

% THEORETICAL BER

M = [4 8];
SNR = 0:2:16;

[X, Y] = meshgrid(M, SNR);
Z = zeros(size(X));

for i=1:size(Z,1)
    for j=1:size(Z,2)
        m = X(i,j);
        snr = Y(i,j);
        if m == 4
            Z(i,j) = 2 * qfunc(sqrt(2 * snr)) * ( 1 - qfunc(sqrt(2 * snr))/2);
        else
            Z(i,j) = 2 * qfunc(sqrt(2 * log2(m) * snr)*sin(pi/m));
        end
    end
end

figure;
a = axes;
s = mesh(X,Y,Z);
title('Theoretical BER for M-PSK');
set(a, 'ZScale', 'log')
s.FaceColor = 'flat';
xlabel('M(bits)');
ylabel('SNR(dB)');
zlabel('BER');