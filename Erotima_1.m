L_b = randsrc(24*4, 1, [0 1]);


% BINARY (UNIPOLAR) ENCODING

M = [4 8];
SNR = 0:2:16;
[X, Y] = meshgrid(M, SNR);
Z = zeros(size(X));

for i=1:size(Z,1)
    for j=1:size(Z,2)
        [~, BER, ~] = M_PSK(L_b, X(i,j), Y(i,j), 'bin', 0);
        %Z(i,j)=log10(BER);
        Z(i,j)=BER;
    end
end

figure;
s = mesh(X,Y,Z);
title('BER for M-PAM with binary encoding');
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
        [~, BER, ~] = M_PSK(L_b, X(i,j), Y(i,j), 'gray', 0);
        %Z(i,j)=log10(BER);
        Z(i,j)=BER;
    end
end

figure;
s = mesh(X,Y,Z);
title('BER for M-PSK with gray encoding');
s.FaceColor = 'flat';
xlabel('M(bits)');
ylabel('SNR(dB)');
zlabel('BER');