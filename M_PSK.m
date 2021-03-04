function [out, BER, SER] = M_PSK(inpt, M, SNR, encoding, show_plot)

E_s=1;

R = log2(M);

% We have 4 samples (T_sample) per carrier period (Tc) and
% 2 carrier periods (Tc) per symbol period (T_sym).

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

if show_plot
    tmp = repelem(sm, T_sym);
    figure('Position', [10 10 900 600])
    subplot(4,1,1);
    plot(1:T_sample:T_sample*T_sym*10, tmp(1:T_sample*T_sym*10));
    title("Subplot 1: Original signal")
end

% ----------------------------Modulation-----------------------------------

smt = zeros(length(sm) * T_sym, 1);
t = (0:T_sample:(T_sym - T_sample))';

for i = 1:length(sm)
    for j = 1:length(t)
        smt((i - 1) * T_sym + j) = ...
            g_t * cos(2 * pi * (sm(i)/M)) * cos(2 * pi * f_c * t(j)) - ...
            g_t * sin(2 * pi * (sm(i)/M)) * sin(2 * pi * f_c * t(j));
    end
end

if show_plot
    subplot(4,1,2);
    plot(1:T_sample:T_sample*T_sym*10, smt(1:T_sample*T_sym*10));
    title("Subplot 2: Transmitter's output");
end

% ------------------------------AWGN---------------------------------------

var = (E_s / (2 * log2(M))) * 10 ^ (- SNR / 10);
noise = sqrt(var) * randn(length(smt), 1);

r = smt + noise; % receiving signal

if show_plot
    subplot(4,1,3);
    plot(1:T_sample:T_sample*T_sym*10, r(1:T_sample*T_sym*10));
    title("Subplot 3: Signal with AWGN");
end

%---------------------------Demodulation-----------------------------------    
% Calculates the norm of the projection of vector r to every basis signal,
% which is the inner product <sm(t), φ_j(t)>

r = repmat(r, 1, 2);

% Calculate product sm(t) * φ_j(t)
for i = 1:length(r)
    r(i,1) = r(i,1) * g_t * cos(2*pi*f_c*(i-1));
    r(i,2) = r(i,2) * g_t * -sin(2*pi*f_c*(i-1));
end

% if show_plot
%     subplot(4,1,4);
%     plot(1:T_sample:T_sample*T_sym*10, r(1:T_sample*T_sym*10));
%     title("Subplot 4: Receiver's signal after filtering")
% end

% Calculate inner product <sm(t), φ_j(t)>
% basically it integrates the previously calculated product
demodulated = zeros(size(r,1)/T_sym, size(r,2)); 
for i = 1:size(demodulated,1)
    from = ((i-1)*T_sym)+1;
    to = from + T_sym -1;
    demodulated(i, :) = sum(r(from:to,:));
end

%-----------------------------Decision------------------------------------- 

decision=zeros(size(sm));
for i=1:length(demodulated)
    dist=abs(sum(sqrt((constellation-demodulated(i,:)).^2),2));
    [~,ind]=min(dist);
    decision(i)=ind-1;
end

% Calculate SER
errorsym=0;
totalsymb=0;
for i=1:length(sm)
    if decision(i)~=sm(i)
        errorsym=errorsym+1;
    end
    totalsymb=totalsymb+1;
end


% ---------------------------Demapper--------------------------------------

if strcmp(encoding,'bin')
    blocks_out=de2bi(decision,R,'left-msb');
else
    decision=gray2bin(decision,'psk',M);
    blocks_out=de2bi(decision,R,'left-msb');
end

output_sequence=reshape(blocks_out.',[],1);

k=inpt(inpt~=output_sequence);
BER=length(k)/length(inpt);
out=output_sequence;
SER=errorsym/totalsymb;

end

