clear;

%------Random 0&1 signal------
N= 10000;
Signal = randi([0 1],1,N);
sig_data = Signal;

%------NRZ NUMBER------
for i = 1:N
    if Signal(i) ==0
        Signal(i) = -1;
    else
        Signal(i) = 1;
    end
end

%------Resamble------
for m = 1:N
    for n = 1:32
        signal(32*(m-1)+n) = Signal(m);
    end
end

%------module------
f0 =2;
fs = 32;
t = 0:f0/fs:1-1/fs;
Carrier = cos(t*2*pi);

for k = 1:2*N
    signal(16*(k-1)+1:16*k) = signal(16*(k-1)+1:16*k).*Carrier;
end
figure
plot(signal)
%------Transmit------
K = 0 :0.2:1;
for x = 1:length(K)
    Sigma(x) = sqrt(0.25*32/(10^K(x)));
    Noise(x,:) = normrnd(0,Sigma(x),1,32*N);
    Tsig(x,:) = signal+Noise(x,:);
end

% SNR = 10;
% Tsig = awgn(signal,SNR,'measured','dB');
% Tsig = signal;


% -------FFT------
for j = 1:length(K)
    for r = 1:N
        Rsig(32*(r-1)+1:32*r) = fft(Tsig(j,32*(r-1)+1:32*r));
        if real(Rsig(32*(r-1)+3)) < 0
            det(j,r) = 0;
        else
            det(j,r) = 1;
        end
    end
%     figure
%     plot(real(Rsig));
    Result(j,:) = xor(sig_data,det(j,:));
    BER(j) = sum(Result(j,:))/N;
end
figure
plot(10*K,log10(BER));
xlabel('S/N');
ylabel('BER');

