clear;

%------Random 0&1 signal------
N= 64000;
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

%------S/P------
Bmod = zeros(64,N/64);
for m = 1:N/64
    for n = 1:64
        Bmod(n,m) = Signal(64*(m-1)+n);
    end
    Ifft_mod(:,m) = ifft(Bmod(:,m));
end

Tran_sig = zeros(1,N);
for x = 1:N/64
    for y = 1:64
        Tran_sig(64*(x-1)+y) = Ifft_mod(y,x);
    end
end


%------AWGN Channel------
K = 0:0.3:3;
for c = 1:length(K)
    Sigma(c) = sqrt(1/(2*(10^K(c))));
    Noise(c,:) = normrnd(0,Sigma(c),1,N);
    Tsig(c,:) = real(Tran_sig)+Noise(c,:)+(imag(Tran_sig)+Noise(c,:))*1j;
%     Tsig(c,:) = Tran_sig;
end


%------Receive signal------
for r = 1:length(K)
    for s = 1:N/64
        for d = 1:64
            Bmod(d,s) = Tsig(r,64*(s-1)+d);
        end
        Fft_dm(:,s,r) = fft(Bmod(:,s));
    end
end

%------detect the bit------
for t = 1:length(K)
    for p = 1:N/64
        for q = 1:64
            if real(Fft_dm(q,p,t))<0
                Resig(t,64*(p-1)+q) = 0;
            else
                Resig(t,64*(p-1)+q) = 1;
            end
        end
    end
    result(t,:) = xor(Resig(t,:),sig_data);
    BER(t) = sum(result(t,:))/N;
end
figure
semilogy(10*K,BER);grid on;
xlabel('Eb/N0 in dB');
ylabel('BER')