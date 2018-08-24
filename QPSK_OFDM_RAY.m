clear
clc


N = 128000;
signal = randi([0 1],1,N);

for k = 1:N/2
    if signal(2*k-1) == 0 && signal(2*k) ==0
        Data(k) = -1-1*1j;
    elseif signal(2*k-1) == 0 && signal(2*k) ==1
        Data(k) = -1+1*1j;
    elseif signal(2*k-1) == 1 && signal(2*k) ==0
        Data(k) = 1-1*1j;
    else
        Data(k) = 1+1*1j;
    end
end

Bmod = zeros(64,N/128);
for m = 1:N/128
    for n = 1:64
        Bmod(n,m) = Data(64*(m-1)+n);
    end
    Ifft_mod(:,m) = ifft(Bmod(:,m));
end

Tran_sig = zeros(1,N/2);
for x = 1:N/128
    for y = 1:64
       Tran_sig(64*(x-1)+y) = Ifft_mod(y,x);
    end
end

%------Rayleigh channel------
Ts = 50*10^(-9);
Fd = 10;
raychn1 = rayleighchan(Ts,Fd);
raysig = filter(raychn1,Tran_sig);
%------AWGN Channel------
K = 0:0.3:3;
for c = 1:length(K)
    Sigma(c) = sqrt(1/(4*(10^K(c))));
    Noise(c,:) = normrnd(0,Sigma(c),1,N/2);
    Tsig(c,:) = real(raysig)+Noise(c,:)+(imag(raysig)+Noise(c,:))*1j;
%     Tsig(c,:) = Tran_sig;
end


%------Receive signal------
for r = 1:length(K)
    for s = 1:N/128
        for d = 1:64
            Bmod(d,s) = Tsig(r,64*(s-1)+d);
        end
        Fft_dm(:,s,r) = fft(Bmod(:,s));
    end
end

%------detect the bit------
for t = 1:length(K)
    for p = 1:N/128
        for q = 1:64
            if real(Fft_dm(q,p,t))<0
                Resig_I(t,64*(p-1)+q) = 0;
            else
                Resig_I(t,64*(p-1)+q) = 1;
            end
            if imag(Fft_dm(q,p,t))<0
                Resig_Q(t,64*(p-1)+q) = 0;
            else
                Resig_Q(t,64*(p-1)+q) = 1;
            end
        end
    end
    for v = 1:N/2
        R_sig(t,2*v-1) = Resig_I(t,v);
        R_sig(t,2*v) = Resig_Q(t,v);
    end
    result(t,:) = xor(R_sig(t,:),signal);
    BER(t) = sum(result(t,:))/N;
end
figure
semilogy(10*K,BER);grid on;
xlabel('Eb/N0 in dB');
ylabel('BER')



