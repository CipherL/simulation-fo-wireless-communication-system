clear;

N = 10000;
signal = randi([0 1],1,N);
sig_data = signal;

for i = 1:N
    if signal(i) == 0
        signal(i) = -1;
    else
        signal(i) = 1;
    end
end
%------generate I&Q then resamble------
sig_I = signal(1:2:N);
sig_Q = signal(2:2:N);

for m = 1:N/2
    Dig_I(32*(m-1)+1:32*m) = sig_I(m);
    Dig_Q(32*(m-1)+1:32*m) = sig_Q(m);
end
%-------carray frequence------
f0 =2;
fs = 32;
t = 0:f0/fs:1-f0/fs;
Carr_I = cos(t*2*pi);
Carr_Q = sin(t*2*pi);

%------module------
for k = 1:N
    Data_I(16*(k-1)+1:16*k) = Dig_I(16*(k-1)+1:16*k).*Carr_I;
    Data_Q(16*(k-1)+1:16*k) = Dig_Q(16*(k-1)+1:16*k).*Carr_Q;
end
%------transimit signal+Noise------
K = 0 :0.2:1;
for x = 1:6
    Sigma(x) = sqrt(0.25*32/(10^K(x)));
    Noise(x,:) = normrnd(0,Sigma(x),1,160000);
%     Tsig(x,:) = Data_I-Data_Q*1j;
    Tsig(x,:) = Data_I+Noise(x,:)-(Data_Q+Noise(x,:))*1j;
end

for n = 1:6
    for r = 1:N/2
         Rsig_I(n,32*(r-1)+1:32*r) = fft(real(Tsig(n,32*(r-1)+1:32*r)));
         Rsig_Q(n,32*(r-1)+1:32*r) = fft(imag(Tsig(n,32*(r-1)+1:32*r)));
        if real(Rsig_I(n,32*(r-1)+3)) < 0
            det_I(n,r) = 0;
        else
            det_I(n,r) = 1;
        end
        if imag(Rsig_Q(n,32*(r-1)+3)) >0
            det_Q(n,r) = 1;
        else
            det_Q(n,r) = 0;
        end

    end
    Result(1:2:N) = det_I(n,:);
    Result(2:2:N) = det_Q(n,:);
    result(n,:) = Result;
    Dires(n,:) = xor(sig_data,result(n,:));
    BER(n) = sum(Dires(n,:))/N;    
end
figure
plot(real(Rsig_I(6,:)));
figure
plot(imag(Rsig_Q(6,:)));
figure
plot(10*K,log10(BER));
xlabel('S/N');
ylabel('BER');