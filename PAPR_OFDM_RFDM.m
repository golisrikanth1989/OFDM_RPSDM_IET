% plot_CCDF.m: % Plot the CCDF curves of Fig. 7.3.
clear all; clc; clf
Ns=2.^[3 6:9]; b=4; M=2^b; Nblk=1e4; zdBs=[0:0.1:15];
N_zdBs=length(zdBs);
CCDF_formula=inline('1-((1-exp(-z.^2/(2*s2))).^(0.001*N))','N','s2','z'); %(7.9)


%R=riemann(N+1);
%%
for n = 1:length(Ns)
N=Ns(n); x = zeros(Nblk,N); sqN=sqrt(N);
%R=riemann(N+1); 
NoD = divisors(N);
Si = round(NoD/2);
Nr = length(Si);
Nod = floor(NoD/2)+1;

N_Mbits = Nr*log2(16)*N;
N_Ibits = sum(log2(Si))*N;


if N == 256
   load Ramanujan_256;
   R_Orth_N  =R;
elseif N==128
    load Ramanujan_128;
   R_Orth_N  =R;
elseif N==512
    load Ramanujan_512;
    R_Orth_N  =R;
else
[nonOrthoBasesMtrx,orthoBasesMtrx,ramSqnce,R_Orth_N] = nonOrtho_ramanujanBaseMtrx(N);
%Scal = diag(ramSqnce(length(NoD)).enerGy).*sqN;
end


%%
x = [];
x1 = [];
CFx=[];
CFx1=[];
X=0;
X1=0;
for k=1:Nblk
% Conventional OFDM
X=mapper(b,N);
%X =Gcd_Mod_r_Sequence(N,X1).';
x(k,:)=ifft(X,N)*sqN; CFx(k)=PAPR(x(k,:));
x1(k,:)=(R_Orth_N*X.').'; CFx1(k)=PAPR(x1(k,:));


% [Data0 Data1] = Shuffle_Index_Ramanujan(N);
% x2(k,:)=ifft((Data0'.*X),N)*sqN; CFx2(k)=PAPR(x2(k,:));
% x3(k,:)=(R_Orth_N*(Data0'.*X).'); CFx3(k)=PAPR(x3(k,:));
% 
% x4(k,:)=ifft((Data1'.*X),N)*sqN; CFx4(k)=PAPR(x4(k,:));
% x5(k,:)=(R_Orth_N*(Data1'.*X).'); CFx5(k)=PAPR(x5(k,:));



%end
fprintf('%d \n',k);
end


s2 = mean(mean(abs(x1)))^2/(pi/2);
CCDF_theoretical=CCDF_formula(N,s2,10.^(zdBs/20));
for i=1:N_zdBs, CCDF_simulated(i)=sum(CFx>zdBs(i))/Nblk; end
for i=1:N_zdBs, CCDF_simulated1(i)=sum(CFx1>zdBs(i))/Nblk; end
% for i=1:N_zdBs, CCDF_simulated2(i)=sum(CFx2>zdBs(i))/Nblk; end
% for i=1:N_zdBs, CCDF_simulated3(i)=sum(CFx3>zdBs(i))/Nblk; end
% for i=1:N_zdBs, CCDF_simulated4(i)=sum(CFx4>zdBs(i))/Nblk; end
% for i=1:N_zdBs, CCDF_simulated5(i)=sum(CFx5>zdBs(i))/Nblk; end
%for i=1:N_zdBs, CCDF_simulated2(jj,i)=sum(CFx2(jj,:)>zdBs(i))/Nblk; end
%semilogy(zdBs,CCDF_simulated2(jj,:),'r:*'); hold on

% [minpapr index] = min(CFx1);
% [minpp jjindex] = min(minpapr);
% jj = index(jjindex);

% [minpapr1 index1] = min(CFx2);
% [minpp1 jjindex1] = min(minpapr1);
% jj1 = index1(jjindex1);


%semilogy(zdBs,CCDF_theoretical,'k-'); hold on; grid on;
semilogy(zdBs(3:3:end),CCDF_simulated(3:3:end),'rO-','LineWidth',1.4);hold on;grid on;
semilogy(zdBs(3:3:end),CCDF_simulated1(3:3:end),'bO-','LineWidth',1.4);
% semilogy(zdBs(1:3:end),CCDF_simulated2(1:3:end),'b+-','LineWidth',1.2);
% semilogy(zdBs(1:3:end),CCDF_simulated3(1:3:end),'r+-','LineWidth',1.2);
% semilogy(zdBs(1:3:end),CCDF_simulated4(1:3:end),'rO--','LineWidth',1.2);
% semilogy(zdBs(1:3:end),CCDF_simulated5(1:3:end),'bO--','LineWidth',1.2);
%semilogy(zdBs,CCDF_simulated2(jj1,:),'r:*');

end
%%
axis([zdBs([1 end-30]) 1e-3 1]); %title('OFDM system with N-point FFT');
xlabel('z[dB]'); ylabel('Pr(\chi_{dB} > z[dB])'); legend('OFDM','RPSDM');
%xlabel('z[dB]'); ylabel('Pr(PAPR > z[dB])'); legend('OFDM','RPSDM','SIM-OFDM-LD','SIM-RPSDM-LD','SIM-OFDM-HD','SIM-RPSDM-HD');
%% Proper and Improper
% close all;
% xx=0;
% xx1=0;
% cxx=0;
% cxx1=0;
% for jj =1:Nblk
% xx = xx+conv(x(jj,:)',fliplr(x(jj,:)));
% xx1 = xx1+ conv(x1(jj,:)',fliplr(x1(jj,:)));
% 
% cxx = cxx+conv(x(jj,:).',fliplr(x(jj,:)));
% cxx1 = cxx1+ conv(x1(jj,:).',fliplr(x1(jj,:)));
% end
% figure;
% %bar3(abs(xx)/Nblk);
% plot(abs(xx)/Nblk);hold on;
% plot(abs(cxx)/Nblk);
% %figure;
% plot(abs(xx1)/Nblk);
% plot(abs(cxx1)/Nblk);
% %bar3(abs(xx1)/Nblk);