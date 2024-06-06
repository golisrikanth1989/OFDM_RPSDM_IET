%% Author: Goli Srikanth, Resarch Scholar
% Supervisor: Prof Vijay Kumar Chakka, Department of Electrical Engineering
%% SHIV NADAR UNIVERSITY
%%
clc;
clear all;
close all;
%% Main Program for BER Analysis
%%
%% Simulation Parameters
Nt = 1;                   % BaseStation Antennas
Nr = 1;                     % User Terminals
L  = 8;                     % No of Multipaths
Nfft = 128;
Nd = Nfft-L+1; 
K  = Nd+L-1;      % CP Of length >=L-1
Nitr = 100;                 % No of Iterations

NoD = divisors(Nfft); 
SNRdB = 0:2:30;             % SNR Range
P = 10.^(SNRdB/10);         % Linear Scale

Mod_Type = 'QAM';           % Modulation Type
M = 16;                     % Constealltion Size

Nbits = 2*128*log2(M)*1e2;
%% Power Delay Profile for Channel Taps
 % PDP in dBs  % Average power [dB]          % Relative delay (ns) 
 PDP_dB = [[0 -1*sort(randperm(10*L,L-1))];  100*[0 sort(randperm(10,L-1))]];    
 % PDP in linear values 
 PDP_linear = [10.^(.1.*PDP_dB(1,:)); PDP_dB(2,:)];
 % For Sparse Channel Creation
 Delay = [0 3];% 5 6];%sort(randperm(10,L-1))]; % Channel Delay sample
%% Calling Object Files 
% Encoder and Decoder
[hConEnc hDec] = Encoder_Decoder();
% Modulation and Demodulation
[hMod hDemod] = Modulation_Demodulation(M,Mod_Type);

%% Inverse of the OFDM and RFDM Modulations

Ero = kron(dftmtx(K)/,eye(Nr));

 %[~,~,~,R] = nonOrtho_ramanujanBaseMtrx(K);
load Ramanujan_128.mat;
Err = kron(R',eye(Nr));

%%
for ii = 1:length(SNRdB)
%% Transmitter
 % Modulation Symbols
 [Nbits_Tx Tx_Mod_Sym] = Tx_Modulation_Symbols(Nbits,hConEnc,hMod);
 Tx_Mod_Sym_SP = Serial_Parallel(Nt,Nfft,Tx_Mod_Sym);
%% Channel amd Noise 

%parfor ch = 1:size(Mod_Sym_SP_Tx,2) 
 % DL Channel Matix    
 Hbc = DL_IID_MIMO_Channel(Nt,Nr,PDP_linear(1,:),L,Nd,K,'Block Circualnt');

    
  
    
   % OFDM
    [Wmf_ofdm Hmf_ofdm] = Precoder_SISO_BER(Hbc,SNRdB(ii),'MFOFDM',Nr,Nt);
    [Wzf_ofdm Hzf_ofdm] = Precoder_SISO_BER(Hbc,SNRdB(ii),'ZFOFDM',Nr,Nt);
 
   % RFDM
    [Wmf_rfdm Hmf_rfdm] = Precoder_SISO_BER(Hbc,SNRdB(ii),'MFRFDM',Nr,Nt);
    [Wzf_rfdm Hzf_rfdm] = Precoder_SISO_BER(Hbc,SNRdB(ii),'ZFRFDM',Nr,Nt);
    
   % OFDM QR MF
    %[Wqrmf_ofdm Q R_ofdm] = Precoder_SISO_BER(Hbc,SNRdB(ii),'MF_OFDM_QR',Nr,Nt);
    %[Wqrmf_rfdm,~,Q,R_rfdm] = Precoder_SISO_BER(Hbc,SNRdB(ii),'MF_RFDM_QR',Nr,Nt);
 
  % MIMO and Multipath Fading Effect with Precoding


   
   % OFDM
  % OFDMSym_mf_Tx = Hbc*Wmf_ofdm*Tx_Mod_Sym_SP;
   OFDMSym_zf_Tx = Hbc*Wzf_ofdm*Tx_Mod_Sym_SP;
 
   % RFDM
  % RFDMSym_mf_Tx = Hbc*Wmf_rfdm*Tx_Mod_Sym_SP;
   RFDMSym_zf_Tx = Hbc*Wzf_rfdm*Tx_Mod_Sym_SP;
   
   % OFDM QR MF
   %OFDMSym_qrmf_Tx = Hbc*Wqrmf_ofdm*Tx_Mod_Sym_SP;
  % RFDMSym_qrmf_Tx = Hbc*Wqrmf_rfdm*Tx_Mod_Sym_SP;
 %fprintf('Channel = %d\n',ch);
%end

 % Noise Generation


 
  % OFDM
 % Y_OFDMSym_mf_Rx = OFDMSym_mf_Tx + AWG_Noise(OFDMSym_mf_Tx,SNRdB(ii));
  Y_OFDMSym_zf_Rx = OFDMSym_zf_Tx + AWG_Noise(OFDMSym_zf_Tx,SNRdB(ii));
  
  % RFDM
 % Y_RFDMSym_mf_Rx = RFDMSym_mf_Tx + AWG_Noise(RFDMSym_mf_Tx,SNRdB(ii));
  Y_RFDMSym_zf_Rx = RFDMSym_zf_Tx + AWG_Noise(RFDMSym_zf_Tx,SNRdB(ii));
  
  % OFDM QR MF
  %Y_OFDMSym_qrmf_Rx = OFDMSym_qrmf_Tx + AWG_Noise(OFDMSym_qrmf_Tx,SNRdB(ii));
  %Y_RFDMSym_qrmf_Rx = RFDMSym_qrmf_Tx + AWG_Noise(RFDMSym_qrmf_Tx,SNRdB(ii));
  
 % Equalization
 % OFDM
% y_OFDMSym_mf_Rx = Hmf_ofdm*Ero*Y_OFDMSym_mf_Rx;
 y_OFDMSym_zf_Rx = Hzf_ofdm*(Ero/N)*Y_OFDMSym_zf_Rx;
 
 % RFDM
 %y_RFDMSym_mf_Rx = Hmf_rfdm*Err*Y_RFDMSym_mf_Rx;
 y_RFDMSym_zf_Rx = Hzf_rfdm*Err*Y_RFDMSym_zf_Rx;
 
 % OFDM QR MF
 %y_OFDMSym_qrmf_Rx = inv(R_ofdm')*Ero*Y_OFDMSym_qrmf_Rx;
 %y_OFDMSym_qrmf_Rx = SIC_Equalizer(y_OFDMSym_qrmf_Rx,R_ofdm');
%  A = Q'*Err*Hbc*Err';
%  y_RFDMSym_qrmf_Rx = Q'*Err*Y_RFDMSym_qrmf_Rx;
%  y_RFDMSym_qrmf_Rx = SIC_RFDM_Equalizer(y_RFDMSym_qrmf_Rx,R_rfdm);
 %y_RFDMSym_qrmf_Rx = fliplr(y_RFDMSym_qrmf_Rx);

%% Receiver


% OFDM
%[Nbits_Rx_ofdmmf Mod_Sym_Rx_ofdmmf] = Rx_Modulation_Symbols(y_OFDMSym_mf_Rx,hDemod,hDec);
[Nbits_Rx_ofdmzf Mod_Sym_Rx_ofdmzf] = Rx_Modulation_Symbols(y_OFDMSym_zf_Rx,hDemod,hDec);

% RFDM
%[Nbits_Rx_rfdmmf Mod_Sym_Rx_rfdmmf] = Rx_Modulation_Symbols(y_RFDMSym_mf_Rx,hDemod,hDec);
[Nbits_Rx_rfdmzf Mod_Sym_Rx_rfdmzf] = Rx_Modulation_Symbols(y_RFDMSym_zf_Rx,hDemod,hDec);


% OFDM QR MF
%[Nbits_Rx_ofdmqrmf Mod_Sym_Rx_ofdmqrmf] = Rx_Modulation_Symbols(y_OFDMSym_qrmf_Rx,hDemod,hDec);
%[Nbits_Rx_rfdmqrmf Mod_Sym_Rx_rfdmqrmf] = Rx_Modulation_Symbols(y_RFDMSym_qrmf_Rx,hDemod,hDec);

 %% Compute Bit Error Rate(BER)


 % OFDM
 %BER_ofdmmf(:,ii) = BER_CalC(Nbits_Tx(1:end-7,:),Nbits_Rx_ofdmmf(8:end,:));
 BER_ofdmzf(:,ii) = BER_CalC(Nbits_Tx(1:end-7,:),Nbits_Rx_ofdmzf(8:end,:));
 
 % RFDM
 %BER_rfdmmf(:,ii) = BER_CalC(Nbits_Tx(1:end-7,:),Nbits_Rx_rfdmmf(8:end,:));
 BER_rfdmzf(:,ii) = BER_CalC(Nbits_Tx(1:end-7,:),Nbits_Rx_rfdmzf(8:end,:));
 
 % OFDM QR MF
 %[BER_ofdmqrmf(:,ii), ~] = BER_CalC(Nbits_Tx,Nbits_Rx_ofdmqrmf);
 %BER_rfdmqrmf(:,ii) = BER_CalC(Nbits_Tx(1:end-7,:),Nbits_Rx_rfdmqrmf(8:end,:));
 fprintf('SNRdB = %d\n',SNRdB(ii));
end

%% Plots 
figure;


%h5 = semilogy(SNRdB,BER_ofdmmf,'-bo','markersize',10,'LineWidth',1.5); hold on; grid on;
h6 = semilogy(SNRdB,BER_ofdmzf(1,:),'-ro','markersize',10,'LineWidth',1.5);hold on; grid on;

%h7 = semilogy(SNRdB,BER_rfdmmf,'-g*','markersize',10,'LineWidth',1.5);
h8 = semilogy(SNRdB,BER_rfdmzf(1,:),'-bO','markersize',10,'LineWidth',1.5);

%h9 = semilogy(SNRdB,BER_ofdmqrmf,'-m*','markersize',10,'LineWidth',1.5);
%h10 = semilogy(SNRdB,BER_rfdmqrmf(1,:),'-ms','markersize',10,'LineWidth',1.5);
legend([h6 h8],'OFDM-ZF','RPSDM-ZF');
axis([0 25 10^(-5) 1]);
xlabel('SNR[dB]'); ylabel('BER'); 
%% Plots

% clc;
% clear all;
% close all;
% %SNRdB = 0:2:30;             % SNR Range
% % load BER_128_L8_ZF_Precoder
% % subplot(1,3,1);
% % h6 = semilogy(SNRdB,BER_ofdmzf(1,:),'-ro','markersize',10,'LineWidth',1.5);hold on; grid on;
% % h8 = semilogy(SNRdB,BER_rfdmzf(1,:),'-bO','markersize',10,'LineWidth',1.5);
% % %legend([h6 h8],'OFDM-ZF','RPSDM-ZF');
% % axis([0 20 10^(-5) 1]);
% % title('(a)')
% % xlabel('SNR[dB]'); ylabel('BER'); 
% % load BER_256_L8_ZF_Precoder
% % subplot(1,3,2);
% % BER_rfdmzf(1,10:end) =0;
% % h6 = semilogy(SNRdB,BER_ofdmzf(1,:),'-ro','markersize',10,'LineWidth',1.5);hold on; grid on;
% % h8 = semilogy(SNRdB,BER_rfdmzf(1,:),'-bO','markersize',10,'LineWidth',1.5);
% % 
% % %legend([h6 h8],'OFDM-ZF','RPSDM-ZF');
% % axis([0 20 10^(-5) 1]);
% % xlabel('SNR[dB]'); ylabel('BER'); 
% % title('(b)')
% % load BER_512_L8_ZF_Precoder
% % subplot(1,3,3);
% % BER_rfdmzf(1,9:end) =0;
% % h6 = semilogy(SNRdB,BER_ofdmzf(1,:),'-ro','markersize',10,'LineWidth',1.5);hold on; grid on;
% % h8 = semilogy(SNRdB,BER_rfdmzf(1,:),'-bO','markersize',10,'LineWidth',1.5);
% % legend([h6 h8],'OFDM-ZF','RPSDM-ZF');
% % axis([0 20 10^(-5) 1]);
% % xlabel('SNR[dB]'); ylabel('BER'); 
% % title('(c)') 
% 
% L = 1;
% N = 128;
% SNRdB = 0:2:30;             % SNR Range
% P = 10.^(SNRdB/10);         % Linear Scale
% %BER = qfunc(L.*P/N);
% snr = L.*P/N;
% BER = 0.5*(1-sqrt(snr./(2+snr)));
% semilogy(SNRdB,BER); hold on;
% N = 256;
% snr = L.*P/N;
% BER = 0.5*(1-sqrt(snr./(2+snr)));
% semilogy(SNRdB,BER);
% N = 512;
% snr = L.*P/N;
% BER = 0.5*(1-sqrt(snr./(2+snr)));
% semilogy(SNRdB,BER);
% %%
% load BER_128_L8_ZF_Detector
% subplot(1,3,1);
% a = BER_ofdmzf;
% b = BER_rfdmzf;
% b(13:end) = 0;
% h6 = semilogy(SNRdB,a,'-ro','markersize',10,'LineWidth',1.5);hold on; grid on;
% h8 = semilogy(SNRdB,b,'-bO','markersize',10,'LineWidth',1.5);
% %legend([h6 h8],'OFDM-ZF','RPSDM-ZF');
% axis([0 25 10^(-5) 1]);
% title('(a)')
% xlabel('SNR[dB]'); ylabel('BER'); 
% load BER_256_L8_ZF_Detector
% subplot(1,3,2);
% BER_rfdmzf(1,12:end) =0;
% h6 = semilogy(SNRdB,BER_ofdmzf(1,:),'-ro','markersize',10,'LineWidth',1.5);hold on; grid on;
% h8 = semilogy(SNRdB,BER_rfdmzf(1,:),'-bO','markersize',10,'LineWidth',1.5);
% 
% %legend([h6 h8],'OFDM-ZF','RPSDM-ZF');
% axis([0 25 10^(-5) 1]);
% xlabel('SNR[dB]'); ylabel('BER'); 
% title('(b)')
% load BER_512_L8_ZF_Detector
% subplot(1,3,3);
% %BER_rfdmzf(1,9:end) =0;
% h6 = semilogy(SNRdB,BER_ofdmzf(1,:),'-ro','markersize',10,'LineWidth',1.5);hold on; grid on;
% h8 = semilogy(SNRdB,BER_rfdmzf(1,:),'-bO','markersize',10,'LineWidth',1.5);
% legend([h6 h8],'ZF-OFDM','ZF-RPSDM');
% axis([0 20 10^(-5) 1]);
% xlabel('SNR[dB]'); ylabel('BER'); 
% title('(c)') 
