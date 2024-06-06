function [W R1] = Precoder_Design(H,SNRdB,Type,nr,nt,L)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[Nr Nt] = size(H);

D = eye(Nt);

Wm = zeros(Nt,Nr);
K = Nt/nt;
switch Type
    case 'MF'
         for ii = 1:Nr
             channelvector = (H(ii,:))'; %Useful channel
             Wm(:,ii) = channelvector/norm(channelvector); %Normalization of useful channel
         end
           beta = 1;%sqrt((nt*L)/trace(Wm*Wm'));
           W = beta*Wm;
    case 'ZF'
        EffChan = (H);  % Effective Chnanels
        ChanInv =  EffChan'/(EffChan*EffChan'); % Zero frocing Channel inversion
        for ii = 1:Nr
           Wm(:,ii) = ChanInv(:,ii)/norm(ChanInv(:,ii)); % Normalization of Zero Forcing
        end
           beta = 1;%sqrt((nt*L)/trace(Wm*Wm'));
           W = beta*Wm;
    case 'MMSE'
         % Signal-to-noise ratio (SNR) manipulation
          sigma_s = 1; % signal square root power
          snr = exp(SNRdB*log(10)/10); % conversion from dB  
          sigma_n = sqrt( (Nr)/snr)*sigma_s; % noise standard deviation
          Wm = H'*inv(H*H' + (sigma_n^2/sigma_s^2)*eye(Nr)); 
          beta = sqrt(size(H,1)/trace(Wm*Wm'));
          W = beta*Wm;
    case 'RCI'
%         % Signal-to-noise ratio (SNR) manipulation
%           sigma_s = 1; % signal square root power
%           snr = exp(SNRdB*log(10)/10); % conversion from dB  
%           sigma_n = sqrt( (Nr)/snr)*sigma_s+rand(1); % noise standard deviation
%           Wm = H'*inv(H*H' + (sigma_n^2/sigma_s^2)*eye(Nr)); 
%           %Wm = H'*inv(H*H' + 0.01*eye(Nr)); 
%           beta = sqrt(size(H,1)/trace(Wm*Wm'));
%           W = beta*Wm;
         Sigma2 = Nt*0.5*SNRdB; Sigma = sqrt(Sigma2);
         alpha = csi;%Nr/SNRdB+sqrt(1/2);%+rand(1); %
         Wm = H'*inv(H*H' + Nt*alpha*eye(Nr));
         beta = sqrt((Nt)/trace(Wm*Wm'));
         c = sqrt((Nr)/trace(H'*(H*H' + Nt*alpha*eye(Nr))^(-2)*H));
         W = c*Wm;
    case 'MFOFDM'
         F = dftmtx(K)/sqrt(K);
         Et = kron(F,eye(nt));
         Er = kron(F,eye(nr));
         Hbd = Er*H*Et';
         Wmf = Precoder_Design(Hbd,SNRdB,'MF',nr,nr,K);
         Wm = Et'*Wmf;
         beta = 1;%sqrt((Nt)/trace(Wm*Wm'));
         W = beta*Wm;
    case 'ZFOFDM'
         F = dftmtx(K)/sqrt(K);
         Et = kron(F,eye(nt));
         Er = kron(F,eye(nr));
         Hbd = Er*H*Et';
         Wzf = Precoder_Design(Hbd,SNRdB,'ZF',nr,nt,K);
         Wm = Et'*Wzf;
         beta = 1;%sqrt((Nt)/trace(Wm*Wm'));
         W = beta*Wm;
     case 'MF_RFDM_QR'
         %[~,~,~,R] = nonOrtho_ramanujanBaseMtrx(K);
         load Ramanujan_32
         Et = kron(R,eye(nt));
         Er = kron(R,eye(nr));
         Hbd = Er'*H*Et;
         Wmf = Precoder_Design(Hbd,SNRdB,'MF',nr,nt,K);
         bb=1;%sqrt((nt*L)/trace(Hbd*Hbd'));
         [Q R1] = qr(bb*Hbd*Wmf);
         Wqr  = bb*Et*Wmf*Q;
         Wzf = Precoder_Design(R1',SNRdB,'ZF',nr,nt,K);
         Wm = Wqr;%*Wzf;
         beta = 1;%sqrt((Nt)/trace(Wm*Wm'));
         W = beta*Wm;    
     case 'MF_OFDM_QR'
         F = dftmtx(K)/sqrt(K);
         Et = kron(F,eye(nt));
         Er = kron(F,eye(nr));
         Hbd = Er*H*Et';
         Wmf = Precoder_Design(Hbd,SNRdB,'MF',nr,nt,K);
         bb=1;%sqrt((nr)/trace(Hbd*Hbd'));
         [Q R1] = qr(bb*Hbd*Wmf);
         Wqr  = bb*Et'*Wmf*Q;
         Wzf = Precoder_Design(R1',SNRdB,'ZF',nr,nt,K);
         Wm = Wqr;%*Wzf;
         beta = 1;%sqrt((Nt)/trace(Wm*Wm'));
         W = beta*Wm; 
         
     case 'MFRFDM'
         %[~,~,~,R] = nonOrtho_ramanujanBaseMtrx(K);
         load Ramanujan_32
         Et = kron(R,eye(nt));
         Er = kron(R,eye(nr));
         Hbd = Er'*H*Et;
         Wmf = Precoder_Design(Hbd,SNRdB,'MF',nr,nt,K);
         Wm =  Et*Wmf;
         beta = 1;%sqrt((Nt)/trace(Wm*Wm'));
         W = beta*Wm;
    case 'ZFRFDM'
         %[~,~,~,R] = nonOrtho_ramanujanBaseMtrx(K);
         load Ramanujan_32
         Et = kron(R,eye(nt));
         Er = kron(R,eye(nr));
         Hbd = Er'*H*Et;
         Wzf = Precoder_Design(Hbd,SNRdB,'ZF',nr,nt,K);
         Wm  = Et*Wzf;
         beta = 1;%sqrt((Nt)/trace(Wm*Wm'));
         W = beta*Wm;
         
         
         

    otherwise
        disp('No Match')
end

