function [W W_eq R1 Q] = Precoder_Design_BER(H,SNRdB,Type,nr,nt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[Nr Nt] = size(H);

D = eye(Nt);
K = Nt/nt;
Wm = zeros(Nt,Nr);

switch Type
    case 'MF'
         for ii = 1:Nr
             channelvector = (H(ii,:)*D)'; %Useful channel
             Wm(:,ii) = channelvector;%/norm(channelvector); %Normalization of useful channel
         end
           beta = 1/nt;%sqrt(size(H,2)/trace(Wm*Wm'));
           W = beta*Wm;
    case 'ZF'
        EffChan = (H*D);  % Effective Chnanels
        ChanInv =  EffChan'/(EffChan*EffChan'); % Zero frocing Channel inversion
        for ii = 1:Nr
            W(:,ii) = ChanInv(:,ii);%/norm(ChanInv(:,ii)); % Normalization of Zero Forcing
        end
    
    case 'MFOFDM'
         F = dftmtx(Nt/nt)/sqrt(Nt/nt);
         Et = kron(F,eye(nt));
         Er = kron(F,eye(nr));
         Hbd = Er*H*Et';
         W_eq = Precoder_Design_BER(Hbd,SNRdB,'MF',nr,nt);
         Wm = Et';
         beta = 1;%sqrt((Nt)/trace(Wm*Wm'));
         W = beta*Wm;
    case 'ZFOFDM'
         F = dftmtx(Nt/nt)/sqrt(Nt/nt);
         Et = kron(F,eye(nt));
         Er = kron(F,eye(nr));
         Hbd = Er*H*Et';
         W_eq = Precoder_Design_BER(Hbd,SNRdB,'ZF',nr,nt);
         Wm = Et';
         beta = 1;%sqrt((Nt)/trace(Wm*Wm'));
         W = beta*Wm;
     case 'MF_RFDM_QR'
         [~,~,~,R] = nonOrtho_ramanujanBaseMtrx(Nt/nt);
         Et = kron(R,eye(nt));
         Er = kron(R,eye(nr));
         
         Hbd = Er'*H*Et;
         %Wmf = Precoder_Design(Hbd,SNRdB,'MF');
         bb=sqrt((Nt)/trace(Hbd*Hbd'));
         [Q R1] = qr(Hbd);
         Wqr  = bb*Et;
         Wm = Wqr;%*Wzf;%*Wmf;
         beta = 1;%sqrt((Nt)/trace(Wm*Wm'));
         W = beta*Wm;         
      
     case 'MFRFDM'
         [~,~,~,R] = nonOrtho_ramanujanBaseMtrx(Nt/nt);
         Et = kron(R,eye(nt));
         Er = kron(R,eye(nr));
         Hbd = Er'*H*Et;
         W_eq = Precoder_Design_BER(Hbd,SNRdB,'MF',nr,nt);
         Wm =  Et;
         beta = 1;%sqrt((Nt)/trace(Wm*Wm'));
         W = beta*Wm;
    case 'ZFRFDM'
         [~,~,~,R] = nonOrtho_ramanujanBaseMtrx(Nt/nt);
         Et = kron(R,eye(nt));
         Er = kron(R,eye(nr));
         Hbd = Er'*H*Et;
         W_eq = Precoder_Design_BER(Hbd,SNRdB,'ZF',nr,nt);
         Wm  = Et;%*Wzf;
         beta = 1;%sqrt((Nt)/trace(Wm*Wm'));
         W = beta*Wm;   

    otherwise
        disp('No Match')
end

