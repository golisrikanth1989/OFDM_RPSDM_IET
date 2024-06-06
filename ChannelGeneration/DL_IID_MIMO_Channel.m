function [H_iid H_Blk] = DL_IID_MIMO_Channel(Nt,Nr,P,L,Nd,K,Type)
%% Wireless Fading Generation
%% Inputs
% Nt    : Number of Transmitting Antennas
% Nr    : Number of Receving Antennas
% L     : Number of Multipath Components due to Delay Spread of the channel 
%        between the MIMO system or  Number of Cells If Type is Multicell scenrio
% Nd    : Number of Data Symbols to Be Transmitted 
% K     : Number of Output Symbols Received
% Type  :  Channel Type
%% Outputs
% H_Blk : iid Channel Matrix in Block Wise
% H_iid : Full iid Channel Matrix 
switch(Type)
case {'MIMO Linear Time Varying'}
%% MIMO Linear Time Varying
    H_Sub = [];
       H_Blk = cell(K,0);
        
       for jj =0:Nd-1
           H = cell(K,1);
       for kk = 0:K-1
           if kk<=L-1
           H{kk+1,1} =  P(kk+1)^(1/2).*sqrt(1/2).*(randn(Nr,Nt)+j*randn(Nr,Nt));
           else
           H{kk+1,1} =  zeros(Nr,Nt);
           end
       end
       %H_Blk=[];
           
               H1 = circshift(H,jj);
               H_Blk = [H_Blk H1];
     
       end
       H_iid = cell2mat(H_Blk);
case {'Block Toeplitz'}
%% MIMO Selective Block Toeplitz
    H_Sub = [];
       H_Blk = cell(K,0);
       H = cell(K,1); 
       for kk = 0:K-1
           if kk<=L-1
           H{kk+1,1} =  P(kk+1)^(1/2).*sqrt(1/2).*(randn(Nr,Nt)+j*randn(Nr,Nt));
           else
           H{kk+1,1} =  zeros(Nr,Nt);
           end
       end
       %H_Blk=[];
       for kk =0:Nd-1
           H1 = circshift(H,kk);
           H_Blk = [H_Blk H1];
       end
       H_iid = cell2mat(H_Blk);
case {'Block Circualnt'}
%% MIMO Selective  Block Circulant
    H_Sub = [];
       H_Blk = cell(K,0);
       H = cell(K,1); 
       for kk = 0:K-1
           if kk<=L-1
           H{kk+1,1} =  P(kk+1)^(1/2).*sqrt(1/2).*(randn(Nr,Nt)+j*randn(Nr,Nt));
           else
           H{kk+1,1} =  zeros(Nr,Nt);
           end
       end
       %H_Blk=[];
       for kk =0:K-1 
           H1 = circshift(H,kk);
           H_Blk = [H_Blk H1];
       end
       H_iid = cell2mat(H_Blk);
case {'Toeplitz Blocks'}
%%  MIMO Selective Toeplitz Blocks
   H_Sub = [];
   H_Blk = cell(Nr,Nt);
   for rr =0:Nr-1
       for tt = 0:Nt-1
           hsub = [P.^(1/2).*sqrt(1/2).*(randn(1,L)+j*randn(1,L)) zeros(1,K-L)].';
           %H = [hsub zeros(1,tones-length(hsub))].';
           H = hsub;
           Hbb=[];
           for kk =0:Nd-1
               H1 = circshift(H,kk);
               Hbb = [Hbb H1];
           end
               H_Blk{rr+1,tt+1} =  Hbb;
       end
   end
   H_iid = cell2mat(H_Blk);
case {'Circualnt Blocks'}
%%  MIMO Selective Circulant Blocks
   H_Sub = [];
   H_Blk = cell(Nr,Nt);
   for rr =0:Nr-1
       for tt = 0:Nt-1
           hsub = [P.^(1/2).*sqrt(1/2).*(randn(1,L)+j*randn(1,L)) zeros(1,K-L)].';
           %H = [hsub zeros(1,tones-length(hsub))].';
           H = hsub;
           Hbb=[];
           for kk =0:K-1
               H1 = circshift(H,kk);
               Hbb = [Hbb H1];
           end
               H_Blk{rr+1,tt+1} =  Hbb;
       end
   end
   H_iid = cell2mat(H_Blk);
   
case {'MIMOMultiCell'}
%% MIMO Spatial Multiplexing Channel
   H_Sub = [];
   H_Blk = cell(Nr,Nt);
     for rr =0:L-1
       for tt = 0:L-1
           H = sqrt(1/2).*(randn(Nr,Nt)+1j*randn(Nr,Nt));
           H_Blk{rr+1,tt+1} =  H;
       end
   end
H_iid = cell2mat(H_Blk);

    otherwise
        disp('No type Selected');
end
end

