function [H_iid H_Blk] = UL_IID_MIMO_Channel(Nt,Nr,P,L,Nd,K,Type)
%% Wireless Fading Generation
%% Inputs
% Nt    : Number of BS Antennas
% Nr    : Number of UT Antennas
% L     : Number of Multipath Components due to Delay Spread of the channel 
%        between the MIMO system or  Number of Cells If Type is Multicell scenrio
% Nd    : Number of Data Symbols to Be Transmitted 
% K     : Number of Output Symbols Received
% Type  :  Channel Type
%% Outputs
% H_Blk : iid Channel Matrix in Block Wise
% H_iid : Full iid Channel Matrix 

[H_iid H_Blk] = Downlink_IID_Channel(Nr,Nt,P,L,Nd,K,Type)
end

