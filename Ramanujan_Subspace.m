%% Author: Goli Srikanth
%% Ph.D Scholar
%% Shiv Nadar University
clc;
clear all;
close all;
%%
N = 8;
D  = dftmtx(N)';
rr=32;
aa=0;
x = exp(1j*2*pi/10);
k1 = linspace(-pi,pi,rr*N);
[A,B,~,R] = nonOrtho_ramanujanBaseMtrx(N);
ii=1;
for jj= floor(divisors(N)/2)+1
signal_R(jj,:) = abs(fft(repmat(B(:,jj)*x,[1,1]),rr*N));
 aa = aa+signal_R(jj,:);
plot(k1,fftshift((signal_R(jj,:)))/N,'LineWidth',1); hold on; %grid on;
if jj == 1
ii = ii;
else
ii=ii+1;    
end
end
set(gca,'XTick',-pi:pi/4:pi);
set(gca,'XTickLabel',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'})
set(gca,'YTick',0:0.5:1);
set(gca,'FontSize',12);
xlabel('Frequency');
ylabel('Magnitude');
xlim([-pi,pi]);
legend('\bf{S}_1','\bf{S}_2','\bf{S}_4','\bf{S}_8');
% plot(10*log10(aa),'LineWidth',1.5);