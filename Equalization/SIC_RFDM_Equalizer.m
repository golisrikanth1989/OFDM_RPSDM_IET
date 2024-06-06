function [X_est] = SIC_RFDM_Equalizer(y,R)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[Nr Nt]=size(R);
[xx yy]=size(y);
X = zeros(size(y));
X_est = zeros(size(y));

for ii = 1:yy
for jj = Nr:-1:1
X(jj,ii)= y(jj,ii)/R(jj,jj);

X_est(jj,ii) = QAM16_slicer(X(jj,ii));

for kk = (jj-1):-1:1
    %for cc = jj
        y(kk,ii) = y(kk,ii)-R(kk,jj)*X_est(jj,ii);
    %end
end
end
end

end

