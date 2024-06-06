function [X_est] = SIC_Equalizer(y,R)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[Nr Nt]=size(R);
[xx yy]=size(y);
for ii = yy:1
for jj = Nr:1
X(jj,ii)= y(jj,ii)/R(jj,jj);

X_est(jj,ii) = QAM16_slicer(X(jj,ii));

for kk = (jj+1):Nr
    %for cc = jj
        y(kk,ii) = y(kk,ii)-R(kk,jj)*X_est(jj,ii);
    %end
end
end
end

end

