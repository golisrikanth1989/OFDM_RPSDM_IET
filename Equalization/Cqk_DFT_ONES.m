function [ A] = Cqk_DFT_ONES( N )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
rr = 1;
q = divisors(N);
 A = zeros(length(q),N);
for qq = q
a = [];
for kk = 0:qq-1
    
if(gcd(kk,qq)==1)
 a = [a kk];
 
 end
end
a = (N/qq)*a;
 a = a+1;
 A(rr,a) = 1/length(a);
 rr = rr+1;
 end


end

