function [expPairSum] = expPairSumCnstrcn(q)
% Author : Goli Srikanth Email ID: gs499@snu.edu.in, Research
% Scholar, Shiv Nadar University. Note: If there are any
% modifications/suggestions, please contact above mail ID.
n = 1:q; gcdIndex = n(gcd(n,q)==1);
expPairs = {};
expPairSum = struct();
Count = 1;
if q==1
    expPairSum(Count).Sum = 1;
    expPairSum(Count).Pairs = 1;
    expPairSum(Count).circMtrx = 1;
    expPairSum(Count).Rank = 1;
    Count = Count+1;
elseif q==2
    expPairSum(Count).Sum = [1;-1];
    expPairSum(Count).Pairs = 1;
    expPairSum(Count).circMtrx = [1 -1;-1 1];
    expPairSum(Count).Rank = 1;
    Count = Count+1;
else
    if (length(gcdIndex)/2)>=1
        for n = 1:length(gcdIndex)/2
            expPairs{n} = [gcdIndex(1,n) gcdIndex(1,end-n+1)];
        end
    else
        expPairs{1} = gcdIndex;
    end
    for k = 1:length(expPairs)
        Temp = []; Temp1=[];circMtrx=[];circMtrx1=[];%circMtrx2=[];circMtrx3=[];
        for n = 0:q-1
            [Temp] = [Temp;exp((1j*2*pi*expPairs{1,k}(1)*n)/q)+exp((1j*2*pi*expPairs{1,k}(2)*n)/q)];
            [Temp1] = [Temp1;(exp((1j*2*pi*expPairs{1,k}(1)*n)/q)-exp((1j*2*pi*expPairs{1,k}(2)*n)/q))/(1j)];
        end
        expPairSum(Count).Sum = Temp;
        expPairSum(Count).Sum1 = Temp1;
%         expPairSum(Count).normSum = Temp/sqrt(norm(Temp));
%         expPairSum(Count).normSum1 = Temp1/sqrt(norm(Temp1));
        expPairSum(Count).Pairs = expPairs{1,k};
        for n = 0:q-1
            [circMtrx] = [circMtrx circshift(Temp,n)];
            [circMtrx1] = [circMtrx1 circshift(Temp1,n)];
%             [circMtrx2] = [circMtrx2 circshift(expPairSum(Count).normSum,n)];
%             [circMtrx3] = [circMtrx3 circshift(expPairSum(Count).normSum1,n)];
        end
        expPairSum(Count).circMtrx = circMtrx;
        expPairSum(Count).Rank = rank(circMtrx);
        expPairSum(Count).circMtrx1 = circMtrx1;
        expPairSum(Count).Rank1 = rank(circMtrx1);
%         expPairSum(Count).normCircMtrx = circMtrx2;
%         expPairSum(Count).normCircMtrxRank = rank(circMtrx2);
%         expPairSum(Count).normCircMtrx1 = circMtrx3;
%         expPairSum(Count).normCircMtrxRank1 = rank(circMtrx3);
        Count = Count+1;
    end
end
end
