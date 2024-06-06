function [Summation] = matrixCnstrcn(Summation)
% Author : Goli Srikanth Email ID: gs499@snu.edu.in, Research
% Scholar, Shiv Nadar University Note: If there are any
% modifications/suggestions, please contact above mail ID 
% Aim: Dictionary construction
trMatrix = [];normTrMatrix = [];
Div = divisors(Summation.q);
Summation.Divisors = Div;
for i = 1:length(Div)
    [expPairSum] = expPairSumCnstrcn(Div(i));
    if Div(i)==1 || Div(i)==2
        Rank = 1;
    else
        Rank = 2;
    end
    n = 1:Div(i); totientVec = n(gcd(n,Div(i))==1);
    Summation.Div(i).expPairSum = expPairSum;
    Summation.Div(i).Divisor = Div(i);
    Summation.Div(i).totientVec = totientVec;
    for j = 1:length(expPairSum)
        Temp = [];Temp1 = [];Temp2 = [];Temp3 = [];subMtrx1 = [];
        subMtrx = [];
        [Temp] = [Temp repmat(expPairSum(j).Sum,Summation.q/Div(i),1)];
        [Temp2] = [Temp2 repmat((expPairSum(j).Sum)./sqrt(Temp.'*Temp),Summation.q/Div(i),1)];
        if Rank ==2
            [Temp1] = [Temp1 repmat(expPairSum(j).Sum1,Summation.q/Div(i),1)];
            [subMtrx] = [subMtrx Temp Temp1];
            [Temp3] = [Temp3 repmat((expPairSum(j).Sum1)./sqrt(Temp1.'*Temp1),Summation.q/Div(i),1)];
            [subMtrx1] = [subMtrx1 Temp2 Temp3];
        else
            [subMtrx] = [subMtrx Temp];
            [subMtrx1] = [subMtrx1 Temp2];
        end
        %         for k = 1:Rank
        %             [subMtrx] = [subMtrx circshift(Temp,k-1)];
        %         end
        [trMatrix] = [trMatrix subMtrx];
        [normTrMatrix] = [normTrMatrix subMtrx1];
    end
end
Summation.transformMatrix = trMatrix;
Summation.normTransformMatrix = normTrMatrix;
end
