%% Author : Goli Srikanth
% Under the guidance of  Prof. Vijay Kumar Chakka
% Mail ID: gs499@snu.edu.in
% Ph.D research Scholar
% Shiv Nadar University
% Project : Implementation of Ramanujan Transform
%%
function [nonOrthoBasesMtrx,orthoBasesMtrx,ramS,neOrthoBasesMtrx] = nonOrtho_ramanujanBaseMtrx(N)
nonOrthoBasesMtrx = zeros(N);
orthoBasesMtrx = zeros(N);
neOrthoBasesMtrx = zeros(N);
a = [];
index = 1;
count = 1;
for q = 1:N
    ramSqnce = zeros(1,q);
% Generating Ramanujan Sequency
    for n = 1:q
       for k = 1:q
          if gcd(k,q) == 1% here I am making very low values to zeros by using int8 
              %then again to double using double keyword, because airthematic
              %operations will not work on integer values
              ramSqnce(n) = ramSqnce(n) + exp(1i*2*pi*k*(n-1)/q);
          end  
       end
       %ramSqnce = double(int8(ramSqnce));
    end
    ramSqnce = double(int64(ramSqnce));
    tempMtrx = zeros(N,1);
    tempMtrx1 = [];
    if mod(N,q) == 0
        ramS(count).q = q;
        ramS(count).bases = ramSqnce;
        for i = 1:N
            if i>q
                tempMtrx(i,1) = ramSqnce(mod(i-1,q)+1);
            else
                tempMtrx(i,1) = ramSqnce(i);
            end
        end
       a = [a tempMtrx.'*tempMtrx];
       ramS(count).enerGy = a;
       tempMtrx2 = tempMtrx.*(1/sqrt(a(count))); %half of the energy per matrix
        for j = 1:tempMtrx(1)
            orthoBasesMtrx(:,index) = circshift(tempMtrx,j-1);
            neOrthoBasesMtrx(:,index) = circshift(tempMtrx2,j-1);
            [tempMtrx1] = [tempMtrx1 orthoBasesMtrx(:,index)];
            index = index+1;
        end
        ramS(count).gQ = tempMtrx1;
        count = count+1;
    end
    
    for i = 1:N
        if i>q
            nonOrthoBasesMtrx(i,q) = ramSqnce(mod(i-1,q)+1);
        else
            nonOrthoBasesMtrx(i,q) = ramSqnce(i); 
        end
    end
   
end
nonOrthoBasesMtrx = nonOrthoBasesMtrx;
orthoBasesMtrx = orthoBasesMtrx;
end
