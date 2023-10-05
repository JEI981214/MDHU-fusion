function imgf = CSR(s1_h, s2_h, D, lambda, flag)

pt = [];
opt.Verbose = 0;
opt.MaxMainIter = 500;
opt.rho = 100*lambda + 0.5;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
[X1, optinf1] = cbpdn(D, s1_h, lambda, opt);
[X2, optinf2] = cbpdn(D, s2_h, lambda, opt);

%activity level measure
A1=sum(abs(X1),3);
A2=sum(abs(X2),3);

if flag == 1  
    r=9;  
else
    r=3; 
end

ker=ones(2*r+1,2*r+1)/((2*r+1)*(2*r+1));
AA1=imfilter(A1,ker);
AA2=imfilter(A2,ker);
decisionMap=AA1>AA2;



%detail layer fusion
[height,width]=size(A1);
X=X1;
for j=1:height
    for i=1:width
        if decisionMap(j,i)==0
            X(j,i,:)=X2(j,i,:);
        end
    end
end
imgf = ifft2(sum(bsxfun(@times, fft2(D, size(X,1), size(X,2)), fft2(X)),3),'symmetric');
