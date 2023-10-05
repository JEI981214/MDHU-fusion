close all;
clear all;
clc;

B=imread('37-2SPECT-011.png');
A=imread('38-1MRI-015.png');

if size(A,3)>1
    A=rgb2gray(A);             
end

A =im2double(A);
B =im2double(B);

B_YUV=ConvertRGBtoYUV(B);
B_Y=B_YUV(:,:,1);
[hei, wid] = size(A);  

%% decomposition
smallNum = 1e-3;
thr = 1; 
lambda = 0.5;  % 0.5 for "01", 1.25 for "02", 0.75 for "03", 0.7 for "04"
rData = 1;
rSmooth = 1;
aData = smallNum;
aSmooth = smallNum;
bData = thr; 
bSmooth = thr;
alpha = 0.5;
stride = 1;
iterNum = 10;
AL = generalized_smooth(A, A, lambda, rData, rSmooth, aData, bData, aSmooth, bSmooth, alpha, stride, iterNum);
BL = generalized_smooth(B_Y, B_Y, lambda, rData, rSmooth, aData, bData, aSmooth, bSmooth, alpha, stride, iterNum);
AH=A-AL;
BH=B_Y-BL;

%% ¸ßÆµ
lambda=0.01; 
flag=2; 
D=cell2mat(struct2cell(load('D1.mat')));  %
FH1=CSR(AH, BH, D, lambda, flag);
FH1=double(FH1);
Diff1=FH1-AH;

%
D=cell2mat(struct2cell(load('D2.mat')));  % 
FH2=CSR(AH, BH, D, lambda, flag);
FH2=double(FH2);
Diff2=FH2-AH;
% % 
% % %
D=cell2mat(struct2cell(load('D3.mat')));  %
FH3=CSR(AH, BH, D, lambda, flag);
FH3=double(FH3);
Diff3=FH3-AH;
% 
% 
% % %
D=cell2mat(struct2cell(load('D4.mat')));  %
FH4=CSR(AH, BH, D, lambda, flag);
FH4=double(FH4);
Diff4=FH4-AH;

% ÅÐ¶Ï
[row,column]=size(AH);
T=0.05;%*************************************************************************************
% 
Diff1(abs(Diff1)<T)=0;
Diff1(abs(Diff1)>=T) = 1;
Diff2(abs(Diff2)<T)=0;
Diff2(abs(Diff2)>=T) = 1;
Diff3(abs(Diff3)<T)=0;
Diff3(abs(Diff3)>=T) = 1;
Diff4(abs(Diff4)<T)=0;
Diff4(abs(Diff4)>=T) = 1;
C= Diff1+Diff2+Diff3+Diff4;
for i=1:row
for j=1:column
if C(i,j)>4/2
FH(i,j)=BH(i,j);
else
FH(i,j)=AH(i,j);
end
end
end
          
%% µÍÆµ
ww=3;
LE1=local_energy(AL,ww);
LE2=local_energy(BL,ww);
%%
s1_or=AL;
s2_or=BL;
[w,h,c] = size(s1_or);
unit = sqrt(w);

s1 = s1_or;
s2 = s2_or;

% norm
s1_nu = zeros(1,h);
s2_nu = zeros(1,h);
for i=1:h
    temp1 = reshape(s1(:,i), [unit unit]);
    [U, S, V] = svd(temp1, 'econ');
    nu_norm1 = sum(diag(S));
    s1_nu(:,i)=nu_norm1;
    temp2 = reshape(s2(:,i), [unit unit]);
    [U, S, V] = svd(temp2, 'econ');
    nu_norm2 = sum(diag(S));
    s2_nu(:,i)=nu_norm2;
end
s1_norm = s1_nu;
s2_norm = s2_nu;
w1_value = s1_norm./(s1_norm+s2_norm ); 
w2_value = s2_norm./(s1_norm+s2_norm );

w1 = repmat(w1_value, w, 1);
w2 = repmat(w2_value, w, 1);
%%
[row,column]=size(A);
TT=1 %************************************************************************************
for i=1:row
    for  j=1:column
        if LE1(i,j)-LE2(i,j)>TT || (abs(LE1(i,j)-LE2(i,j))<TT  &&  w1(i,j)>w2(i,j));
            FL(i,j)=AL(i,j);
        else LE2(i,j)-LE1(i,j)>TT || (abs(LE2(i,j)-LE1(i,j))<TT  &&  w2(i,j)>w1(i,j));
            FL(i,j)=BL(i,j);
        end
    end
end

%% ÈÚºÏ
F=FL+FH;

F=double(F);
imgf_YUV=zeros(hei,wid,3);                                                 
imgf_YUV(:,:,1)=F;
imgf_YUV(:,:,2)=B_YUV(:,:,2);
imgf_YUV(:,:,3)=B_YUV(:,:,3);
FINAL=ConvertYUVtoRGB(imgf_YUV);

figure,imshow(FINAL);

