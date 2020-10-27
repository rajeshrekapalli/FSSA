function [X_R,Ev] = FSSA(X,L)
%UNTITLED Summary of this function goes here
% Inputs:
%         X----> input data vector
%         L---->windowlength, which is less than or equals to the half of the data points
%         
% Outputs:
%         X_R ----> Reconstructed series
%         Ev  ----> Eigen value percentages of Eigen modes
% This function allow the user to compute the Principal compoents
% de-noised signals etc., using Eigen spectrum and Weighted correlations
% Code testing:Run the following
% x=rand(100,1);
%[X_R,Ev] = FSSA(x,30);
% when prompted " Enter Eigen modes for signal reconstrution as [ 1:2 4
% 5]" enter the triplet group. More details on parameter selection are
% discussed in: Tiwari, R. K., & Rekapalli, R. (2020). Singular Spectrum Analysis with MATLAB®. 
% In Modern Singular Spectral-Based Denoising and Filtering Techniques for 2D and 3D Reflection Seismic Data 
% (pp. 125-138). Springer, Cham.

if nargin<1
    disp('No input data');
    pause();
    return;
end
X=X(:);
N=length(X);
N0=N;
if mod(N,2)>1
    X(nextpow2(n),1)=0;
    N=length(X);
end
disp(['Data contains' num2str(N) 'data points']);

if nargin<2
L=input('Enter window length (L)=');
end
X=fft(X);
K=N-L+1;
for ii=1:K
T(:,ii)=X(ii:ii+L-1,1);
end
[U S V] =svd(T);
Ev=diag(S)*100/sum(diag(S));
plot(1:L,Ev,'-rd');
title ('Eigen Spectrum', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('Eigen Triplet Number', 'FontName', 'Times New Roman', 'FontWeight', 'bold' ,'FontSize',12);
ylabel('Eiegn Value Percentage','FontName','Times New Roman', 'FontWeight', 'bold', 'FontSize',12);
set(gca,'FontName', 'TimeNewRoman','FontSize',12, 'FontWeight', 'bold');
disp('Press Enter to continue');
pause
for p=1:L
    T=U(:,p)*V(:,p)'*S(p,p);
    PC(:,p)=real(ifft(digavg(T,L,K,N)));
end
figure();
weightedcorr(PC);
clc
G = input(' Enter Eigen modes for signal reconstrution as [ 1:2 4 5]: ');
G=G(:);
X_R=sum(PC(1:N0,G),2);
end
function [X_R]=digavg(T,L,K,N)
X_R=zeros(N,1);  
Lp=min(L,K); 
Kp=max(L,K);
for s=0:Lp-2
for m=1:s+1
          X_R(s+1)=X_R(s+1)+(1/(s+1))*T(m,s-m+2);
end
end
 
for s=Lp-1:Kp-1
for m=1:Lp
         X_R(s+1)=X_R(s+1)+(1/(Lp))*T(m,s-m+2);
end
end
for s=Kp:N
for m=s-Kp+2:N-Kp+1
           X_R(s+1)=X_R(s+1)+(1/(N-s))*T(m,s-m+2);
end
end
end
function weightedcorr(PC)
%PC----> L number of Eigen modes computed from data (size of N×L)
N=size(PC,1);
L=size(PC,2);
K=N-L+1;
for m=1:N
if m<L+1
        w(m,1)=m;
elseif m>K+1
            w(m,1)=N-m;
else
        w(m,1)=L;
end
end
cc=zeros(L,L);
for m=1:L
for j=1:L
cc(m,j)=sum(w.*PC(:,m).*PC(:,j))/sqrt((sum(w.*PC(:,m).*PC(:,m))*sum(w.*PC(:,j).*PC(:,j))));
end
end
imagesc(cc);
title ('Weighted Correlation', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('Eigen Triplet', 'FontName', 'Times New Roman', 'FontWeight', 'bold' ,'FontSize',12);
ylabel('Eiegn Triplet ','FontName','Times New Roman', 'FontWeight', 'bold', 'FontSize',12);
set(gca,'FontName', 'TimeNewRoman','FontSize',12, 'FontWeight', 'bold');
set(gca, 'xdir','normal'); 
set(gca, 'ydir','normal'); 
colormap('gray');
caxis([0 1]);
colorbar;
end

