function C = nuft_back_gridding(V0, px, py, Mr, E2x, E2y, E3x , E3y)
% Backward fast Gaussian gridding.
% Consider using the MEX version which is much faster.
%
% Copyright, Matthieu Guerquin-Kern, 2012

Nb  = numel(V0);
Msp = numel(E3x)-1;
C   = zeros(Mr);
px  = double(px);
py  = double(py);

%% First dummy version (Works)
% for i = 1:Nb
%     for l2=-Msp+1:Msp
%         indy = mod(py(i)+l2,Mr(2))+1;
%         Ey=V0(i)*E2y(i)^l2*E3y(abs(l2)+1);
%         for l1=-Msp+1:Msp
%             indx = mod(px(i)+l1,Mr(1))+1;
%             C(indx,indy)=C(indx,indy)+Ey*E2x(i)^l1*E3x(abs(l1)+1);
%         end
%     end
% end

%% Matlab spirit implementation (Works)

l2=-Msp+1:Msp;
l1=-Msp+1:Msp;
for i = 1:Nb
    Ey=(E2y(i).^l2)'.*E3y(abs(l2)+1);
    Ex=(E2x(i).^l1)'.*E3x(abs(l1)+1);
    C(mod(px(i)+l1,Mr(1))+1,mod(py(i)+l2,Mr(2))+1)=...
        C(mod(px(i)+l1,Mr(1))+1,mod(py(i)+l2,Mr(2))+1)+...
        V0(i)*(Ex*Ey');
end

%% Version close to MEX implementation (Works)
% for i = 1:Nb
%     E2ypl2 = E2y(i).^(1:Msp);
%     E2yml2 = E2y(i).^(-(0:Msp));
%     E2xpl1 = E2x(i).^(1:Msp);
%     E2xml1 = E2x(i).^(-(0:Msp));
%     % 1<= l2+1 < Msp and -Msp < -l2 <= 0
%     for l2=0:Msp-1
%         indyp = mod(py(i)+l2+1,Mr(2))+1;
%         indym = mod(py(i)-l2,Mr(2))+1;
%         Eyp=V0(i)*E2ypl2(l2+1)*E3y(abs(l2+1)+1);
%         Eym=V0(i)*E2yml2(l2+1)*E3y(abs(l2)+1);
%         % 1<= ll+1 < Msp and -Msp < -ll <= 0
%         for l1=0:Msp-1
%             indxp = mod(px(i)+l1+1,Mr(1))+1;
%             indxm = mod(px(i)-l1,Mr(1))+1;
%             C(indxp,indyp)=C(indxp,indyp)+Eyp*E2xpl1(l1+1)*E3x(abs(l1+1)+1);
%             C(indxm,indyp)=C(indxm,indyp)+Eyp*E2xml1(l1+1)*E3x(abs(l1)+1);
%             C(indxp,indym)=C(indxp,indym)+Eym*E2xpl1(l1+1)*E3x(abs(l1+1)+1);
%             C(indxm,indym)=C(indxm,indym)+Eym*E2xml1(l1+1)*E3x(abs(l1)+1);
%         end
%     end
% end