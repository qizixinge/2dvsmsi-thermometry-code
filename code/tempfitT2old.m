function [dMNS,SDS,SDT,M2,DT2MAProi,M2m,M2mperc,T2base] = tempfitT2old(T2MAP,ROI1,ROI2,ROI3,ROI4,tempr,reft)
% Purpose: T2 calibration in four 3x3 ROIs by averaging nine voxels before regression against probe temperature.
% D=round(2);
% DD=2*D+1;
%
%
% ROI=cat(3,ROI1,ROI2,ROI3,ROI4);
% num=size(T2MAP,3);
% T2MAProi=zeros(DD,DD,num,4);
%
% for id4=1:4
%     for idx=1:num%especially for heating
%         T2MAProi(:,:,idx,id4)=T2MAP((ROI(idx,1,id4)-D):(ROI(idx,1,id4)+D),(ROI(idx,2,id4)-D):(ROI(idx,2,id4)+D),idx);
%         %MNS(idx)=;
%     end
% end
% MNS=zeros((DD-2),(DD-2),num,4);
% M2=zeros((DD-2),(DD-2),4);
%
%
% an = mean(tempr, 1); % 1×4
%
% for idx=1:(DD-2)
%     for idy=1:(DD-2)
%         xyT2MAProi=T2MAProi(idy:(idy+2),idx:(idx+2),:,:);
%         MNS(idy,idx,:,:)=mean(xyT2MAProi,[1 2]);
%         tMNS=squeeze(MNS(idy,idx,:,:));%num*4
%
%         aMNS = mean(tMNS, 1); % 1×4
%
%         dx = tempr - repmat(an, num, 1); % num×4
%         dMNS = tMNS - repmat(aMNS, num, 1); % num×4
%
%         M2(idy,idx,:) = sum(dx .* dMNS, 1) ./ sum(dx.^2, 1); % 1×4
%     end
% end
%
% tMNS=squeeze(mean(MNS,[1 2]));%num*4 %squeeze(MNS(D,D,:,:));
% AMNS=reshape(mean(tMNS,1),1,4);% 1×4
%
% M2m=reshape(mean(M2,[1 2]),1,4);
% M2mstd=reshape(std(M2,0,[1 2]),1,4);
% M2mperc=M2m./AMNS*100;
% M2stdperc=M2mstd./AMNS*100;
% T2base=AMNS;
%
% DT2MAProi=MNS-reshape(mean(MNS,3),3,3,1,4);
% tMNS=squeeze(mean(DT2MAProi,[1 2]));
% SDS=squeeze(std(DT2MAProi,0,[1 2]));
%
%
% cMNS=squeeze(MNS(round((DD-1)/2),round((DD-1)/2),:,:));% num×4
%         tempm2=(squeeze(M2(round((DD-1)/2),round((DD-1)/2),:)))';
%         jieju = mean(cMNS,1) - tempm2 .* an; % 1×4
%         troi=(cMNS-jieju)./tempm2;% num×4
% %SDT=sqrt(sum((troi-tempr).^2,1));
% SDT=std(troi-tempr,1,1);

D = round(1);
DD = 2*D + 1;  % DD = 3

ROI = cat(3, ROI1, ROI2, ROI3, ROI4);
num = size(T2MAP, 3);
% reft: time points whose per-voxel mean is the reference of dMNS / DT2MAProi / SDS
% (default: all).  The reference is a per-ROI constant, so slopes and LOOCV do not depend on it.
if nargin < 7 || isempty(reft), reft = 1:num; end
T2MAProi = zeros(DD, DD, num, 4);

% 3×3 voxels for each ROI
for id4 = 1:4
    for idx = 1:num
        T2MAProi(:,:,idx,id4) = T2MAP((ROI(idx,1,id4)-D):(ROI(idx,1,id4)+D), ...
                                       (ROI(idx,2,id4)-D):(ROI(idx,2,id4)+D), idx);
    end
end

MNS = zeros((DD-2), (DD-2), num, 4);  %  1×1×num×4
M2 = zeros((DD-2), (DD-2), 4);        %  1×1×4

an = mean(tempr, 1); % 1×4 vector

for idx = 1:(DD-2)  % （idx=1）
    for idy = 1:(DD-2)  % （idy=1）
        xyT2MAProi = T2MAProi(idy:(idy+2), idx:(idx+2), :, :);  % 3×3×num×4
        MNS(idy, idx, :, :) = mean(xyT2MAProi, [1 2]);
        tMNS = squeeze(MNS(idy, idx, :, :));  % num×4

        aMNS = mean(tMNS, 1); % 1×4 vector

        % calculate the delta
        dx = tempr - repmat(an, num, 1);       % num×4 matrix
        dMNS = tMNS - repmat(aMNS, num, 1);    % num×4 matrix

        % calculate the slope M2
        M2(idy, idx, :) = sum(dx .* dMNS, 1) ./ sum(dx.^2, 1); % 1×4 vector
    end
end

tMNS = squeeze(mean(MNS, [1 2]));  % num×4
AMNS = reshape(mean(tMNS, 1), 1, 4);  % 1×4 vector

M2m = reshape(mean(M2, [1 2]), 1, 4);
%M2mstd = reshape(std(M2, 0, [1 2]), 1, 4);
M2mperc = M2m ./ AMNS * 100;
%M2stdperc = M2mstd ./ AMNS * 100;
T2base = AMNS;

%DT2MAProi = MNS - reshape(mean(MNS, 3), 1, 1, 1, 4);
DT2MAProi = T2MAProi - reshape(mean(T2MAProi(:,:,reft,:), 3), 3, 3, 1, 4);
dMNS = squeeze(mean(DT2MAProi, [1 2]));
SDS = squeeze(std(DT2MAProi, 0, [1 2]));
%SDS = squeeze(std(T2MAProi, 0, [1 2]));

cMNS = squeeze(MNS(round((DD-1)/2), round((DD-1)/2), :, :)); % num×4， (DD-1)/2 = 1
% calculate the intercept
tempm2 = (squeeze(M2(round((DD-1)/2), round((DD-1)/2), :)))'; % 1×4
jieju = mean(cMNS, 1) - tempm2 .* an; % 1×4 vector
troi = (cMNS - jieju) ./ tempm2;      % num×4 matrix
SDT = std(troi - tempr, 1, 1);


end
