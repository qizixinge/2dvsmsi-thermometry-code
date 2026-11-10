function [MNS,SDS,SDT,alpha,Dsroi,alpham] = PRFreconnew(phasemapROI,ROI1,ROI2,ROI3,ROI4,Tstep,oil1,oil2,oil3,oil4,reft)
% Purpose: PRF thermometry with phase unwrapping, four-oil-reference drift correction, and probe-temperature calibration in four 3x3 ROIs.
% checked at 20260815 by DSpro+me
D=round(1);
DD=2*D+1;
sz=128;
oil=cat(3,oil1,oil2,oil3,oil4);
mask=max(oil,[],3);mask=double(mask);
ROI=cat(3,ROI1,ROI2,ROI3,ROI4); % numx2x4
num=size(phasemapROI,3);
% reft: time points whose per-voxel mean is the phase reference (default: all).
% The reference is a per-ROI constant after the oil-drift correction, so the
% calibration slopes and LOOCV do not depend on it.
if nargin < 11 || isempty(reft), reft = 1:num; end
phasemaproi=zeros(DD,DD,num,4);
gp=gre_params();lamda=gp.lam;% B0*gamma*dTE, dTE=2.46 ms (phase-difference image of the dual-echo GRE; see gre_params.m)

id=find(mask);
phasemapROI=permute(phasemapROI,[3 1 2]);
phasemap=permute(unwrap(phasemapROI),[2 3 1]);
phasemapROIref=mean(phasemap(:,:,reft),3);

dphasemap=phasemap-phasemapROIref;
poils=zeros(1,1,num);
for idx=1:num
    temp=dphasemap(:,:,idx);
    poils(idx)=mean(temp(id),'all');
end
dphasemap=dphasemap-poils;

phasemap=dphasemap+phasemapROIref;

for id4=1:4
    for idx=1:num%especially for heating
        phasemaproi(:,:,idx,id4)=phasemap((ROI(idx,1,id4)-D):(ROI(idx,1,id4)+D),(ROI(idx,2,id4)-D):(ROI(idx,2,id4)+D),idx);
    end
end

mns=zeros((DD-2),(DD-2),num,4);
alpha=zeros((DD-2),(DD-2),4);
an=mean(Tstep, 1);%sum(Tstep)/num;

for idx=1:(DD-2)
    for idy=1:(DD-2)
        xyphaseroi=phasemaproi(idy:(idy+2),idx:(idx+2),:,:);
        mns(idy,idx,:,:)=mean(xyphaseroi,[1 2]);%MNS(idy,idx,:,:)=mean(xyT2MAProi,[1 2]);
        tmns=squeeze(mns(idy,idx,:,:));%num*4
        amns = mean(tmns, 1); % 1×4 vector

        dx = Tstep - repmat(an, num, 1); % num×4 matrix
        dmns = tmns - repmat(amns, num, 1); % num×4 matrix
        % calculate the slope
        alpha(idy,idx,:) = sum(dx .* dmns, 1) ./ (lamda*sum(dx.^2, 1)); % 1×4 vector
    end
end



alpham=reshape(mean(alpha,[1 2]),1,4);
%alphastd=reshape(std(alpha,0,[1 2]),1,4);


sroi=(phasemaproi(round((DD+1)/2-1):round((DD+1)/2+1),round((DD+1)/2-1):round((DD+1)/2+1),:,:))./(lamda);%3*3*num*4%(phasemaproi(D:(D+2),D:(D+2),:,:))./(lamda);%
Amns=reshape(mean(sroi(:,:,reft,:),3),3,3,4);% 3x3×4 ;
Dsroi=sroi-reshape(Amns,3,3,1,4);%/lamda;%3*3*8*4
MNS=squeeze(mean(Dsroi,[1 2]));%(mns-jieju)/lamda;%num*4

        dtroi=MNS./reshape(alpha(1,1,:),1,4);%(tmns-jieju)./alpham;% num×4
        SDT=std(dtroi-(Tstep-an),1,1);%sqrt(sum((dtroi-(Tstep-an)).^2,1));
SDS=(squeeze(std(Dsroi,0,[1 2])));
%SDS=0;

end
