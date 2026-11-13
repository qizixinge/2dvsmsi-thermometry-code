function [T2map] = RDNLS(roimg4i,ctrech,tesp,T2approx)
% Purpose: voxelwise mono-exponential fitting of multi-echo magnitudes using a 0-400 ms grid with 0.05 ms steps.
% No intensity threshold is applied here; display masks are applied in one place
% (statistics_analysis_final.m, DISPLAY_FRAC), identically for 2DVSMSI and GRE.

%img4i default n*n*4
%T2approx need to be larger than T2real
% h=ones(3,3)/9;
% roimg4i=imfilter(roimg4i,h,'symmetric');
A=size(roimg4i,1);
B=size(roimg4i,2);
C=size(roimg4i,3);
nbr=A*B;
rimg4i=reshape(roimg4i,nbr,C);
rimg4i=rimg4i';
t2srch=0:0.05:400;
D=size(t2srch,2);
t2srch=reshape(t2srch,1,1,D);




%weight=exp(2*(ctrech-ctrech(1))*tesp./t2srch);
weight=exp(2*(ctrech-ctrech(1))*tesp/T2approx);
weight=reshape(weight,C,1);
weight=1;
n=ctrech;n=repmat(n',1,nbr);
J=squeeze(sum(weight.*rimg4i.^2,1)-(sum(weight.*rimg4i.*exp(-n*tesp./t2srch),1)).^2./sum(weight.*exp(-2*n*tesp./t2srch),1));


J=J';
[~,idT2]=min(J);
optT2=t2srch(idT2);


T2map=reshape(optT2,A,B);


end
