function [T2map] = RDNLS(roimg4i,ctrech,tesp,T2approx)
% Purpose: voxelwise mono-exponential fitting of multi-echo magnitudes using a 0-400 ms grid with 0.05 ms steps; zero voxels below three times background noise.

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


% Extract 10 rows near the top, bottom, left, and right edges
top_rows = roimg4i(1:10, :,1);
bottom_rows = roimg4i(A-9:A, :,1);
left_rows = roimg4i(:, 1:10,1);
right_rows = roimg4i(:, B-9:B,1);
% Combine the extracted edge rows into a vector
edge_rows = [top_rows(:); bottom_rows(:); left_rows(:); right_rows(:)];
% Calculate the max
max_intensity = max(edge_rows,[],'all');%median(edge_rows);
% Calculate background noise intensity
noise = max_intensity;

T2map=reshape(optT2,A,B);

T2map(roimg4i(:,:,1)<3*noise)=0;

end
