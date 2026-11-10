function [img4is,img4isoss,T2MAPs,T2MAPsoss] = MAP(necho,bins,pff,file,filt)
% Read one raw acquisition, reconstruct MIP/SOS images, and fit native-grid apparent T2.
%filt=1 when mean filter, filt=0 when no mean filter

nsa=1;
twxh054=mapVBVD(file);%for example 'meas_MID00613_FID02871_pulseq.dat'
rawdatah054=twxh054{end}.image();
img4i = dataprocessing(rawdatah054,necho,nsa,bins,pff);
img4isos=dataSOS(rawdatah054,necho,nsa,bins,pff);
h=ones(3,3)/9;
% if filt==1
%     for idx=1:4
%         img4i(:, :, idx)=imfilter(img4i(:, :, idx), h, 'same', 'replicate');
%         img4isos(:, :, idx)=imfilter(img4isos(:, :, idx), h, 'same', 'replicate');
%     end
% end





% roimg4i = img4i .* double(B);
roimg4i = img4i;


% T2MAPs=RDNLS(roimg4i,[1 2 3 4],8,100);%T2MAP(:,:,8,2)
% T2MAPsoss=RDNLS(img4isos,[1 2 3 4],8,100);%T2MAPsos(:,:,8,2)
roimg4i(:,:,1)=[];
img4isos(:,:,1)=[];


fit_timer = tic;
T2MAPs=RDNLS(roimg4i,[1 2 3],8,60);%T2MAP(:,:,8,2)
fit_mip_seconds = toc(fit_timer);
fit_timer = tic;
T2MAPsoss=RDNLS(img4isos,[1 2 3],8,60);%T2MAPsos(:,:,8,2)
fit_sos_seconds = toc(fit_timer);
fprintf('RDNLS_TIMING pff=%.1f MIP=%.3f s SOS=%.3f s (100x100, three echoes)\n', ...
    pff,fit_mip_seconds,fit_sos_seconds);
img4is=roimg4i(:,:,1);
img4isoss=img4isos(:,:,1);
% for TE=32ms
% img4is=roimg4i(:,:,3);
% img4isoss=img4isos(:,:,3);


if filt==1
        T2MAPs=imfilter(T2MAPs, h, 'same', 'replicate');
        T2MAPsoss=imfilter(T2MAPsoss, h, 'same', 'replicate');

end

end
