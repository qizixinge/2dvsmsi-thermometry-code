function [img4i] = dataprocessing(rawdatahd,necho,nsa,bins,pff)
% Reorder raw multi-bin echoes, reconstruct partial Fourier, and combine bins by maximum intensity projection.

%assume Nx=Ny

Nxy=size(rawdatahd,1);
Nx=Nxy;iter=30;ave=4;
channel=size(rawdatahd,2);
if pff==1

    Ny=Nxy;
    Nyo=Nxy;
else

    Nyo=Nxy;
    Ny=round(Nyo*pff);
end
rawdatahd=permute(rawdatahd,[2 1 3]);%channel, Nx, Ny
%Nx=100; Ny=100;
nex=floor(Ny/necho);
%tesp=8;
%channel=20;
%bins=12;
%T2approx=100;
rawdatahd=reshape(rawdatahd,[channel,Nx,necho+3,bins,nex,nsa]);

if pff==1%without pocs, no problem

    pe_steps=(1:(necho*nex))-0.5*necho*nex-1;
    PEorder = zeros(necho, nex);
    for k = 1:necho
        % Calculate the start and end indices of the negative part.
        negStart = Ny/2+1 - nex/2*k;
        negEnd = negStart + nex/2-1;

        % Calculate the start and end indices of the positive part.
        posStart = Ny/2+1 + nex/2*(k-1);%0 position
        posEnd = posStart + nex/2-1;

        % fill the negative part
        PEorder(k, 1:floor(nex/2)) = pe_steps(negStart:negEnd);

        % fill the positive part
        PEorder(k, (floor(nex/2)+1):nex) = pe_steps(posStart:posEnd);
    end
else% with pocs
    btm=-floor(Ny-Nyo/2);%
    pe_steps=btm:(Nyo/2-1);
    % if 0==mod(necho,2)
    %     pe_steps=circshift(pe_steps,[0,-round(nex/2)]); % for odd number of echoes we have to apply a shift to avoid a contrast jump at k=0
    % end
    pe_steps=circshift(pe_steps,[0,-1]);%floor(nex/2)]);%Not universally applicable
    [~,iPEmin1]=min(abs(pe_steps));
    k0curr1=floor((iPEmin1-1)/nex)+1; % calculate the 'native' central echo index
    k0prescr1=1;%max(round(TEeff/TE1),1); % echo to be aligned to the k-space center
    PEorder=circshift(reshape(pe_steps,[nex,necho])',k0prescr1-k0curr1);
    neworder=[necho,2:(necho-1)];
    PEorder(2:end,:)=PEorder(neworder,:);%Not universally applicable
    PEorder(3:end,end)=PEorder([necho,3:(necho-1)],end);%Not universally applicable
    %PEorder(2:end,ceil(nex/2):end)=PEorder(neworder,ceil(nex/2):end);


end



PEordernatr=PEorder+Nyo/2+1;% id of Ny
ctrech=[1 2 3 4];
imageset=zeros(channel,Nx,Nyo,bins,ave);
sliceid=[5,3,1,11,9,7,4,2,0,10,6,8];   % bin acquisition order of the sequence writers (unused here: MIP/SOS do not depend on bin order)

for id=1:ave
    for bi=1:bins
        %for NSA=1:nsa
            k5d=squeeze(rawdatahd(:,:,:,bi,:,:));
            %temp0=rawdatahd;
            temp0=k5d;
            temp0(:,:,ctrech,:, :) = [];
            ceh=ctrech(id);

            temp0=cat(3,k5d(:,:,ceh,:, :), temp0);
            %temp0=cat(3,k4d(:,:,ceh, :), temp0*exp(-(ceh-1)*tesp/T2approx));
            %Amplitude compensation
            temp1=zeros(channel,Nx,Nyo,nsa);

            for ne=1:necho
                for nx=1:nex
                    temp1(:,:,PEordernatr(ne,nx),:)=reshape(temp0(:,:,ne,nx,:),[channel,Nx,nsa]);
                end
            end
            temp2=sum(temp1,4)/nsa;
            if pff==1
                imageset(:,:,:,bi,id)=reshape(ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(temp2, 3),2), [],3),[],2), 3),2),channel,Nx,Nyo,1,1);
            else
                imageset(:,:,:,bi,id)=reshape(pocsnew( temp2, iter, false ),channel,Nx,Nyo,1,1);
            end
        %end
    end
end
    imagemag=squeeze(sqrt(sum( abs(imageset).^2, 1)));
    img4=squeeze(max(imagemag,[],3));

% Images are returned at the NATIVE acquired matrix (Nx x Nyo, here 100x100).
% Interpolation to the 128x128 display grid was moved downstream (see MAP.m /
% statistics_analysis_final.m) so that T2 fitting, ROI placement and
% calibration are all performed on the native voxel grid before interpolation.
% Interpolating before fitting correlates neighbouring voxels and invalidates
% the 3x3 ROI error bars; do not reinstate an interp2 call here.
img4i=img4;





end
