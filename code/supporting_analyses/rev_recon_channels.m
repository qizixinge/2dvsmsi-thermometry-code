function imageset = rev_recon_channels(rawdatahd,necho,nsa,bins,pff)
% Per-channel complex images with the k-space ordering of dataSOS.m/dataprocessing.m.
% Output: channel x Nx x Nyo x bins x 4 (effective-TE index). Bin and coil
% combination are left to the caller so that SOS and MIP share one reconstruction.
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
rawdatahd=permute(rawdatahd,[2 1 3]);
nex=floor(Ny/necho);
rawdatahd=reshape(rawdatahd,[channel,Nx,necho+3,bins,nex,nsa]);

if pff==1
    pe_steps=(1:(necho*nex))-0.5*necho*nex-1;
    PEorder = zeros(necho, nex);
    for k = 1:necho
        negStart = Ny/2+1 - nex/2*k;
        negEnd = negStart + nex/2-1;
        posStart = Ny/2+1 + nex/2*(k-1);
        posEnd = posStart + nex/2-1;
        PEorder(k, 1:floor(nex/2)) = pe_steps(negStart:negEnd);
        PEorder(k, (floor(nex/2)+1):nex) = pe_steps(posStart:posEnd);
    end
else
    btm=-floor(Ny-Nyo/2);
    pe_steps=btm:(Nyo/2-1);
    pe_steps=circshift(pe_steps,[0,-1]);
    [~,iPEmin1]=min(abs(pe_steps));
    k0curr1=floor((iPEmin1-1)/nex)+1;
    k0prescr1=1;
    PEorder=circshift(reshape(pe_steps,[nex,necho])',k0prescr1-k0curr1);
    neworder=[necho,2:(necho-1)];
    PEorder(2:end,:)=PEorder(neworder,:);
    PEorder(3:end,end)=PEorder([necho,3:(necho-1)],end);
end

PEordernatr=PEorder+Nyo/2+1;
ctrech=[1 2 3 4];
imageset=zeros(channel,Nx,Nyo,bins,ave);
for id=1:ave
    for bi=1:bins
        k5d=squeeze(rawdatahd(:,:,:,bi,:,:));
        temp0=k5d;
        temp0(:,:,ctrech,:, :) = [];
        ceh=ctrech(id);
        temp0=cat(3,k5d(:,:,ceh,:, :), temp0);
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
    end
end
end
