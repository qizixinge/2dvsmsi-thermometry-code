% Export the five matched-field-of-view Figure 3 panels as PNG and PDF.
%% paper_figure.m
% Figure 3 of the manuscript, as FIVE separate image files so the LaTeX
% subfigure block (\includegraphics{MIP.png} ... ) picks them up unchanged.
%
%   A  MIP.png     2DVSMSI-PFA, MIP  -- two large yellow boxes (artefact zones)
%   B  SOS.png     2DVSMSI-PFA, SOS  -- same two boxes, for A/B comparison
%   C  GRE.png     GRE (128 grid)    -- ROI I-IV
%   D  BVSMSI.png  2DVSMSI-PFA, SOS  -- ROI I-IV
%   E  AVSMSI.png  2DVSMSI-Full, SOS -- ROI I-IV
%
% FRAMES: POCS and Full now share the standard centered-FFT origin.
% ROI overlays follow the actual statistical coordinates from analysis_geometry.
% PFA indices are transformed +1 native column with the corrected image origin;
% Full indices stay unchanged. The archived ROI file is not overwritten.
% Rod size from the physical geometry (8 cm x 1 cm) converted per grid:
%   2DVSMSI 100 grid, 2.56 mm/voxel : 31 x 4 voxels -> rows 40-70
%   GRE     128 grid, 2.00 mm/voxel : 40 x 5 voxels -> rows 51-91, centred on
%   the measured dark-bore axis (row 71, col 64.5).  The dark bore extends past
%   the metal (rows 46-96) because of the susceptibility artefact; that is NOT
%   the rod.
%
% Writes five PNG/PDF pairs to project.figures (output/figures); the TBME tex
% uses the PNG files.

clear; close all;
project = project_paths();
cd(project.code);

OUT_DIR = project.figures;
TP = 3;     % time point of the illustrated images (= TP_CMP of statistics_analysis_final.m)

S1 = load(fullfile(project.data,'T2MAP_and_img4i_native100.mat'), 'img4i','img4isos','reconstruction_frame');
R  = load(fullfile(project.data,'ROI_and_tempr_resolved.mat'));
geometry = analysis_geometry(R,S1.reconstruction_frame);
C2D = geometry.pfa_roi;
C2D_full = geometry.full_roi;
%     % 4x2 [row col], PFA frame, native 100
CPRF = geometry.prf_roi;              % 4x2 [row col], GRE 128 grid

imgMIP_pfa  = S1.img4i(:,:,TP,1);     % native 100x100, time point TP
imgSOS_pfa  = S1.img4isos(:,:,TP,1);
imgSOS_full = S1.img4isos(:,:,TP,3);
G = load(fullfile(project.data,'phaseMAPr_and_magmap.mat'),'magmap');
magmap = double(G.magmap);   % series 0016 (acquisition 2): the only GRE magnitude exported; no acquisition-3 GRE magnitude exists

ROD_PFA  = geometry.rod_100([3 4 1 2]); % [xlo xhi ylo yhi]
ROD_FULL = ROD_PFA;
ROD_GRE  = geometry.rod_128([3 4 1 2]);

% Two large boxes for panels A/B: 24x24 native voxels (= the original 31x31
% on the 128 grid), centred on the rod axis at the rod's two ends, where MIP
% shows the ring artefact that SOS suppresses.
rod_cx = mean(ROD_PFA(1:2));
BOX_A  = {[rod_cx-12, ROD_PFA(3)-12, 24, 24], ...  % around the top end   (ROI III side)
          [rod_cx-12, ROD_PFA(4)-12, 24, 24]};     % around the bottom end (ROI II side)

files = { ...
  {'MIP.png',    imgMIP_pfa,  ROD_PFA,  [],        BOX_A}, ...
  {'SOS.png',    imgSOS_pfa,  ROD_PFA,  [],        BOX_A}, ...
  {'GRE.png',    magmap,      ROD_GRE,  CPRF,      {}   }, ...
  {'BVSMSI.png', imgSOS_pfa,  ROD_PFA,  C2D,       {}   }, ...
  {'AVSMSI.png', imgSOS_full, ROD_FULL, C2D_full, {}   } };

for f = 1:numel(files)
    nm = files{f}{1}; im = files{f}{2}; rod = files{f}{3}; roi = files{f}{4}; boxes = files{f}{5};
    fig = figure('Color','w','Position',[80 80 420 420]);
    ax = axes('Position',[0 0 1 1]);
    imagesc(im); colormap(gray); axis image off; hold on;
    patch([rod(1) rod(1) rod(2) rod(2)], [rod(3) rod(4) rod(4) rod(3)], ...
        'r','FaceColor','red','EdgeColor','none');
    for b = 1:numel(boxes)
        rectangle('Position',boxes{b},'EdgeColor','y','LineWidth',1.2);
    end
    for k = 1:size(roi,1)
        rr = roi(k,1);  cc = roi(k,2);
        rectangle('Position',[cc-1.5, rr-1.5, 3, 3],'EdgeColor','y','LineWidth',1.0);
        text(cc, rr-4, roman(k), 'Color','y','FontSize',10,'FontWeight','bold','HorizontalAlignment','center');
    end
    hold off;
    % Lock the entire native FOV; annotation extents must not change cropping.
    xlim(ax,[0.5 size(im,2)+0.5]); ylim(ax,[0.5 size(im,1)+0.5]);
    set(ax,'Position',[0 0 1 1],'DataAspectRatio',[1 1 1]);
    set(fig,'Color','black','InvertHardcopy','off','PaperUnits','inches', ...
        'PaperPosition',[0 0 4.2 4.2],'PaperSize',[4.2 4.2]);
    print(fig,fullfile(OUT_DIR,nm),'-dpng','-r300');
    print(fig,fullfile(OUT_DIR,strrep(nm,'.png','.pdf')),'-dpdf','-painters');
    close(fig);
    fprintf('wrote %s\n', nm);
end

function s = roman(k), r = {'I','II','III','IV'}; s = r{k}; end
