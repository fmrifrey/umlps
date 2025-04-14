%% load the image from .h5 file

nvol = size(img_lps,4);

%% set plotting parameters
slices2show = 20:2:80;
off = [0,-0.33,-0.12]; % ROI offset
sma = [0.5,0.2,0.15]; % ROI semi-major axes

% create the ROI mask
roi = fmriutl.ellipsoid_mask(seq_args.N*ones(1,3),off,sma);
msk = mean(abs(img_lps),4) > std(mean(abs(img_lps),4),[],1:3);

%% create the design matrix
tdur = 20; % task duration (s)
toff = 0; % task onset (s)

tr = seq_args.tr*1e-3*recon_args.volwidth; % volume TR (s)
t = tr * (0:nvol-1); % time array (s)
A = [ones(nvol,1), cos(pi*t(:)/tdur), sin(pi*t(:)/tdur)];

%% perform UNFOLD filtering on the object
vol_per_cyc = seq_args.nprj/recon_args.volwidth;
[img_lps_unfold,unfold_filt] = fmriutl.unfold3d(img_lps,vol_per_cyc);

%% calculate tscores
[tscore, beta] = fmriutl.fmri_tscore(A, abs(img_lps));
df_unfold = size(A,1) - size(A,2) - sum(unfold_filt == 0); % subtract out notch filtered points
[tscore_unfold, beta_unfold] = fmriutl.fmri_tscore(A, abs(img_lps_unfold), df_unfold);

%% plot images with and without UNFOLD
figure 

subplot(1,2,1)
im(img_lps(:,:,slices2show,1));
axis off
title('first frame pre-UNFOLD');

subplot(1,2,2)
im(img_lps_unfold(:,:,slices2show,1));
axis off
title('first frame post-UNFOLD');

%% plot tSNR with and without UNFOLD
tSNR = mean(img_lps,4) ./ std(img_lps,[],4);
tSNR_unfold = mean(img_lps_unfold,4) ./ std(img_lps_unfold,[],4);

figure

subplot(1,2,1)
im(tSNR(:,:,slices2show));
axis off
colorbar
title('tSNR pre-UNFOLD');

subplot(1,2,2)
im(tSNR_unfold(:,:,slices2show));
axis off
colorbar
title('tSNR post-UNFOLD');

%% view timeseries at ROI with and without UNFOLD
f = (-nvol/2:nvol/2-1) * 1/(tr * nvol);
f_task = 1/(tdur*2); % task frequency (Hz)

off = [0.1,-0.65,-0.21]; % ROI offset
sma = [0.18,0.12,0.1]; % ROI semi-major axes
roi = fmriutl.ellipsoid_mask(seq_args.N*ones(1,3),off,sma);

figure(1)

subplot(4,2,[1,3]);
ov1 = fmriutl.overlayview;
ov1.addlayer(abs(img_lps(:,:,slices2show,1))); % add underlay
ov1.addlayer(roi(:,:,slices2show).*abs(img_lps(:,:,slices2show,1))/max(abs(img_lps(:,:,slices2show,1)),[],'all'), ...
    'caxis',[eps,0.5],'cmap',[0;1].*[0,0,1]); % add ROI
ov1.show;
colorbar off
axis off
title('pre-UNFOLD');

subplot(4,2,5)
ts = reshape(img_lps,[],size(img_lps,4));
plot(t,squeeze(mean(abs(ts(msk,:)),1))), hold on
plot(t,squeeze(mean(abs(ts(roi,:)),1))), hold off
ylabel('mean magnitude signal')
xlabel('time (s)')
title('timeseries')
legend('whole brain','ROI')

subplot(4,2,7)
spec = fftshift(fft(ts,[],2),2);
spec(:,floor(end/2)+1) = 0; % zero out DC term
plot(f, squeeze(mean(abs(spec(msk,:)),1))), hold on
plot(f, squeeze(mean(abs(spec(roi,:)),1))), hold off
xline(f_task, '--k', 'label', 'task'); xline(-f_task, '--k');
ylabel('mean magnitude signal')
xlabel('freq (Hz)')
title('spectrum (minus DC)')

subplot(4,2,[2,4]);
ov2 = fmriutl.overlayview;
ov2.addlayer(abs(img_lps_unfold(:,:,slices2show,1))); % add underlay
ov2.addlayer(roi(:,:,slices2show).*abs(img_lps_unfold(:,:,slices2show,1))/max(abs(img_lps_unfold(:,:,slices2show,1)),[],'all'), ...
    'caxis',[eps,0.5],'cmap',[0;1].*[0,0,1]); % add ROI
ov2.show;
colorbar off
axis off
title('post-UNFOLD');

subplot(4,2,6)
ts_unfold = reshape(img_lps_unfold,[],size(img_lps,4));
plot(t,squeeze(mean(abs(ts_unfold(msk,:)),1))), hold on
plot(t,squeeze(mean(abs(ts_unfold(roi,:)),1))), hold off
ylabel('mean magnitude signal')
xlabel('time (s)')
title('timeseries')

subplot(4,2,8)
spec_unfold = fftshift(fft(ts_unfold,[],2),2);
spec_unfold(:,floor(end/2)+1) = 0; % zero out DC term
plot(f, squeeze(mean(abs(spec_unfold(msk,:)),1))), hold on
plot(f, squeeze(mean(abs(spec_unfold(roi,:)),1))), hold off
xline(f_task, '--k', 'label', 'task'); xline(-f_task, '--k');
ylabel('mean magnitude signal')
xlabel('freq (Hz)')
title('spectrum (minus DC)')

%% plot the activation maps
figure

subplot(1,2,1);
tscore_range = [1,5];
ov3 = fmriutl.overlayview;
ov3.addlayer(abs(img_lps(:,:,slices2show,1))); % add underlay
ov3.addlayer(-sum(tscore(:,:,slices2show,2:3),4), 'caxis', tscore_range, ...
    'cmap', 'winter'); % add neg act
ov3.addlayer(sum(tscore(:,:,slices2show,2:3),4), 'caxis', tscore_range, ...
    'cmap', 'hot', 'name', 'tscore'); % add pos act
ov3.show;
axis off
title('activation maps pre-UNFOLD');

subplot(1,2,2);
tscore_range_unfold = [2,7];
ov4 = fmriutl.overlayview;
ov4.addlayer(abs(img_lps_unfold(:,:,slices2show,1))); % add underlay
ov4.addlayer(-sum(tscore_unfold(:,:,slices2show,2:3),4), 'caxis', tscore_range_unfold, ...
    'cmap', 'winter'); % add neg act
ov4.addlayer(sum(tscore_unfold(:,:,slices2show,2:3),4), 'caxis', tscore_range_unfold, ...
    'cmap', 'hot', 'name', 'tscore'); % add pos act
ov4.show;
axis off
title('activation maps post-UNFOLD');