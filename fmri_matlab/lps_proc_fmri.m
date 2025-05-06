%% load the image from .h5 file
fname_recon = "./lps_fmri_unfold_recon2.h5";
seq_args = fmriutl.loadh5struct(fname_recon, '/seq_args');
recon_args = fmriutl.loadh5struct(fname_recon, '/recon_args');
img_ts = abs(fmriutl.loadh5struct(fname_recon,'/sol').real + ...
    1i*fmriutl.loadh5struct(fname_recon,'/sol').imag);

%% remove discarded frames
ndiscard = 10;
nvol = size(img_ts,4);
img_ts = img_ts(:,:,:,ndiscard+1:end);

%% detrend the data
img_ts = fmriutl.poly_detrend(img_ts,1);

%% set general plotting parameters
slices2show = [15:62];

% create brain mask
msk = mean(abs(img_ts),4) > std(mean(abs(img_ts),4),[],1:3);

%% create the design matrix
tdur = 20; % task duration (single "on" period, s)
tdel = 0; % task onset delay (s)

tr = seq_args.tr*1e-3*recon_args.volwidth; % volume TR (s)
t = tr * (ndiscard:nvol-1); % time array (s)
ref = fmriutl.fmri_act(t,tdur,tdur,tdel);
A = ref(:).^[0,1];

figure
imagesc(A)

%% perform UNFOLD filtering on the object
vol_per_cyc = seq_args.nprj/recon_args.volwidth;
[img_ts_unfold,unfold_filt] = fmriutl.unfold3d(img_ts,vol_per_cyc);

%% calculate tscores
[tscore, beta] = fmriutl.fmri_tscore(A, abs(img_ts));
df_unfold = size(A,1) - size(A,2) - sum(unfold_filt == 0); % subtract out notch filtered points
[tscore_unfold, beta_unfold] = fmriutl.fmri_tscore(A, abs(img_ts_unfold), df_unfold);

%% plot images with and without UNFOLD
figure 

subplot(1,2,1)
im(img_ts(:,:,slices2show,1));
axis off
title('first frame pre-UNFOLD');

subplot(1,2,2)
im(img_ts_unfold(:,:,slices2show,1));
axis off
title('first frame post-UNFOLD');

%% plot tSNR with and without UNFOLD
tSNR = mean(img_ts,4) ./ std(img_ts,[],4);
tSNR_unfold = mean(img_ts_unfold,4) ./ std(img_ts_unfold,[],4);

figure

subplot(1,2,1)
im(tSNR(:,end:-1:1,slices2show(end:-2:1)));
axis off
colormap jet
cb = colorbar;
cb.FontSize = 12;
title('tSNR pre-UNFOLD');

subplot(1,2,2)
im(tSNR_unfold(:,end:-1:1,slices2show(end:-2:1)));
axis off
colormap jet
cb = colorbar;
cb.FontSize = 12;
title('tSNR post-UNFOLD');

%% view timeseries at ROI with and without UNFOLD
f = (-(nvol-ndiscard)/2:(nvol-ndiscard)/2-1) * 1/(tr * (nvol-ndiscard));
f_task = 1/(tdur*2); % task frequency (Hz)

% create ROI mask
off = [0.45,-0.1,0.3]; % ROI offset
sma = [0.15,0.15,0.15]; % ROI semi-major axes
roi = fmriutl.ellipsoid_mask(recon_args.M*ones(1,3),off,sma);

normabsmean = @(x,x_ref) (squeeze(mean(abs(x),1)) - min(mean(abs(x_ref),1))) ...
    ./ (max(mean(abs(x_ref),1)) - min(mean(abs(x_ref),1)));

subplot(4,2,[1,3]);
ov1 = fmriutl.overlayview;
ov1.addlayer(abs(img_ts(:,end:-1:1,slices2show,1))); % add underlay
ov1.addlayer(roi(:,end:-1:1,slices2show).*abs(img_ts(:,end:-1:1,slices2show,1))/max(abs(img_ts(:,:,slices2show,1)),[],'all'), ...
    'caxis',[eps,0.5],'cmap',[0;1].*[0,0,1]); % add ROI
ov1.show;
colorbar off
axis off
title('pre-UNFOLD');

subplot(4,2,5)
ts = reshape(img_ts,[],size(img_ts,4));
plot(t,ref,'--k'), hold on
plot(t,normabsmean(ts(msk,:),ts(msk,:)), 'color', 'r')
plot(t,normabsmean(ts(roi,:),ts(roi,:)), '-', 'color', [0,0,1], 'linewidth', 2), hold off
ylabel('signal')
xlabel('time (s)')
ylim([0,1])
title('timeseries')
legend('task','whole brain','ROI')

subplot(4,2,7)
spec = fftshift(fft(ts,[],2),2);
spec(:,floor(end/2)+1) = 0; % zero out DC term
plot(f,normabsmean(spec(msk,:),spec(msk,:)), 'color', 'r'), hold on
plot(f,normabsmean(spec(roi,:),spec(roi,:)), '-', 'color', [0,0,1], 'linewidth', 2), hold off
xline(f_task, '--k', 'label', 'task'); xline(-f_task, '--k');
ylabel('signal')
xlabel('freq (Hz)')
ylim([0,1])
title('spectrum (minus DC)')

subplot(4,2,[2,4]);
ov2 = fmriutl.overlayview;
ov2.addlayer(abs(img_ts_unfold(:,end:-1:1,slices2show,1))); % add underlay
ov2.addlayer(roi(:,end:-1:1,slices2show).*abs(img_ts_unfold(:,end:-1:1,slices2show,1))/max(abs(img_ts_unfold(:,:,slices2show,1)),[],'all'), ...
    'caxis',[eps,0.5],'cmap',[0;1].*[0,0,1]); % add ROI
ov2.show;
colorbar off
axis off
title('post-UNFOLD');

subplot(4,2,6)
ts_unfold = reshape(img_ts_unfold,[],size(img_ts,4));
plot(t,ref,'--k'), hold on
plot(t,normabsmean(ts_unfold(msk,:),ts(msk,:)), 'color', 'r')
plot(t,normabsmean(ts_unfold(roi,:),ts(roi,:)), '-', 'color', [0,0,1], 'linewidth', 2), hold off
ylabel('signal')
xlabel('time (s)')
ylim([0,1])
title('timeseries')

subplot(4,2,8)
spec_unfold = fftshift(fft(ts_unfold,[],2),2);
spec_unfold(:,floor(end/2)+1) = 0; % zero out DC term
plot(f,normabsmean(spec_unfold(msk,:),spec(msk,:)), 'color', 'r'), hold on
plot(f,normabsmean(spec_unfold(roi,:),spec(roi,:)), '-', 'color', [0,0,1], 'linewidth', 2), hold off
xline(f_task, '--k', 'label', 'task'); xline(-f_task, '--k');
ylabel('signal')
xlabel('freq (Hz)')
ylim([0,1])
title('spectrum (minus DC)')

%% plot the activation maps
slices2show = 20:23;

figure

subplot(2,1,1);
tscore_range = [0.5,2.5];
ov3 = fmriutl.overlayview;
ov3.addlayer(abs(img_ts(:,end:-1:1,slices2show,1))); % add underlay
% ov3.addlayer(tscore(:,:,slices2show,2), 'caxis', tscore_range, ...
%     'cmap', 'winter'); % add neg act
ov3.addlayer(tscore(:,end:-1:1,slices2show,2), 'caxis', tscore_range, ...
    'cmap', 'hot', 'name', 'tscore'); % add pos act
ov3.show('viewargs',{'row',1});
axis off
title('visual cortex pre-UNFOLD');

subplot(2,1,2);
tscore_range_unfold = [3,8];
ov4 = fmriutl.overlayview;
ov4.addlayer(abs(img_ts_unfold(:,end:-1:1,:,1))); % add underlay
% ov4.addlayer(-tscore_unfold(:,:,slices2show,2), 'caxis', tscore_range_unfold, ...
%     'cmap', 'winter'); % add neg act
ov4.addlayer(tscore_unfold(:,end:-1:1,:,2), 'caxis', tscore_range_unfold, ...
    'cmap', 'hot', 'name', 'tscore'); % add pos act
ov4.show('viewargs',{'mid3'},'shift',[0,0,10]);
axis off
title('visual cortex post-UNFOLD');