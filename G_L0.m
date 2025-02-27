clc
clear all
close all
warning off
%% data type
acquisition = 'simulation';     % simulation || experiments || in_vivo
phantom = 'contrast_speckle';   % resolution_distorsion || contrast_speckle || carotid_cross || carotid_long
suffix = 'simu';             % simu || expe
path_dataset = ['./PICMUS/database/',acquisition,'/',phantom,'/',phantom,'_',suffix,'_dataset_rf','.hdf5'];% update the folder directory
path_scan = ['./PICMUS/database/',acquisition,'/',phantom,'/',phantom,'_',suffix,'_scan','.hdf5'];% update the folder directory
addpath('./minFunc_2012')% update the folder directory
%% loading dataset
%-- load scan and dataset
scan = linear_scan();
scan.read_file(path_scan);
dataset = us_dataset();
dataset.read_file(path_dataset);
%-- define scan based on time axis
time = (150:1400).'/dataset.sampling_frequency+dataset.initial_time;
z_axis= time*dataset.c0/2;
scan = linear_scan(scan.x_axis,z_axis);
%% Image reconstruction using DAS
%-- receive apodization
%-- dynamically expanding receive aperture with hanning apodization
rx_f_number = 0.5;
rx_aperture = scan.z/rx_f_number;
rx_aperture_distance = abs(scan.x*ones(1,dataset.channels)-ones(scan.pixels,1)*dataset.probe_geometry(:,1).');
receive_apodization = tools.apodization(rx_aperture_distance,rx_aperture*ones(1,dataset.channels),'hanning');
%-- angular apodization
angular_apodization = ones(scan.pixels,dataset.firings);
%-- beamforming loop
beamformed_data = zeros(scan.pixels,1);
time_vector = dataset.initial_time+(0:(dataset.samples-1))/dataset.sampling_frequency;
w0 = 2*pi*dataset.modulation_frequency;
wb = waitbar(0,'DAS beamforming');
%
tgc1 = time_vector'./max(time_vector);
tgc1 = exp(4*tgc1);
%tgc1 = ones([1891,1]);


for pw=57
    %-- transmit delay
    transmit_delay = scan.z*cos(dataset.angles(pw))+scan.x*sin(dataset.angles(pw));
    for nrx=1:dataset.channels
        %-- progress bar
        step=(nrx + (pw-1)*dataset.channels)/length(dataset.angles)/dataset.channels;
        waitbar(step,wb,sprintf('DAS-RF beamforming %0.0f%%',step*100));
        %-- receive delay
        receive_delay = sqrt((dataset.probe_geometry(nrx,1)-scan.x).^2+(dataset.probe_geometry(nrx,3)-scan.z).^2);
        %-- total delay
        delay = (transmit_delay+receive_delay)/dataset.c0;
        %-- phase shift
        phase_shift = exp(1i.*w0*(delay-2*scan.z/dataset.c0));
        %-- beamformed data
        beamformed_data = beamformed_data + phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,tgc1.*dataset.data(:,nrx,pw),delay,'spline',0);
    end
end

close(wb);
beamformed_data(isnan(beamformed_data))=0;
DAS_RF = beamformed_data;
DAS_RF1 = DAS_RF./max(abs(DAS_RF));
DAS_ENV = abs(hilbert(reshape(DAS_RF1,[scan.Nz scan.Nx])));
%% Inverse beamforming step
load('Phi.mat');%loading the weighting matrix
b = dataset.data(:,:,pw); 
b = tgc1.*b; 
b = double(b(:));
eps = 1.5e-6;
mu = 0.1;       
lambda = 1.5;  
beta1 =0.5;    
beta2 =1;    
beta3 = 0.5;      
maxiter = 200;
%[x, objective] = G_L0_v0(Phi, b, eps, mu, lambda, beta1, beta2, beta3, maxiter);
tic
[x, objective,criterion] = G_L0_v0(Phi,DAS_RF, b, eps, mu, lambda, beta1, beta2, beta3, maxiter);%子函数
toc
L1_RF = x(:,end)./max(x(:,end));
L1_ENV = abs(hilbert(reshape(L1_RF,[scan.Nz scan.Nx])));
%% Show the corresponding beamformed images
dynamic_range = 60;
%-- setting axis limits (mm)
x_lim = [min(scan.x_axis) max(scan.x_axis)]*1e3;
z_lim = [min(scan.z_axis) max(scan.z_axis)]*1e3;
%-- compute dB values
DAS_Bmode = 20*log10(DAS_ENV./max(DAS_ENV(:))+0.001);
B_Bmode = 20*log10(L1_ENV./max(L1_ENV(:))+0.001);
vrange = [-dynamic_range 0];
%-- display image
figure
subplot(121)
imagesc((scan.x_axis)*1e3,(scan.z_axis)*1e3,DAS_Bmode);
shading flat; colormap gray; caxis(vrange); colorbar; hold on;
axis equal manual;
xlabel('x [mm]');
ylabel('z [mm]');
set(gca,'YDir','reverse');
set(gca,'fontsize',16);
axis([x_lim z_lim]);
title('DAS')
drawnow; hold off;
pause(0.5);
subplot(122)
imagesc((scan.x_axis)*1e3,(scan.z_axis)*1e3,B_Bmode);
shading flat; colormap gray; caxis(vrange); colorbar; hold on;
axis equal manual;
xlabel('x [mm]');
ylabel('z [mm]');
set(gca,'YDir','reverse');
set(gca,'fontsize',16);
axis([x_lim z_lim]);
title('L0-LOWRANK')
drawnow; hold off;