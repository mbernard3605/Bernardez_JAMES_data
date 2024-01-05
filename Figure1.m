addpath('mexcdf/mexnc/')
addpath('mexcdf/snctools/')

R=287.04; kap=2/7; p0=1000e2; cp=R/kap;  G=9.81;      sclht = R*256./G;
Lv = 2.5e6;
Rd = 287.06;
Cp = 1005; 
Lv = 2.5e6;
Lf = 3.33e5; 
lvlSize = 60;
g = 9.80665;
Tr = 273.15;
Pr = 1e5;
rhov=1000;
earthRad=6.371e6;

ERARad = '../ERA/radiation_new.nc';
ERAHFX = '../ERA/flux.nc';
ERASurf = '../ERA/sst.nc';
ERAVert = '../ERA/Omega.nc';
ERAprecip= '../ERA/rain.nc';
landFile = '../ERA/ERA5_land.nc';
land = nc_varget(landFile,'lsm');
land=cat(3,land(:,:,73:end),land(:,:,1:72));
lons = nc_varget(landFile,'longitude');
lats = nc_varget(landFile,'latitude');
lons=[lons(72:end); lons(1:71)+360];


SST = nc_varget(ERASurf,'sst');
time = nc_varget(ERARad,'time');
t=datestr(datenum('1900-01-01 00:00:0.0')+double(time)/24);
tInd = find(mod(time-time(1),12)==0);
lon = nc_varget(ERARad,'longitude');
lat = nc_varget(ERARad,'latitude');
west=[5 7.5 140 160];
east=[7.5 10 220 240];
lonSurf = nc_varget(ERAHFX,'longitude');
latSurf = nc_varget(ERAHFX,'latitude');
%sp=nc_varget(ERASurf,'SP_GDS4_SFC_S123');

[Y,X]=meshgrid(lat,lon);
x=BlockMean(X,2);
y=BlockMean(Y,2);
dy=(y(:,3:end)-y(:,1:end-2))/2;
dx=(x(3:end,:)-x(1:end-2,:))/2;
dxkm=deg2rad(mod(dx,360)).*earthRad.*cosd(y(2:end-1,:));
dykm=deg2rad(mod(-dy,360))*earthRad;


SST_mean=squeeze(nanmean(SST));
SST_small=interp2(Y,X,SST_mean',y,x)';

sst_dx=diff(SST_small,2,2)./dxkm'.^2;
sst_dy=diff(SST_small,2,1)./dykm'.^2;
sst_grad=sst_dx(2:end-1,:)+sst_dy(:,2:end-1);

N=13;
avg_filt=ones(N,N)./(N.^2);
sst_grad_smooth=filter2(avg_filt,sst_grad);

pressure = nc_varget(ERAVert,'level')*100;
dpressure = diff(pressure);
pressureshort = (pressure(1:end-1)+pressure(2:end))/2;

HFX = nc_varget(ERAHFX,'sshf')/86400;
LH = nc_varget(ERAHFX,'slhf')/86400;
Rain = nc_varget(ERAHFX,'tp')/86400;
Rain1 = Rain*1000;
SHflux = -squeeze(mean(HFX));
LHflux = -squeeze(mean(LH));
Rflux = squeeze(mean(Rain1));

precip_ls=nc_varget(ERAprecip,'mlspr');
precip_ls = squeeze(mean(precip_ls*86400));
precip_total=nc_varget(ERAprecip,'mtpr');
precip_total = squeeze(mean(precip_total*86400));
precip_conv=nc_varget(ERAprecip,'mcpr');
precip_conv = squeeze(mean(precip_conv*86400));
evap=nc_varget(ERAprecip,'mer');
evap = -squeeze(mean(evap*86400));

SSR = nc_varget(ERARad,'ssr')/86400;

STR = nc_varget(ERARad,'str')/86400;

TTR = nc_varget(ERARad,'ttr')/86400;

TSR = nc_varget(ERARad,'tsr')/86400;


radFlux = -SSR - STR + TSR + TTR;
radFlux = squeeze(mean(radFlux));


SST_mean=squeeze(nanmean(SST));
sst_land=isnan(SST_mean);
Omega=zeros([37,size(SST_mean)]);
for i = 1:length(time)/12
   Omega=Omega+squeeze(sum(nc_varget(ERAVert,'w',[(i-1)*12, 0, 0, 0],[12, inf, inf, inf])));

end

Omega = Omega/length(time);
%W=reshape(permute(Omega,[2 3 4 1]),37,[]);

OmegaShort = (Omega(1:end-1,:,:)+Omega(2:end,:,:))/2;
Omega_tmp=reshape(Omega,size(Omega,1),[]);

[pcs,lds]=Calc_TH_angle(Omega_tmp,pressure);
o1=reshape(-pcs(1,:),size(Omega,2),size(Omega,3));o2=reshape(pcs(2,:),size(Omega,2),size(Omega,3));

angles_map=atan2d(o2,o1);

%save('ERA5_angle_map.mat','angles_map','lon','lat');
% 
% figure,
% % plot(dSdp.*OmegaShort/g,pressureshort,'linewidth',2);
% set(gca,'ydir','reverse','ylim',[10000 100000]);
% xlabel('Omega ds/dp (J/s)');
% ylabel('Pressure (pa)');
% title('Omega ds/dp')
anglecolors=interp1(angles,testmap,angles_map,[],'extrap');

eastBox = [219.5 241.5 7 10.5]; %[east west sout north]
westBox = [139.5 161.5 4.5 8];

[Meshx,Meshy]=meshgrid(0:.25:360,-20:.25:20);
[eastBoxpointsx,eastBoxpointsy]=find((Meshx == eastBox(2) | Meshx == eastBox(1)) & (Meshy == eastBox(4) | Meshy == eastBox(3)));
eastPlotBox=zeros(size(Meshx));%-10;
eastPlotBox(eastBoxpointsx(1),eastBoxpointsy(1):eastBoxpointsy(3))=1;
eastPlotBox(eastBoxpointsx(2),eastBoxpointsy(1):eastBoxpointsy(3))=1;
eastPlotBox(eastBoxpointsx(1):eastBoxpointsx(2),eastBoxpointsy(1))=1;
eastPlotBox(eastBoxpointsx(1):eastBoxpointsx(2),eastBoxpointsy(3))=1;

[westBoxpointsx,westBoxpointsy]=find((Meshx == westBox(2) | Meshx == westBox(1)) & (Meshy == westBox(4) | Meshy == westBox(3)));
westPlotBox=zeros(size(Meshx));%-10;
westPlotBox(westBoxpointsx(1),westBoxpointsy(1):westBoxpointsy(3))=1;
westPlotBox(westBoxpointsx(2),westBoxpointsy(1):westBoxpointsy(3))=1;
westPlotBox(westBoxpointsx(1):westBoxpointsx(2),westBoxpointsy(1))=1;
westPlotBox(westBoxpointsx(1):westBoxpointsx(2),westBoxpointsy(3))=1;




figure,
ax(2)=subplot(3,1,1);
[c,h]=contourf([lon(length(lon)/2+1:end); lon(1:length(lon)/2)+360],lat,[SST_mean(:,length(lon)/2+1:end) SST_mean(:,1:length(lon)/2)],[296:303]);
axis equal;
set(h,'linecolor','none');set(gca,'fontsize',20);
colormap(ax(2),redblue);colorbar;
caxis([296 303]);
hold on;
contour([lon(length(lon)/2:end-1); lon(1:length(lon)/2)+360],lat,[sst_land(:,length(lon)/2+1:end) sst_land(:,1:length(lon)/2)],[1 1],'k','linewidth',2);
set(gca,'ylim',[-20 20],'xlim',[0 360]);
title('SST');
contour(Meshx,Meshy,eastPlotBox,[1 1],'k','linewidth',2);
contour(Meshx,Meshy,westPlotBox,[1 1],'k','linewidth',2);
ax(3)=subplot(3,1,3);
[c,h]=contourf([lon(length(lon)/2:end-1); lon(1:length(lon)/2)+360],lat,[LHflux(:,length(lon)/2+1:end) LHflux(:,1:length(lon)/2)],[80:10:160]);
axis equal;
set(h,'linecolor','none');set(gca,'fontsize',20);
colormap(ax(3),redblue);colorbar;
caxis([50 200]);
hold on;
contour([lon(length(lon)/2:end-1); lon(1:length(lon)/2)+360],lat,[sst_land(:,length(lon)/2+1:end) sst_land(:,1:length(lon)/2)],[1 1],'k','linewidth',2);
set(gca,'ylim',[-20 20],'xlim',[0 360]);
title('Surface Latent Heat Flux');
contour(Meshx,Meshy,eastPlotBox,[1 1],'k','linewidth',2);
contour(Meshx,Meshy,westPlotBox,[1 1],'k','linewidth',2);
ax(4)=subplot(3,1,2);
[c,h]=contourf([x(length(x(2:end-1,1))/2+1:end-1,1); x(2:length(x(2:end-1,1))/2,1)+360],y(1,2:end-1),[sst_grad_smooth(:,length(x(2:end-1,1))/2+1:end) sst_grad_smooth(:,1:length(x(2:end-1,1))/2)],linspace(-1e-11,0,8));
axis equal;
set(h,'linecolor','none');set(gca,'fontsize',20);
colormap(ax(4),flip(redblue));cbh=colorbar;
caxis([-10e-12 0]);
hold on;
contour([lon(length(lon)/2:end-1); lon(1:length(lon)/2)+360],lat,[sst_land(:,length(lon)/2+1:end) sst_land(:,1:length(lon)/2)],[1 1],'k','linewidth',2);
set(gca,'ylim',[-20 20],'xlim',[0 360]);
title('Laplacian of SST');
contour(Meshx,Meshy,eastPlotBox,[1 1],'k','linewidth',2);
contour(Meshx,Meshy,westPlotBox,[1 1],'k','linewidth',2);
%save('ERA5_sst_land.mat','sst_land');


