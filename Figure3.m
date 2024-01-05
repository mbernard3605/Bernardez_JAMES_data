load('anglecolormap5.mat');


east='../ERA/ERA5/ERA5_force_east.nc';
west='../ERA/ERA5/ERA5_force_west.nc';

ple2=nc_varget(east,'P_FORCE');
wle2=nc_varget(east,'W_SUBS');
Te2=nc_varget(east,'TH_LARGESCALE');
QVHe2=nc_varget(east,'QV_LARGESCALE_TEND')*2.5e6;
QVe2=nc_varget(east,'QV_LARGESCALE');
the2=Te2.*(ple2./1000).^(2/7);
rhoe2=ple2*100./the2/287;
omegae2=-wle2.*9.81.*rhoe2;

plw2=nc_varget(west,'P_FORCE');
wlw2=nc_varget(west,'W_SUBS');
Tw2=nc_varget(west,'TH_LARGESCALE');
QVHw2=nc_varget(west,'QV_LARGESCALE_TEND')*2.5e6;
QVw2=nc_varget(west,'QV_LARGESCALE');
thw2=Tw2.*(plw2./1000).^(2/7);
rhow2=plw2*100./thw2/287;
omegaw2=-wlw2.*9.81.*rhow2;

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


figure('units','normalized','outerposition',[1.01 0.01 .5 .6])

subplot('position',[.1 .45 .3 .5]);
ax1=plot(nanmean(omegae2),nanmean(ple2),'b',nanmean(omegaw2),nanmean(plw2),'r','linewidth',2);
legend(ax1,{'Bottom-heavy','Top-Heavy'},'location','northeast','fontsize',8);
grid on;
title('Vertical motion profiles');
xlabel('Pressure Velocity [Pa/s]','fontsize',8)
set(gca,'ylim',[100 1000],'ydir','reverse','fontsize',12,'xdir','reverse');

load('ERA5_angle_map.mat');
load('ERA5_EOFS.mat');
load('ERA5_sst_land.mat','sst_land');
weights_old=weights;

eof=lds./repmat(weights_old.^0.5,1,size(lds,2)); 

angleIndex=linspace(-180,180,60);
r=1;
angles = -180:45:180 ;
o1_plot = r*cosd(angles);
o2_plot = r*sind(angles);
omega_legend_plot=-eof(:,1)*o1_plot+eof(:,2)*o2_plot;
colors_plot=interp1(angleIndex',testmap,angles,[],'extrap');
only_legend=1;


subplot('position',[.52 .6 .13 .18]);
ax(1)=plot(omega_legend_plot(:,1),level,'color',colors_plot(1,:),'linewidth',2);
set(gca,'ylim',[100 1000]*100,'ydir','reverse','xlim',[-2 2],'yticklabels',[],'xticklabels',[],'xdir','reverse','YAxisLocation','right');
grid on;
ylabel('\bf Descent','fontsize',10,'rotation',-90,'VerticalAlignment','bottom');
subplot('position',[.52 .8 .13 .18]);
ax(2)=plot(omega_legend_plot(:,end-1),level,'color',colors_plot(end-1,:),'linewidth',2);
set(gca,'ylim',[100 1000]*100,'ydir','reverse','xlim',[-2 2],'xticklabels',[],'xdir','reverse');
grid on;
subplot('position',[.67 .8 .13 .18]);
ax(3)=plot(omega_legend_plot(:,end-2),level,'color',colors_plot(end-2,:),'linewidth',2);
set(gca,'ylim',[100 1000]*100,'ydir','reverse','xlim',[-2 2],'yticklabels',[],'xticklabels',[],'xdir','reverse');
grid on;
xlabel('\bf Top-Heavy','fontsize',10);
subplot('position',[.82 .8 .13 .18]);
ax(4)=plot(omega_legend_plot(:,end-3),level,'color',colors_plot(end-3,:),'linewidth',2);
set(gca,'ylim',[100 1000]*100,'ydir','reverse','xlim',[-2 2],'yticklabels',[],'xticklabels',[],'xdir','reverse');
grid on;
subplot('position',[.82 .6 .13 .18]);
ax(5)=plot(omega_legend_plot(:,5),level,'color',colors_plot(5,:),'linewidth',2);
set(gca,'ylim',[100 1000]*100,'ydir','reverse','xlim',[-2 2],'yticklabels',[],'xticklabels',[],'xdir','reverse');
grid on;
ylabel('\bf Ascent','fontsize',10);
subplot('position',[.82 .4 .13 .18]);
ax(6)=plot(omega_legend_plot(:,4),level,'color',colors_plot(4,:),'linewidth',2);
set(gca,'ylim',[100 1000]*100,'ydir','reverse','xlim',[-2 2],'yticklabels',[],'xticklabels',[],'xdir','reverse');
grid on;
subplot('position',[.67 .4 .13 .18]);
ax(7)=plot(omega_legend_plot(:,3),level,'color',colors_plot(3,:),'linewidth',2);
set(gca,'ylim',[100 1000]*100,'ydir','reverse','xlim',[-2 2],'yticklabels',[],'xticklabels',[],'xdir','reverse');
grid on;
title('Bottom-Heavy','fontsize',10);
subplot('position',[.52 .4 .13 .18]);
ax(8)=plot(omega_legend_plot(:,2),level,'color',colors_plot(2,:),'linewidth',2);
set(gca,'ylim',[100 1000]*100,'ydir','reverse','xlim',[-2 2],'yticklabels',[],'xticklabels',[],'xdir','reverse');
grid on;
h1=subplot('position',[.67 .6 .11 .16]);
colormap(testmap);
phasebar('deg','size',.9);
set(h1,'visible','off')

angles_map(sst_land)=NaN;
subplot('position',[.05 .06 .9 .28]);
[h]=pcolor([lon(length(lon)/2+1:end); lon(1:length(lon)/2)+360],lat,[angles_map(:,length(lon)/2+1:end) angles_map(:,1:length(lon)/2)]);
%axis equal;
set(h,'edgecolor','none')
colormap(testmap);
hold on;
contour([lon(length(lon)/2:end-1); lon(1:length(lon)/2)+360],lat,[sst_land(:,length(lon)/2+1:end) sst_land(:,1:length(lon)/2)],[1 1],'k','linewidth',2);
set(gca,'ylim',[-20 20],'xlim',[0 360]);
%title('Top-Heaviness Angle','fontsize',24);
hold on;
contour(Meshx,Meshy,eastPlotBox,[1 1],'k','linewidth',1);
contour(Meshx,Meshy,westPlotBox,[1 1],'k','linewidth',1);

print(gcf,'-dpng','../PaperPlots/map_and_legend.png');