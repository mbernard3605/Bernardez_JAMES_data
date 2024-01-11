directory = '../WRF_OUT/average/ew_experiment/basic/';

east1='ERA5_force_east.nc';
west1='ERA5_force_west.nc';
dyn1='Dynamo_force_conwind.nc';



ple=nc_varget(east1,'P_FORCE');
wle=nc_varget(east1,'W_SUBS');
Te=nc_varget(east1,'TH_LARGESCALE');
the=Te.*(ple./1000).^(2/7);
rhoe=ple*100./the/287;
omegae=-wle.*9.81.*rhoe;

plw=nc_varget(west1,'P_FORCE');
wlw=nc_varget(west1,'W_SUBS');
Tw=nc_varget(west1,'TH_LARGESCALE');
thw=Tw.*(plw./1000).^(2/7);
rhow=plw*100./thw/287;
omegaw=-wlw.*9.81.*rhow;

pld=nc_varget(dyn1,'P_FORCE');
wld=nc_varget(dyn1,'W_SUBS');
Td=nc_varget(dyn1,'TH_LARGESCALE');
thd=Td.*(pld./1000).^(2/7);
rhod=pld*100./thd/287;
omegad=-wld.*9.81.*rhod;

directory = '';


east6=[directory 'wrfout_east_ERA5_DE_qvhadv_imprad_wtgopt10_avg'];

west6=[directory 'wrfout_west_ERA5_DE_qvhadv_imprad_wtgopt10_avg'];

pe = nc_varget(east6,'P')+nc_varget(east6,'PB');
pplot=squeeze(mean(pe(end-120:end,:)))/100;

omegaEast6 = nc_varget(east6,'OMEGA_WTG');

omegaWest6 = nc_varget(west6,'OMEGA_WTG');


figure('units','normalized','outerposition',[1.25 .1 .6 .7]),
ax1=subplot(1,2,1);
plot(omegaw(2,:),plw(2,:),'r',mean(omegaWest6(end-120:end,:)),pplot,'r--','linewidth',3.5);
title('Top-Heavy');
xlabel('Vertical Velocity [Pa/s]');
ylabel('Pressure [hPa]');
set(gca,'ylim',[100 1000],'xlim',[-.1 .05],'ydir','reverse');grid on;
lg1=legend('Reanalysis','SWTG','location',[.286 .666 .214 .103]);
ax2=subplot(1,2,2);
plot(omegae(2,:),ple(2,:),'b',mean(omegaEast6(end-120:end,:)),pplot,'b--','linewidth',3.5);
title('Bottom-Heavy');
xlabel('Vertical Velocity [Pa/s]');
set(gca,'ylim',[100 1000],'xlim',[-.1 .05],'ydir','reverse');
set([ax1 ax2],'fontsize',16);grid on;
% print(gcf,'-djpeg','../PaperPlots/ERA5_swtg_basic_w.jpg');
% print(gcf,'-dpng','../PaperPlots/ERA5_swtg_basic_w.png');