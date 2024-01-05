directory = '../WRF_OUT/average/ew_experiment/basic/';

east1='ERA5_force_east.nc';
west1='ERA5_force_west.nc';
dyn1='WRF_OUT/Dynamo_force_conwind.nc';



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

east1=[directory 'wrfout_east_DE_cor_qvhadv_intrad_avg'];
east2=[directory 'wrfout_east_DE_cor_qvhadv_imposedrad_avg'];

west1=[directory 'wrfout_west_DE_cor_qvhadv_intrad_avg'];
west2=[directory 'wrfout_west_DE_cor_qvhadv_imposedrad_avg'];

east3=[directory 'wrfout_east_DE_TRMM_qvhadv_intrad_avg'];
east4=[directory 'wrfout_east_DE_TRMM_qvhadv_imposedrad_avg'];

west3=[directory 'wrfout_west_DE_TRMM_qvhadv_intrad_avg'];
west4=[directory 'wrfout_west_DE_TRMM_qvhadv_imposedrad_avg'];

east5=[directory 'wrfout_east_ERA5_DE_qvhadv_intrad_wtgopt10_avg'];
east6=[directory 'wrfout_east_ERA5_DE_qvhadv_imprad_wtgopt10_avg'];

west5=[directory 'wrfout_west_ERA5_DE_qvhadv_intrad_wtgopt10_avg'];
west6=[directory 'wrfout_west_ERA5_DE_qvhadv_imprad_wtgopt10_avg'];

dyn1=[directory 'wrfout_Dynamo_DE_fullwind_qvhadv_intrad_avg'];
dyn2=[directory 'wrfout_Dynamo_DE_qvhadv_imposedrad_avg'];
dyn3=[directory 'wrfout_Dynamo_DE_qvhadv_intrad_avg'];


pe = nc_varget(east1,'P')+nc_varget(east1,'PB');
pplot=squeeze(mean(pe(end-120:end,:)))/100;

omegaEast4 = nc_varget(east4,'OMEGA_WTG');
omegaEast5 = nc_varget(east5,'OMEGA_WTG');
omegaEast6 = nc_varget(east6,'OMEGA_WTG');

omegaWest4 = nc_varget(west4,'OMEGA_WTG');
omegaWest5 = nc_varget(west5,'OMEGA_WTG');
omegaWest6 = nc_varget(west6,'OMEGA_WTG');

omegaDyn1 = nc_varget(dyn1,'OMEGA_WTG');
omegaDyn2 = nc_varget(dyn2,'OMEGA_WTG');
omegaDyn3 = nc_varget(dyn3,'OMEGA_WTG');

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
print(gcf,'-djpeg','../PaperPlots/ERA5_swtg_basic_w.jpg');
print(gcf,'-dpng','../PaperPlots/ERA5_swtg_basic_w.png');