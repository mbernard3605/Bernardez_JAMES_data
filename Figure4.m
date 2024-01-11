
east1='wrfout_east_ERA5_DE_corr_qvhadv_avg';
east2='ERA5_force_east.nc';
west1='wrfout_west_ERA5_DE_corr_qvhadv_avg';
west2='ERA5_force_west.nc';



ple=nc_varget(east1,'P')+nc_varget(east1,'PB');
wle=nc_varget(east1,'OMEGA_WTG');
Te=nc_varget(east1,'T')+300;
QVe=nc_varget(east1,'QVAPOR');
the=Te.*(ple./100000).^(2/7);
QVSe=calc_qvstar(the,ple);

plw=nc_varget(west1,'P')+nc_varget(west1,'PB');
wlw=nc_varget(west1,'OMEGA_WTG');
Tw=nc_varget(west1,'T')+300;
QVw=nc_varget(west1,'QVAPOR');
thw=Tw.*(plw./100000).^(2/7);
QVSw=calc_qvstar(thw,plw);


ple2=nc_varget(east2,'P_FORCE');
wle2=nc_varget(east2,'W_SUBS');
Te2=nc_varget(east2,'TH_LARGESCALE');
QVHe2=nc_varget(east2,'QV_LARGESCALE_TEND')*2.5e6;
QVe2=nc_varget(east2,'QV_LARGESCALE');
the2=Te2.*(ple2./1000).^(2/7);
QVSe2=calc_qvstar(the2,ple2*100);
rhoe2=ple2*100./the2/287;
omegae2=-wle2.*9.81.*rhoe2;

plw2=nc_varget(west2,'P_FORCE');
wlw2=nc_varget(west2,'W_SUBS');
Tw2=nc_varget(west2,'TH_LARGESCALE');
QVHw2=nc_varget(west2,'QV_LARGESCALE_TEND')*2.5e6;
QVw2=nc_varget(west2,'QV_LARGESCALE');
thw2=Tw2.*(plw2./1000).^(2/7);
QVSw2=calc_qvstar(thw2,plw2*100);
rhow2=plw2*100./thw2/287;
omegaw2=-wlw2.*9.81.*rhow2;


figure(1),
set(gcf,'units','normalized','position',[1.1 .2 .8 .6]);
ax1=subplot(1,4,3);
plot(omegae2(2,:),ple2(2,:),'b',omegaw2(2,:),plw2(2,:),'r','linewidth',2);
set(gca,'ylim',[100 1000],'ydir','reverse','xlim',[-.1 0]);
title('Vertical Velocity');
grid on;
ax2=subplot(1,4,4);
plot(QVHe2(2,:),ple2(2,:),'b',QVHw2(2,:),plw2(2,:),'r','linewidth',2);
set(gca,'ylim',[100 1000],'ydir','reverse','xlim',[-0.05 0.01]);
legend('Bottom-heavy','Top-heavy','location','northwest');
title('Horizontal Advection');
grid on;
ax3=subplot(1,4,1);
plot(Te2(2,:)-Tw2(2,:),plw2(2,:),'k',nanmean(Te(end-119:end,:))-nanmean(Tw(end-119:end,:)),plw(2,:)/100,'k--','linewidth',2);
set(gca,'ylim',[100 1000],'ydir','reverse','xlim',[-2.5 0.5]);
title('T difference');
legend('Reanalysis','DE');
grid on;
ax4=subplot(1,4,2);
plot(QVe2(2,:)./QVSe2(2,:)-QVw2(2,:)./QVSw2(2,:),plw2(2,:),'k',nanmean(QVe(end-110:end,:)./QVSe(end-110:end,:))-nanmean(QVw(end-110:end,:)./QVSw(end-110:end,:)),plw(2,:)/100,'k--',nanmean(QVe(end-119:end,:))./nanmean(QVSe(end-119:end,:)),ple(2,:),'b--',nanmean(QVw(end-119:end,:))./nanmean(QVSw(end-119:end,:)),plw(2,:),'r--','linewidth',2);
set(gca,'ylim',[100 1000],'ydir','reverse');
title('RH difference');
grid on;
set([ax1 ax2 ax3 ax4],'fontsize',16);
