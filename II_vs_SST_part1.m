%addpath('mexcdf/mexnc/')
%addpath('mexcdf/snctools/')
clear variables;
close all;
clc;

load('anglecolormap4.mat');
angles=linspace(-180,180,62);

R=287.04; kap=2/7; p0=1000e2; cp=R/kap;  G=9.81;      sclht = R*256./G;
Rd = 287.06;
earthRad=6.371e6;

Cp = 1005; 
Lv = 2.5e6;
Lf = 3.33e5;
lvlSize = 60;
g = 9.80665;
Tr = 273.15;
Pr = 1e5;
rhov=1000;

ERARad = '../ERA/radiation.nc';
ERAHFX = '../ERA/flux.nc';
ERASurf = '../ERA/surface.nc';
ERAVert = '../ERA/Pressure.nc';

SST = nc_varget(ERASurf,'SSTK_GDS4_SFC_S123');
time = nc_varget(ERARad,'initial_time0');
t=datestr(datenum('1900-01-01 00:00:0.0')+double(time)/24);
tInd = find(mod(time-time(1),12)==0);
lon = nc_varget(ERARad,'g4_lon_2');
lat = nc_varget(ERARad,'g4_lat_1');
latin = find(lat <= 10 & lat >= 7.5);
lonin = find(lon >= 220 & lon <= 240);
lonSurf = nc_varget(ERAHFX,'longitude');
latSurf = nc_varget(ERAHFX,'latitude');
latinVert = latin;
latinSurf=find(latSurf <= 10 & latSurf >= 7.5);
latinHorz = [latinVert(1)-1; latinVert; latinVert(end)+1];
loninVert = lonin;
loninSurf=find(lonSurf <= -120 & lonSurf >= -140);
loninHorz = [loninVert(1)-1; loninVert; loninVert(end)+1];
sp=nc_varget(ERASurf,'SP_GDS4_SFC_S123');
a=sp(:,latin,lonin);

[Y,X]=meshgrid(lat,lon);
dy=Y(:,3:end)-Y(:,1:end-2);
dx=X(3:end,:)-X(1:end-2,:);
dxkm=deg2rad(mod(dx,360)).*earthRad.*cosd(Y(2:end-1,:));
dykm=deg2rad(mod(-dy,360))*earthRad;


sst_dx=diff(SST,2,3)./permute(repmat(dxkm'.^2,1,1,132),[3 1 2]);
sst_dy=diff(SST,2,2)./permute(repmat(dykm'.^2,1,1,132),[3 1 2]);
sst_grad=squeeze(nanmean(nanmean(sst_dx(:,latin,lonin+1)+sst_dy(:,latin+1,lonin),2),3));
       

pressure = nc_varget(ERAVert,'lv_ISBL1')*100;
dpressure = diff(pressure);
pressureshort = (pressure(1:end-1)+pressure(2:end))/2;

HFX = nc_varget(ERAHFX,'sshf');
HFX = HFX(:,latinSurf,loninSurf)/86400;
LH = nc_varget(ERAHFX,'slhf');
LH = LH(:,latinSurf,loninSurf)/86400;
Rain = nc_varget(ERAHFX,'tp');
Rain1 = Rain(:,latinSurf,loninSurf)*1000;
SHflux = -mean(mean(HFX,2),3);
LHflux = -mean(mean(LH,2),3);
Rflux = mean(mean(Rain1,2),3);

SSR = nc_varget(ERARad,'SSR_GDS4_SFC_120');
SSR = SSR(:,latin,lonin)/86400;

STR = nc_varget(ERARad,'STR_GDS4_SFC_120');
STR = STR(:,latin,lonin)/86400;

TTR = nc_varget(ERARad,'TTR_GDS4_SFC_120');
TTR = TTR(:,latin,lonin)/86400;

TSR = nc_varget(ERARad,'TSR_GDS4_SFC_120');
TSR = TSR(:,latin,lonin)/86400;

radFlux = -SSR - STR + TSR + TTR;
radFlux = mean(mean(radFlux,2),3);

D=radFlux+SHflux(1:length(radFlux));

T = nc_varget(ERAVert,'T_GDS4_ISBL_S123');
SST_box=squeeze(nanmean(nanmean(SST(:,latin,lonin),2),3));

Q = nc_varget(ERAVert,'Q_GDS4_ISBL_S123');
q=squeeze(mean(mean(Q(:,:,latinVert,loninVert),3),4));
tv = T.*(1+.608*Q);

rho = pressure'./(tv*Rd);

geopotential = nc_varget(ERAVert,'Z_GDS4_ISBL_S123');
gp=squeeze(mean(mean(geopotential(:,:,latinVert,loninVert),3),4));
S = cp*T + geopotential;
H = S + Q*Lv;
Sp=squeeze(mean(mean(S(:,:,latinVert,loninVert),3),4));
Hp=squeeze(mean(mean(H(:,:,latinVert,loninVert),3),4));
Qp=squeeze(mean(mean(Lv*Q(:,:,latinVert,loninVert),3),4));

Omega = nc_varget(ERAVert,'W_GDS4_ISBL_S123');
Omega = Omega(:,:,latinVert,loninVert);
W=reshape(permute(Omega,[2 3 4 1]),37,[]);

Omega = squeeze(mean(mean(Omega,3),4));
OmegaShort = (Omega(:,1:end-1)+Omega(:,2:end))/2;



[pcs,lds]=Calc_TH_angle(Omega,pressure);
[pcs_norm,lds_norm]=Calc_TH_angle_norm(Omega,pressure);
o1=-pcs(1,:);o2=pcs(2,:);
lds(:,1)=-lds(:,1);
o1_norm=-pcs_norm(1,:);o2_norm=pcs_norm(2,:);
lds_norm(:,1)=-lds_norm(:,1);
angle=atan2d(o2,o1);
radius_norm = sqrt(o1_norm.^2+o2_norm.^2);
angle_norm=atan2d(o2_norm,o1_norm);
radius = sqrt(o1_norm.^2+o2_norm.^2);
radius_mean=sqrt(mean(o1)^2+mean(o2)^2);
angle_mean=atan2d(mean(o2),mean(o1));


U = nc_varget(ERAVert,'U_GDS4_ISBL_S123');
U=squeeze(mean(mean(U(:,:,latinVert,loninVert),3),4));
Usfc = nc_varget(ERASurf,'10U_GDS4_SFC_S123');
Usfc=mean2(Usfc(:,latinVert,loninVert));
V = nc_varget(ERAVert,'V_GDS4_ISBL_S123');
V=squeeze(mean(mean(V(:,:,latinVert,loninVert),3),4));
Vsfc = nc_varget(ERASurf,'10V_GDS4_SFC_S123');
Vsfc=mean2(Vsfc(:,latinVert,loninVert));
spd=sqrt(Usfc.^2+Vsfc.^2);
disp(spd);

dSdp = (Sp(:,1:end-1)-Sp(:,2:end))./(dpressure');
dHdp = (Hp(:,1:end-1)-Hp(:,2:end))./(dpressure');
dQdp = (Qp(:,1:end-1)-Qp(:,2:end))./(dpressure');

dSdv = squeeze(mean(mean(S(:,:,latinHorz(3:end),loninVert)-S(:,:,latinHorz(1:end-2),loninVert),3),4))/(2*78.875*1000);
dSdu = squeeze(mean(mean(S(:,:,latinVert,loninHorz(3:end))-S(:,:,latinVert,loninHorz(1:end-2)),3),4))/(2*78.875*1000);
dHdv = squeeze(mean(mean(H(:,:,latinHorz(3:end),loninVert)-H(:,:,latinHorz(1:end-2),loninVert),3),4))/(2*78.875*1000);
dHdu = squeeze(mean(mean(H(:,:,latinVert,loninHorz(3:end))-H(:,:,latinVert,loninHorz(1:end-2)),3),4))/(2*78.875*1000);

MSE=zeros(size(dSdp,1),size(lds,2));
MQE=zeros(size(dSdp,1),size(lds,2));
MHE=zeros(size(dSdp,1),size(lds,2));
MSE_norm=zeros(size(dSdp,1),size(lds,2));
MQE_norm=zeros(size(dSdp,1),size(lds,2));
MHE_norm=zeros(size(dSdp,1),size(lds,2));

for i = 1:size(lds,2)
    
   MSE(:,i)=-trapz(pressureshort',(dSdp.*repmat(midData(lds(:,i))',size(dSdp,1),1))'/9.81);
   MQE(:,i)=-trapz(pressureshort',(dQdp.*repmat(midData(lds(:,i))',size(dHdp,1),1))'/9.81);
   MHE(:,i)=-trapz(pressureshort',(dHdp.*repmat(midData(lds(:,i))',size(dHdp,1),1))'/9.81);
   MSE_norm(:,i)=-trapz(pressureshort',(dSdp.*repmat(midData(lds_norm(:,i))',size(dSdp,1),1))'/9.81);
   MQE_norm(:,i)=-trapz(pressureshort',(dQdp.*repmat(midData(lds_norm(:,i))',size(dHdp,1),1))'/9.81);    
   MHE_norm(:,i)=-trapz(pressureshort',(dHdp.*repmat(midData(lds_norm(:,i))',size(dHdp,1),1))'/9.81);    
    
end


MSE1=trapz(pressureshort',-(dSdp.*repmat(midData(lds(:,1))',size(dSdp,1),1))'/9.81);
MSE2=trapz(pressureshort',-(dSdp.*repmat(midData(lds(:,2))',size(dSdp,1),1))'/9.81);
MSE1_norm=trapz(pressureshort',-(dSdp.*repmat(midData(lds_norm(:,1))',size(dSdp,1),1))'/9.81);
MSE2_norm=trapz(pressureshort',-(dSdp.*repmat(midData(lds_norm(:,2))',size(dSdp,1),1))'/9.81);
MHE1=trapz(pressureshort',-(dHdp.*repmat(midData(lds(:,1))',size(dSdp,1),1))'/9.81);
MHE2=trapz(pressureshort',-(dHdp.*repmat(midData(lds(:,2))',size(dSdp,1),1))'/9.81);
MHE1_norm=trapz(pressureshort',-(dHdp.*repmat(midData(lds_norm(:,1))',size(dSdp,1),1))'/9.81);
MHE2_norm=trapz(pressureshort',-(dHdp.*repmat(midData(lds_norm(:,2))',size(dSdp,1),1))'/9.81);
disp(mean(MSE1)/mean(MSE2))
disp(mean(MSE1_norm)/mean(MSE2_norm))
disp(mean(MHE1)/mean(MHE2))
disp(mean(MHE1_norm)/mean(MHE2_norm))

% figure,
% % plot(dSdp.*OmegaShort/g,pressureshort,'linewidth',2);
% set(gca,'ydir','reverse','ylim',[10000 100000]);
% xlabel('Omega ds/dp (J/s)');
% ylabel('Pressure (pa)');
% title('Omega ds/dp')
for i = 1:size(dSdp,1)
int_S(i) = trapz(pressureshort(11:end),-dSdp(i,11:end).*OmegaShort(i,11:end)/g);
int_Sh(i) = trapz(pressure(11:end),-(dSdu(i,11:end).*U(i,11:end)+dSdv(i,11:end).*V(i,11:end))/g);
int_H(i) = trapz(pressureshort(11:end),-dHdp(i,11:end).*OmegaShort(i,11:end)/g);
int_Hh(i) = trapz(pressure(11:end),-(dHdu(i,11:end).*U(i,11:end)+dHdv(i,11:end).*V(i,11:end))/g);
int_Q(i) = trapz(pressureshort(11:end),-dQdp(i,11:end).*OmegaShort(i,11:end)/g);
end
GMS = int_H./int_S;
GMS_mean=mean(int_H)./mean(int_S);
GMS_tot=(int_H+int_Hh)./(int_S+int_Sh);

gamma=int_S(1:length(D))'\D;
disp(gamma);

es =  610.78*exp(17.269388*(T-273.16)./(T-35.86));
es2=611.2*exp((23.036-(T-273.15)/333.7).*(T./(279.82-273.15+T)));
QVS = .622*es./(pressure'-es);
QVS2 = .622*es2./(pressure'-es2);
QVS(QVS<0)=0;
qvs=squeeze(mean(mean(QVS(:,:,latinVert,loninVert),3),4));
ents = Cp*log(T/273.15)-Rd*log(pressure'/100000)+Lv*QVS/273.15;
entS=squeeze(mean(mean(ents(:,:,latinVert,loninVert),3),4));
II = mean(entS(:,(gp(1,:)/9.81<=3000 & gp(1,:)/9.81>=1000)),2)-mean(entS(:,(gp(1,:)/9.81<=7000 & gp(1,:)/9.81>=5000)),2);
for i = 1:size(q,1)
satfrac(i) = trapz(pressure(11:end),q(i,11:end))./trapz(pressure(11:end),qvs(i,11:end));
end

precipitation = (int_S(1:length(radFlux))'-radFlux-SHflux(1:length(radFlux)))*(86400*1000)/(Lv*rhov);

anglecolors=interp1(angles,testmap,angle,[],'extrap');
mean_color=interp1(angles,testmap,angle_mean,[],'extrap');


% figure(2),
% ax1w=subplot(2,2,1);
% hold on;
% scatter(satfrac,II);
% ax2w=subplot(2,2,2);
% hold on;
% scatter(satfrac,GMS_tot);
% set(gca,'ylim',[-1 1]);
% ax3w=subplot(2,2,3);
% hold on;
% scatter(II,GMS_tot);
% set(gca,'ylim',[-1 1]);
% ax4w=subplot(2,2,4);
% hold on;
% scatter(satfrac(1:length(precipitation)),precipitation);
% hold on;






figure(2),
set(gcf,'units','pixels','position',[2088 217 1020 705]);
hold on;
scatter(SST_box,II,108,anglecolors,'d','filled');
hold on;
scatter(mean(SST_box),mean(II),324,mean_color,'*');
set(gca,'fontsize',16,'linewidth',2)
title('ERA5 Monthly mean SST vs II');
xlabel('SST (K)');
ylabel('Instability Index');


ERA5_west_precip;