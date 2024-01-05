clear variables;
%close all;
load('anglecolormap2.mat');
testmap=[anglemap(33:end,:); anglemap(1:32,:)];
angles=linspace(-180,180,64);

pr=100000;



allDirectory = '../WRF_OUT/average/ew_experiment/';
 
%dynamoNames={};
     
dynamoFiles = dir([allDirectory 'dynamo/wrfout_*']);

 dynamofullname=[allDirectory 'dynamo/' dynamoFiles(1).name];
 th_comp= mean(nc_varget(dynamofullname,'TH_LARGESCALE'));

dynamoFiles = dynamoFiles(2:end-1);
dynamoTherms=zeros(length(dynamoFiles),7);


%1-ii,2-DCIN,3-SF,4-gms,5-angle_mean,6-rr,7-hfx
 for i = 1:length(dynamoFiles)
 dynamoTherms(i,:)=calc_thermodynamics([allDirectory 'dynamo/' dynamoFiles(i).name],1);
 end
 


T=zeros(length(dynamoFiles),59); 
TH=zeros(length(dynamoFiles),59);
qv=zeros(length(dynamoFiles),59);
qvs=zeros(length(dynamoFiles),59);
z=zeros(length(dynamoFiles),59);
p=zeros(length(dynamoFiles),59);
tp1 = zeros(length(dynamoFiles),59);
tp2 = zeros(length(dynamoFiles),59);
zindex=zeros(length(dynamoFiles),2);
w=zeros(length(dynamoFiles),59);
t_start=zeros(length(dynamoFiles),40);th_start=zeros(length(dynamoFiles),40);

z_start=mean(nc_varget([allDirectory 'dynamo/' dynamoFiles(i).name],'Z_FORCE'));



 for i = 1:length(dynamoFiles)
     
 %dynamoNames{i}=dynamoFiles(i).name;
 dynamofullname=[allDirectory 'dynamo/' dynamoFiles(i).name];
 temp=nc_varget(dynamofullname,'T')+300;
 TH(i,:)=mean(temp(end-120:end,:));
 th_start(i,:)=mean(nc_varget(dynamofullname,'TH_LARGESCALE'));
 temp=nc_varget(dynamofullname,'QVAPOR');
 qv(i,:)=mean(temp(end-120:end,:));
 
 temp=midData(nc_varget(dynamofullname,'PHB')'+nc_varget(dynamofullname,'PH')')';
 z(i,:)=mean(temp(end-120:end,:))/9.81;
 temp=nc_varget(dynamofullname,'PB')+nc_varget(dynamofullname,'P');
 p(i,:)=mean(temp(end-120:end,:));
 temp=nc_varget(dynamofullname,'OMEGA_WTG');
 w(i,:)=mean(temp(end-120:end,:));    
% cwv(i)=trapz(p(i,:)',-qv(i,:)'/9.81); %units in mm
 T(i,:)=TH(i,:).*(p(i,:)./pr).^(2/7);
 %t_start(i,:)=th_start(i,:).*(p(i,:)./pr).^(2/7);
 
 qvs(i,:) = calc_qvstar(T(i,:),p(i,:));
 [tp1temp,tp2temp] = zero_plume_buoyancy_mse(T(i,:),qv(i,:),z(i,:),p(i,:));%,w(i,:)); 

     
     tp1(i,:)=tp1temp;
     tp2(i,:)=tp2temp;
     
 end
 
 p_start=interp1(z(1,:),p(1,:),z_start,[],'extrap');
 
      RH=qv./qvs;

 
sf=trapz(p(1,:)',qv')./trapz(p(1,:)',qvs');

 %precip_equil = [1,2,3,5];
lgn_entries =  {'interactive';'stable';'unstable'};%;'T=-40C';'T=-10C';'Z=500m'};
%lgn_entries = {'1';'2';'3';'4';'5'};
 
figure,
set(gcf,'position',[73 1 1725 960]);
subplot(1,4,2),
plot(w,p(1,:),'linewidth',3);
% yline(11075,'--');
% yline(6795,'--');
% yline(500,'--');
ylim([10000 100000]);
set(gca,'ydir','reverse','fontsize',16,'linewidth',3)
xlabel('Vertical Velocity (pa/s)');

subplot(1,4,1)
plot(th_start(1,:)-th_start(1,:),p_start,'linewidth',3)
hold on;
plot(th_start(2:3,:)-th_comp,p_start,'linewidth',3)
set(gca,'fontsize',16,'linewidth',3)
% yline(11075,'--');
% yline(6795,'--');
% yline(500,'--');
ylim([10000 100000]);
xlim([-3 3]);
set(gca,'ylim',[10000 100000],'ydir','reverse')
xlabel({'Environmental Temp Anomaly (k)'});

subplot(1,4,3),
plot(RH*100,p(1,:),'linewidth',3)
set(gca,'fontsize',16,'linewidth',3)
% yline(11075,'--');
% yline(6795,'--');
% yline(500,'--');
ylim([10000 100000]);
xlim([0 100]);
set(gca,'ylim',[10000 100000],'ydir','reverse');
xlabel({'Relative Humidity (%)'});
%legend(strcat('SF=',num2str(sf(1),'%4.2f')),strcat('SF=',num2str(sf(2),'%4.2f')),strcat('SF=',num2str(sf(3),'%4.2f')),'location','southwest');
%lgnd1=legend(lgn_entries);


subplot(1,4,4)
plot(tp1-T,p(1,:),'linewidth',3)
set(gca,'fontsize',16,'linewidth',3)
% yline(11075,'--');
% yline(6795,'--');
% yline(500,'--');
ylim([10000 100000]);
set(gca,'ylim',[10000 100000],'ydir','reverse')
xlim([-5 7.5]);
xlabel({'Plume Temp difference (k)'});
lgnd2=legend(lgn_entries);
set(lgnd2,'position',[.135 .82 .10 .10]);


sgtitle('Dynamo WTG simulations','fontsize',24);
saveas(gcf,'../posterplots/Dynamo_Buoyancy_anomalies_qv_normal_entrain.png');

%%

