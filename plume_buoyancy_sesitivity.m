close all;
clear variables;

load('anglecolormap5.mat');


plotDirectory = '../WRF_OUT/Plot/';

 allDirectory = '../WRF_OUT/average/ew_experiment/';

 
dynNames={};
     
 dynPlot=[plotDirectory 'ew_experiment/switch/'];
 dynFiles = dir([allDirectory 'switch_radiation/wrfout_*']);
 
 for i = 1:length(dynFiles)
 dynNames{i}=dynFiles(i).name;
 end
 
dynstats=simple_wtg_plot_tropolimit([allDirectory 'switch_radiation/'],dynNames,dynPlot,1,0);
w1 = dynstats.omega{1};
tp1 = dynstats.tp{1};
hp1 = dynstats.hp{1};
h1 = dynstats.h{1};
hstar1 = dynstats.hstar{1};
t1 = dynstats.th{1};
qv1 = dynstats.qv{1};
w2 = dynstats.omega{4};
tp2 = dynstats.tp{4};
hp2 = dynstats.hp{4};
h2 = dynstats.h{4};
hstar2 = dynstats.hstar{4};
t2 = dynstats.th{4};
qv2 = dynstats.qv{4};
w3 = dynstats.omega_rad{5}';
tp3 = dynstats.tp{5};
hp3 = dynstats.hp{5};
h3 = dynstats.h{5};
hstar3 = dynstats.hstar{5};
t3 = dynstats.th{5};
qv3 = dynstats.qv{5};
w4 = dynstats.omega{3};
tp4 = dynstats.tp{3};
hp4 = dynstats.hp{3};
h4 = dynstats.h{3};
hstar4 = dynstats.hstar{3};
t4 = dynstats.th{3};
qv4 = dynstats.qv{3};
w5 = dynstats.omega{6};
tp5 = dynstats.tp{6};
hp5 = dynstats.hp{6};
h5 = dynstats.h{6};
hstar5 = dynstats.hstar{6};
t5 = dynstats.th{6};
qv5 = dynstats.qv{6};
w6 = dynstats.omega_rad{9}';
tp6 = dynstats.tp{9};
hp6 = dynstats.hp{9};
h6 = dynstats.h{9};
hstar6 = dynstats.hstar{9};
t6 = dynstats.th{9};
qv6 = dynstats.qv{9};
w7 = dynstats.omega{10};
tp7 = dynstats.tp{10};
hp7 = dynstats.hp{10};
h7 = dynstats.h{10};
hstar7 = dynstats.hstar{10};
t7 = dynstats.th{10};
qv7 = dynstats.qv{10};
w8 = dynstats.omega{8};
tp8 = dynstats.tp{8};
hp8 = dynstats.hp{8};
h8 = dynstats.h{8};
hstar8 = dynstats.hstar{8};
t8 = dynstats.th{8};
qv8 = dynstats.qv{8};
w9 = dynstats.omega{2};
tp9 = dynstats.tp{2};
hp9 = dynstats.hp{2};
h9 = dynstats.h{2};
hstar9 = dynstats.hstar{2};
t9 = dynstats.th{2};
qv9 = dynstats.qv{2};
w10 = dynstats.omega{7};
tp10 = dynstats.tp{7};
hp10 = dynstats.hp{7};
h10 = dynstats.h{7};
hstar10 = dynstats.hstar{7};
t10 = dynstats.th{7};
qv10 = dynstats.qv{7};

p=mean(dynstats.p,2);
z=dynstats.z{4};
z=mean((z(2:end,:)+z(1:end-1,:))/2,2);
qvs1=calc_qvstar(mean(t1(:,end-120:end),2),p);RH1 = mean(qv1(:,end-120:end),2)./qvs1;
qvs2=calc_qvstar(mean(t2(:,end-120:end),2),p);RH2 = mean(qv2(:,end-120:end),2)./qvs2;
qvs3=calc_qvstar(mean(t3(:,1),2),p);RH3 = mean(qv3(:,1),2)./qvs3;
qvs4=calc_qvstar(mean(t4(:,end-120:end),2),p);RH4 = mean(qv4(:,end-120:end),2)./qvs4;
qvs5=calc_qvstar(mean(t5(:,end-120:end),2),p);RH5 = mean(qv5(:,end-120:end),2)./qvs5;
qvs6=calc_qvstar(mean(t6(:,1),2),p);RH6 = mean(qv6(:,1),2)./qvs6;
qvs7=calc_qvstar(mean(t7(:,end-120:end),2),p);RH7 = mean(qv7(:,end-120:end),2)./qvs7;
qvs8=calc_qvstar(mean(t8(:,end-120:end),2),p);RH8 = mean(qv8(:,end-120:end),2)./qvs8;
qvs9=calc_qvstar(mean(t9(:,end-120:end),2),p);RH9 = mean(qv9(:,end-120:end),2)./qvs9;
qvs10=calc_qvstar(mean(t10(:,end-120:end),2),p);RH10 = mean(qv10(:,end-120:end),2)./qvs10;
[Se1,He1,HSATe1,Hpe1,Tpe1]= zero_plume_buoyancy_multi_plume(mean(t1(:,end-120:end),2),mean(qv1(:,end-120:end),2),z,p,4000);
ze1=Hpe1>=HSATe1;Tdiff1=Tpe1-mean(t1(:,end-120:end),2);
[Se2,He2,HSATe2,Hpe2,Tpe2]= zero_plume_buoyancy_multi_plume(mean(t2(:,end-120:end),2),mean(qv2(:,end-120:end),2),z,p,4000);
ze2=Hpe2>=HSATe2;Tdiff2=Tpe2-mean(t2(:,end-120:end),2);
[Se3,He3,HSATe3,Hpe3,Tpe3]= zero_plume_buoyancy_multi_plume(mean(t3(:,1),2),mean(qv3(:,1),2),z,p,4000);
ze3=Hpe3>=HSATe3;Tdiff3=Tpe3-mean(t3(:,end-120:end),2);
[Se4,He4,HSATe4,Hpe4,Tpe4]= zero_plume_buoyancy_multi_plume(mean(t4(:,end-120:end),2),mean(qv4(:,end-120:end),2),z,p,4000);
ze4=Hpe4>=HSATe4;Tdiff4=Tpe4-mean(t4(:,end-120:end),2);
[Se5,He5,HSATe5,Hpe5,Tpe5]= zero_plume_buoyancy_multi_plume(mean(t5(:,end-120:end),2),mean(qv5(:,end-120:end),2),z,p,4000);
ze5=Hpe5>=HSATe5;Tdiff5=Tpe5-mean(t5(:,end-120:end),2);
[Se6,He6,HSATe6,Hpe6,Tpe6]= zero_plume_buoyancy_multi_plume(mean(t6(:,1),2),mean(qv6(:,1),2),z,p,4000);
ze6=Hpe6>=HSATe6;Tdiff6=Tpe6-mean(t6(:,end-120:end),2);
[Se7,He7,HSATe7,Hpe7,Tpe7]= zero_plume_buoyancy_multi_plume(mean(t7(:,end-120:end),2),mean(qv7(:,end-120:end),2),z,p,4000);
ze7=Hpe7>=HSATe7;Tdiff7=Tpe7-mean(t7(:,end-120:end),2);
[Se8,He8,HSATe8,Hpe8,Tpe8]= zero_plume_buoyancy_multi_plume(mean(t8(:,end-120:end),2),mean(qv8(:,end-120:end),2),z,p,4000);
ze8=Hpe8>=HSATe8;Tdiff8=Tpe8-mean(t8(:,end-120:end),2);
[Se9,He9,HSATe9,Hpe9,Tpe9]= zero_plume_buoyancy_multi_plume(mean(t9(:,end-120:end),2),mean(qv9(:,end-120:end),2),z,p,4000);
ze9=Hpe9>=HSATe9;Tdiff9=Tpe9-mean(t9(:,end-120:end),2);
[Se10,He10,HSATe10,Hpe10,Tpe10]= zero_plume_buoyancy_multi_plume(mean(t10(:,end-120:end),2),mean(qv10(:,end-120:end),2),z,p,4000);
ze10=Hpe10>=HSATe10;Tdiff10=Tpe10-mean(t10(:,end-120:end),2);


e1e=-((mean(w1(end-120:end,:),1)+mean(w4(end-120:end,:),1)+mean(w9(end-120:end,:),1)+w3+w6)/5);
e2e=-((mean(w1(end-120:end,:),1)+mean(w2(end-120:end,:),1)+mean(w7(end-120:end,:),1)+mean(w4(end-120:end,:),1)+mean(w9(end-120:end,:),1))/5);
e3e=-((mean(w1(end-120:end,:),1)+mean(w2(end-120:end,:),1)+w3+mean(w8(end-120:end,:),1)+mean(w9(end-120:end,:),1))/5);
e4e=-((mean(w1(end-120:end,:),1)+mean(w10(end-120:end,:),1)+mean(w2(end-120:end,:),1)+mean(w4(end-120:end,:),1)+w3)/5);
e1w=-((mean(w5(end-120:end,:),1)+mean(w2(end-120:end,:),1)+mean(w7(end-120:end,:),1)+mean(w8(end-120:end,:),1)+mean(w10(end-120:end,:),1))/5);
e2w=-((mean(w5(end-120:end,:),1)+mean(w8(end-120:end,:),1)+mean(w10(end-120:end,:),1)+w3+w6)/5);
e3w=-((mean(w5(end-120:end,:),1)+mean(w7(end-120:end,:),1)+mean(w4(end-120:end,:),1)+mean(w10(end-120:end,:),1)+w6)/5);
e4w=-((mean(w5(end-120:end,:),1)+mean(w7(end-120:end,:),1)+mean(w9(end-120:end,:),1)+mean(w8(end-120:end,:),1)+w6)/5);

Temp1e=((mean(t1(:,end-120:end),2)+mean(t4(:,end-120:end),2)+mean(t9(:,end-120:end),2)+mean(t3(:,1),2)+mean(t6(:,1),2))/5);
Temp2e=((mean(t1(:,end-120:end),2)+mean(t2(:,end-120:end),2)+mean(t7(:,end-120:end),2)+mean(t4(:,end-120:end),2)+mean(t9(:,end-120:end),2))/5);
Temp3e=((mean(t1(:,end-120:end),2)+mean(t2(:,end-120:end),2)+mean(t8(:,end-120:end),2)+mean(t9(:,end-120:end),2)+mean(t3(:,1),2))/5);
Temp4e=((mean(t1(:,end-120:end),2)+mean(t10(:,end-120:end),2)+mean(t2(:,end-120:end),2)+mean(t4(:,end-120:end),2)+mean(t3(:,1),2))/5);
Temp1w=((mean(t5(:,end-120:end),2)+mean(t2(:,end-120:end),2)+mean(t7(:,end-120:end),2)+mean(t8(:,end-120:end),2)+mean(t10(:,end-120:end),2))/5);
Temp2w=((mean(t5(:,end-120:end),2)+mean(t8(:,end-120:end),2)+mean(t10(:,end-120:end),2)+mean(t3(:,1),2)+mean(t6(:,1),2))/5);
Temp3w=((mean(t5(:,end-120:end),2)+mean(t7(:,end-120:end),2)+mean(t4(:,end-120:end),2)+mean(t10(:,end-120:end),2)+mean(t6(:,1),2))/5);
Temp4w=((mean(t5(:,end-120:end),2)+mean(t7(:,end-120:end),2)+mean(t9(:,end-120:end),2)+mean(t8(:,end-120:end),2)+mean(t6(:,1),2))/5);

QVs1e=calc_qvstar(Temp1e,p);
QVs2e=calc_qvstar(Temp2e,p);
QVs3e=calc_qvstar(Temp3e,p);
QVs4e=calc_qvstar(Temp4e,p);
QVs1w=calc_qvstar(Temp1w,p);
QVs2w=calc_qvstar(Temp2w,p);
QVs3w=calc_qvstar(Temp3w,p);
QVs4w=calc_qvstar(Temp4w,p);



QV1e=RH1.*QVs1e;
QV2e=RH1.*QVs2e;
QV3e=RH1.*QVs3e;
QV4e=RH1.*QVs4e;
QV1w=RH5.*QVs1w;
QV2w=RH5.*QVs2w;
QV3w=RH5.*QVs3w;
QV4w=RH5.*QVs4w;

% Tdiff1e=((Tdiff1+Tdiff4+Tdiff9)/3);
% Tdiff2e=((Tdiff1+Tdiff2+Tdiff7+Tdiff4+Tdiff9)/5);
% Tdiff3e=((Tdiff1+Tdiff2+Tdiff8+Tdiff9)/4);
% Tdiff4e=((Tdiff1+Tdiff10+Tdiff2+Tdiff4)/4);
% Tdiff1w=((Tdiff5+Tdiff2+Tdiff7+Tdiff8+Tdiff10)/5);
% Tdiff2w=((Tdiff5+Tdiff8+Tdiff10)/3);
% Tdiff3w=((Tdiff5+Tdiff7+Tdiff4+Tdiff10)/4);
% Tdiff4w=((Tdiff5+Tdiff7+Tdiff9+Tdiff8)/4);
% 
% HSAT1e=((HSATe1+HSATe4+HSATe9)/3);
% HSAT2e=((HSATe1+HSATe2+HSATe7+HSATe4+HSATe9)/5);
% HSAT3e=((HSATe1+HSATe2+HSATe8+HSATe9)/4);
% HSAT4e=((HSATe1+HSATe10+HSATe2+HSATe4)/4);
% HSAT1w=((HSATe5+HSATe2+HSATe7+HSATe8+HSATe10)/5);
% HSAT2w=((HSATe5+HSATe8+HSATe10)/3);
% HSAT3w=((HSATe5+HSATe7+HSATe4+HSATe10)/4);
% HSAT4w=((HSATe5+HSATe7+HSATe9+HSATe8)/4);1
% 
% H1e=((He1+He4+He9)/3);
% H2e=((He1+He2+He7+He4+He9)/5);
% H3e=((He1+He2+He8+He9)/4);
% H4e=((He1+He10+He2+He4)/4);
% H1w=((He5+He2+He7+He8+He10)/5);
% H2w=((He5+He8+He10)/3);
% H3w=((He5+He7+He4+He10)/4);
% H4w=((He5+He7+He9+He8)/4);
% 
% Hp1e=((Hpe1+Hpe4+Hpe9)/3);
% Hp2e=((Hpe1+Hpe2+Hpe7+Hpe4+Hpe9)/5);
% Hp3e=((Hpe1+Hpe2+Hpe8+Hpe9)/4);
% Hp4e=((Hpe1+Hpe10+Hpe2+Hpe4)/4);
% Hp1w=((Hpe5+Hpe2+Hpe7+Hpe8+Hpe10)/5);
% Hp2w=((Hpe5+Hpe8+Hpe10)/3);
% Hp3w=((Hpe5+Hpe7+Hpe4+Hpe10)/4);
% Hp4w=((Hpe5+Hpe7+Hpe9+Hpe8)/4);

[S1e,H1e,HSAT1e,Hp1e,Tp1e]= zero_plume_buoyancy_multi_plume(Temp1e,QV1e,z,p,4000);
z1e=Hp1e>=HSAT1e;Tdiff1e=Tp1e-Temp1e;
[S2e,H2e,HSAT2e,Hp2e,Tp2e]= zero_plume_buoyancy_multi_plume(Temp2e,QV2e,z,p,4000);
z2e=Hp2e>=HSAT2e;Tdiff2e=Tp2e-Temp2e;
[S3e,H3e,HSAT3e,Hp3e,Tp3e]= zero_plume_buoyancy_multi_plume(Temp3e,QV3e,z,p,4000);
z3e=Hp3e>=HSAT3e;Tdiff3e=Tp3e-Temp3e;
[S4e,H4e,HSAT4e,Hp4e,Tp4e]= zero_plume_buoyancy_multi_plume(Temp4e,QV4e,z,p,4000);
z4e=Hp4e>=HSAT4e;Tdiff4e=Tp4e-Temp4e;
[S1w,H1w,HSAT1w,Hp1w,Tp1w]= zero_plume_buoyancy_multi_plume(Temp1w,QV1w,z,p,4000);
z1w=Hp1w>=HSAT1w;Tdiff1w=Tp1w-Temp1w;
[S2w,H2w,HSAT2w,Hp2w,Tp2w]= zero_plume_buoyancy_multi_plume(Temp2w,QV2w,z,p,4000);
z2w=Hp2w>=HSAT2w;Tdiff2w=Tp2w-Temp2w;
[S3w,H3w,HSAT3w,Hp3w,Tp3w]= zero_plume_buoyancy_multi_plume(Temp3w,QV3w,z,p,4000);
z3w=Hp3w>=HSAT3w;Tdiff3w=Tp3w-Temp3w;
[S4w,H4w,HSAT4w,Hp4w,Tp4w]= zero_plume_buoyancy_multi_plume(Temp4w,QV4w,z,p,4000);
z4w=Hp4w>=HSAT4w;Tdiff4w=Tp4w-Temp4w;
    


    figure('units','normalized','outerposition',[0.05 0.05 .7 .9]),
    subplot(1,2,1),
    plot(e1e-e1w,p,'k--',-e2e+e2w,p,'ko-',-e3e+e3w,p,'k:',-e4e+e4w,p,'k-.','linewidth',3)
    set(gca,'ylim',[10000 100000],'ydir','reverse','xlim',[-0.25 0.1],'fontsize',24);
    legend('SST','Stability','Qrad','Qhadv','location','southwest');
    xlabel('Vertical Velocity Difference');
    subplot(1,2,2),
    plot(Tdiff1e-Tdiff1w,p,'k--',-Tdiff2e+Tdiff2w,p,'ko-',-Tdiff3e+Tdiff3w,p,'k:',-Tdiff4e+Tdiff4w,p,'k-.','linewidth',3)
    set(gca,'ylim',[10000 100000],'ydir','reverse','xlim',[-2 2],'fontsize',24,'yticklabel',[],'xdir','reverse');
    xlabel('Buoyancy Difference');
    sgtitle('Vertical Motion and Buoyancy Sensitivity','fontsize',24);
    print(gcf,'-dpng',[plotDirectory 'Sensitivity_Buoyacy_diagram.png']);