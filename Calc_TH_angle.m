function [pcs_ERA,lds_ERA_interp] = Calc_TH_angle(omega,P)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if size(omega,1)<size(omega,2)
   omega=omega'; 
end
if min(size(P))==1
    Pmean=P;
elseif size(P,1)>size(P,2)
    Pmean=mean(P);
else
   Pmean=mean(P,2); 
end


load('ERA5_EOFS.mat');
weights_old=weights;

lds_ERA=lds./repmat(weights_old.^0.5,1,size(lds,2)); 
lds_ERA_interp=interp1(level,lds_ERA,Pmean,'linear','extrap');
lds_ERA_interp(Pmean<=level(1),:)=0;
omega(:,Pmean<=level(1))=0;
if size(omega,1)==size(lds_ERA_interp,1)
    omega=omega';
end
pcs_ERA=lds_ERA_interp\omega';



end

