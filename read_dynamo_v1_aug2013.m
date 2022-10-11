function [dout,Times,nlev,nt ] = read_dynamo_v1_aug2013(fpath)

if(nargin < 1)
    fpath = 'lsf_dec13_2012/';
    fpath = 'lsf_v1_agu_2013/';
end


%  format for data:
% for each 6 hour period there is one 1 with date/time information followed by 41 with data as a
% function of pressure
% 
% order of fields and their units are given below
% 
%  yy mm dd hh
%     p          z        u        v      omega     T       theta     wmr     rh        div      vort
%     hPa        m       m/s       m/s     mb/hr    C        K        g/kg     (%)     x10^6s-1  x10^5s-1
% 
% 
% to convert omega (mb/hr) to w (mm/s) use the following approximate eqn: omega/(-g*rho) = w
% 
% g = 9.8 m/s
% rho = p/RT can be computed from fields provided (for example near 1000 mb, rho  ~ 1.15 kg/m^3,
%                                                              near  100 mb, rho  ~ 0.18 kg/m^3)
% So that 1 mb/hr is ~2.5 mm/s  near 1000 mb, while 1 mb/hr near 100 mb is ~16 mm/s


file = [fpath '/fields.lsan_v1'];

fid = fopen(file, 'r');
Times = []; nlev = 40; nt=301; ncol=11;
pres=zeros(nt,nlev); z=pres; u=pres; v=pres; omega=pres; T=pres; theta=pres; qv=pres; div=pres;  rh=pres;
adate = zeros(nt,19);
for i=1:nt    
    [aa,dn]=fscanf(fid,'%3i',4);
    if(isempty(aa))        
        disp(['End of data: ' rdate]);
        break;
    end
    rdate = datestr(datenum(aa(1)+2000, aa(2), aa(3), aa(4), 0, 0),'yyyy-mm-dd_HH:MM:SS');
    Times=[ Times; rdate];       
    [data,dn]=fscanf(fid,'%f',ncol*nlev);
    dtmp=reshape(data,ncol, nlev);    
    dtmp(dtmp==-999)=NaN;
    
    adate(i,:)=rdate;
    pres(i,:)=dtmp(1,:);
    z(i,:)=dtmp(2,:); 
    T(i,:)=dtmp(6,:); 
    qv(i,:)=dtmp(8,:); 
    rh(i,:)=dtmp(9,:); 
    omega(i,:) = dtmp(5,:);
    theta(i,:)=dtmp(7,:);
    div(i,:) = dtmp(10,:);
    vor(i,:) = dtmp(11,:);
    u(i,:)=dtmp(3,:); v(i,:)=dtmp(4,:); 
end
fclose(fid);
dout.pres = pres;
dout.z = z; 
dout.omg = omega; 
dout.w = -omega*(100/3600)./(9.81*1e2*pres./(287*(T+273.15)));
dout.div=div; dout.vor=vor;
dout.theta=theta;
dout.u=u; dout.v=v;  
dout.T=T;   dout.qv=qv/1e3;  dout.rh=rh;
dout.date=adate;


% lsf_flds.ifa
% File "lsf_flds.ifa" contains 480 periods of six-hourly data. For each six hour period there are 42 lines of data.
% 
% Line 1 contains: [year, month, day, hour]
%                  written with format (4i3)
% Lines 2-42 contain: [p(mb), hT(C/s), vT(C/s), hq(gr/(gr*s)), vq(gr/(gr*s))]
%                           written with format (f8.2,1p,8e11.3)
%                           where hT - horizontal advection of T
%                           where vT -   vertical advection of T
%                           where hq - horizontal advection of q
%                           where vq -   vertical advection of q
% These horizontal and vertical advection terms were computed using centered differences as follows:
% horizontal advection of "f": h(f) = u*df/dx + v*df/dy
%                       where: dx = a cos(phi)*d(lambda)
%                              dy = a d(phi)
%                              phi    - latitiude
%                              lambda - longitude
% vertical advection of "f": v(f) = omega*df/dp

file = [fpath '/lsf.lsan_v1']; % horizontal and vertical advection of temperature and mositure
fid = fopen(file, 'r');
%Times = []; nlev = 41; nt=480;   ncol=5;
Times = []; nlev = 40; nt=301; ncol=5;
pres=zeros(nt,nlev);   ht=pres; vt=pres; hq=pres; vq=pres;  
for i=1:nt    
    [aa,dn]=fscanf(fid,'%3i',4);
    rdate = datestr(datenum(aa(1)+2000, aa(2), aa(3), aa(4), 0, 0),'yyyy-mm-dd_HH:MM:SS');
    Times=[ Times; rdate];        
    [data,dn]=fscanf(fid,'%f',ncol*nlev);
    dtmp=reshape(data,ncol, nlev);
    dtmp(dtmp==-999)=NaN;
    pres(i,:)=dtmp(1,:);
    ht(i,:)=dtmp(2,:); vt(i,:)=dtmp(3,:); 
    hq(i,:)=dtmp(4,:); vq(i,:)=dtmp(5,:); 
end
fclose(fid);
dout.pres2 = pres; 
dout.hT=ht; dout.vT=vt;  dout.vq=vq; dout.hq=hq; 


%-------------------------------------
load([fpath '/eopo.lsan_v1']);
dout.rain = eopo(:,6);
dout.evap = eopo(:,5);

%-------------------------------------
load([ fpath '/trmm_rain']);
dat2 = trmm_rain(:,1);
dat2 = datenum(num2str(dat2),'yymmddHH');
rr2 = trmm_rain(:,2);

dout.rain_trmm =rr2;
dout.date_trmm =dat2;
rain6_t = BlockMean(rr2',1, 2);
dout.rain6_trmm = rain6_t; 

raind_t = BlockMean(rr2',1, 8);
dout.raind_trmm = raind_t(1:76); 

%-------------------------------------

%---------------------------
file = [fpath '/nsa.dfluxes']; % daily flux from whoi oaflux
fid = fopen(file, 'r');
a=fgetl(fid);a=fgetl(fid);

[aa,dn]=fscanf(fid,'%f\n',92*10);
data=reshape(aa,[10 92])';
datn = datenum(data(:,1), data(:,2), data(:,3));
sst = data(:,9);
lh = data(:,4);
sh = data(:,5);

aa=repmat(sst(1:75),[1 4])'; bb=[aa(:); sst(76)];
dout.sst=bb;

aa=repmat(lh(1:75),[1 4])'; bb=[aa(:); lh(76)];
dout.lh=bb;


aa=repmat(sh(1:75),[1 4])'; bb=[aa(:); sh(76)];
dout.sh=bb;

return





%---------------------------
file = 'lsf_dec13_2012/nesa.dfluxes'; % daily flux from whoi oaflux
fid = fopen(file, 'r');
%Times = []; nlev = 41; nt=480;   ncol=5;
Times = []; nlev = 40; nt=301; ncol=5;
pres=zeros(nt,nlev);   ht=pres; vt=pres; hq=pres; vq=pres;  
for i=1:nt    
    [aa,dn]=fscanf(fid,'%3i',4);
    rdate = datestr(datenum(aa(1)+2000, aa(2), aa(3), aa(4), 0, 0),'yyyy-mm-dd_HH:MM:SS');
    Times=[ Times; rdate];        
    [data,dn]=fscanf(fid,'%f',ncol*nlev);
    dtmp=reshape(data,ncol, nlev);
    dtmp(dtmp==-999)=NaN;
    pres(i,:)=dtmp(1,:);
    ht(i,:)=dtmp(2,:); vt(i,:)=dtmp(3,:); 
    hq(i,:)=dtmp(4,:); vq(i,:)=dtmp(5,:); 
end
fclose(fid);
dout.pres2 = pres; 
dout.hT=ht; dout.vT=vt;  dout.vq=vq; dout.hq=hq; 







dout.hfx = misc_flds_ifa_v2(:,7)*28.9;
dout.lh = misc_flds_ifa_v2(:,8)*28.9;
dout.rain = misc_flds_ifa_v2(:,9);
dout.radcool = misc_flds_ifa_v2(:,10);


 
return;


% Line 1 contains: [year, month, day, hour]
%                   written with format (4i3)
% Lines 2-42 contain: [p(mb), divergence(1/s)*10e6, vertical p-velocity (mb/hr),
%                      Q1 (K/day), Q2(K/day)]
%                   written with format (f7.1,1x,4f8.2)
% The apparent heat source, Q1, and moisture sink, Q2 (Yanai et al. 1973)
% were computed as:
% Q1/cp = [dT/dt + h(T) + (p/po)**kappa * omega * d(theta)/dp]
% Q2/cp = -Lv/cp * [dq/dt + h(q) + v(q)]
%       where dt = 12 hours
%             po = 1000 mb
%             cp = 1004
%             Lv = 2.5e6
%              g = 9.8

file = 'deriv_flds.ifa_v2.1';
fid = fopen(file, 'r');
Times = []; nlev = 41; nt=480;  ncol=5;
pres=zeros(nt,nlev);   omega=pres; div=pres; q1=pres; q2=pres;  
for i=1:nt    
    [aa,dn]=fscanf(fid,'%3i',4);
    rdate = datestr(datenum(aa(1)+1900, aa(2), aa(3), aa(4), 0, 0),'yyyy-mm-dd_HH:MM:SS');
    Times=[ Times; rdate];        
    [data,dn]=fscanf(fid,'%f',ncol*nlev);
    dtmp=reshape(data,ncol, 41);
    pres(i,:)=dtmp(1,:);
    div(i,:)=dtmp(2,:); omega(i,:)=dtmp(3,:); 
    q1(i,:)=dtmp(4,:); q2(i,:)=dtmp(5,:); 
end
fclose(fid);
dout.pres2 = pres; 
dout.omega=omega; dout.div=div;  dout.q1=q1; dout.q2=q2; 




% misc_flds.ifa
% File "misc_flds.ifa" contains values of IFA averaged miscellaneous data (IR brightness temperature, sea surface temperature and sensible and latent heat fluxes, budget-derived rainfall and net-radiative heating) The sensible and latent heat fluxes represent the average from several buoys in the IFA. Once surface evaporation is known, rainfall rates can be computed from the moisture budget by integrating the equation for Q2 (shown above) from 1000 mb to 100 mb as follows:
% po = eo + 1./(g*Lv) * [integral(Q2*dp) from 1000mb to 100mb]
% 
% Line 1-478 contains: [year, month, day, hour, IFA average brightness temp (C), 
% IFA average SST(C), IFA average sensible and latent heat flux (mm/day), 
% rainfall (mm/day), net radiative heating rate (mm/day)]
%                      written with format (4i3,6f8.2)
                     
load misc_flds.ifa_v2.1
dout.sst = misc_flds_ifa_v2(:,6);
dout.hfx = misc_flds_ifa_v2(:,7)*28.9;
dout.lh = misc_flds_ifa_v2(:,8)*28.9;
dout.rain = misc_flds_ifa_v2(:,9);
dout.radcool = misc_flds_ifa_v2(:,10);



if(1==2)
    dfld.omega=dfld.omg;
    dname=fieldnames(dfld);
    for ifld=1:numel(dname)
        %disp(dname{ifld})
        eval(['tmp=dfld.' dname{ifld} ';'])
        if(size(tmp,2)>2)
            nanend=find(tmp(:,end-1:end)==-999);
            if(~isempty(nanend));
                disp(dname{ifld})  ;
                nanend,
                eval((['dfld.' dname{ifld} '(nanend)=-99999999;']))
            end;
        end
    end
end

    
return


