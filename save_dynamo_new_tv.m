%clear;  
close all;

addpath('/home/wangs/bin/mexcdf/mexnc')
addpath('/home/wangs/bin/mexcdf/snctools/')
addpath('/home/wangs/bin/mexcdf/netcdf_toolbox/')
addpath('../');


iloaddata=2; % processing data for 1: v1, 2: v2, 3: sunysb

if(iloaddata==1)% read data from ifa version 1
     

    
elseif(iloaddata==2) % read data from ifa version 2
    %[dfld,ifaTimes,ifaNz_tot,ifaNt_tot ]=read_fields_ifa_v2;
    [dfld,ifaTimes,ifaNz_tot,ifaNt_tot ]=read_dynamo_v1_aug2013('lsf_v1_agu_2013/');
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
    
    
    sst=zeros(ifaNt_tot,1); 
    hfx=dfld.sh;
    lhf=dfld.lh;   
    r_d          = 287.;
    cp           = 7.*r_d/2.;
    g            = 9.81;
    nanid = find(dfld.T==-99999999);
    dfld.theta=(dfld.T+273.15).*(1000./dfld.pres).^(r_d/cp);
    dfld.theta(nanid)=-99999999;
    
    dadv.pres=dfld.pres;
    dadv.hT=dfld.hT; dadv.vT=dfld.vT; 
    dadv.hq=dfld.hq; dadv.vq=dfld.vq; 
    fout = 'dyno_force_v2.nc';
    fout = 'dyno_force_aug13_tv.nc';
 
    
elseif(iloaddata==3) % read data from ifa sunysb
    
end


ifaNt=480;
ifaNz=39;
ifaNt = ifaNt_tot;
ifaNz = ifaNz_tot;

fname = './force_ideal.nc';


ninfo = nc_info(fname);

nc_create_empty ( fout );
for ii=1:length(ninfo.Dimension)
    if(strcmp(ninfo.Dimension(ii).Name,'Time'));  % Special for wrf output
        nc_add_dimension ( fout, ninfo.Dimension(ii).Name, 0 );
    elseif(strcmp(ninfo.Dimension(ii).Name,'force_layers')); 
        nc_add_dimension ( fout, ninfo.Dimension(ii).Name, ifaNz );  
    else; nc_add_dimension ( fout, ninfo.Dimension(ii).Name, ninfo.Dimension(ii).Length ); end
end

for ii=1:length(ninfo.Attribute)
    %nc_attput  ( fout, nc_global,ninfo.Attribute(ii).Name, ninfo.Attribute(ii).Value)
    value = ninfo.Attribute(ii).Value;
    if(ninfo.Attribute(ii).Nctype == 1); % nc_byte = 1;
    elseif(ninfo.Attribute(ii).Nctype == 2); value = char(ninfo.Attribute(ii).Value); % nc_char = 2;
    elseif(ninfo.Attribute(ii).Nctype == 3); % nc_short = 3;
    elseif(ninfo.Attribute(ii).Nctype == 4); value = int32(ninfo.Attribute(ii).Value); % nc_int = 4
    elseif(ninfo.Attribute(ii).Nctype == 5); value = single(ninfo.Attribute(ii).Value);
    elseif(ninfo.Attribute(ii).Nctype == 6); value = double(ninfo.Attribute(ii).Value);
    end
    try
    nc_attput(fout, nc_global, ninfo.Attribute(ii).Name,value);
    nc_attput  ( fout, nc_global, ninfo.Attribute(ii).Nctype, ninfo.Attribute(ii).Nctype)
    catch
    end
end


for ii=1:2 %length(ninfo.Dataset);
    disp([num2str(ii) 'th variable ' ninfo.Dataset(ii).Name]);
    %if (length(ninfo.Dataset(ii).Dimension)>3); continue; end;

    varstruct.Name = ninfo.Dataset(ii).Name;
    varstruct.Nctype = ninfo.Dataset(ii).Nctype;
    varstruct.Dimension = [ ninfo.Dataset(ii).Dimension ];
    %vardatao = nc_varget(fname,ninfo.Dataset(ii).Name,zeros(1,length(varstruct.Dimension)), [-1*ones(1,length(varstruct.Dimension))]);
    vardata= nc_varget(fname,ninfo.Dataset(ii).Name,zeros(1,length(varstruct.Dimension)));
    %vardata = mean(vardatao,1);
    
    switch(varstruct.Name)
        case 'Z_FORCE'
            %vardata=dfld.pres(1:ifaNt,:);
            vardata=dfld.z(1:ifaNt,1:ifaNz);
    end
    
    nc_addvar ( fout, varstruct )
    %if(nc_iscoordvar(fname,varstruct.Name))
    if(strcmp(varstruct.Name,'Times'));  % Special treatment for wrf output        
        %nc_varput ( fout, varstruct.Name, ifaTimes(1:7,:),[0 0],[7 19]);
        nc_varput ( fout, varstruct.Name, ifaTimes(1:ifaNt,:),[0 0],[ifaNt 19]);
    elseif (length(size(vardata)) < length(varstruct.Dimension) )
        nc_varput ( fout, varstruct.Name, reshape(vardata,[1 size(vardata)]),zeros(1,length(varstruct.Dimension)),[1 size(vardata)] );
    elseif (length(size(vardata)) == length(varstruct.Dimension) )
        nc_varput ( fout, varstruct.Name, vardata,zeros(1,length(varstruct.Dimension)),[ size(vardata)] );
    else
        nc_varput ( fout, varstruct.Name, reshape(vardata,[1 size(vardata)]),zeros(1,length(varstruct.Dimension)),[1] );%zeros(1,length(varstruct.Dimension));
    end

    
    for jj=1:length(ninfo.Dataset(ii).Attribute)
        value = ninfo.Dataset(ii).Attribute(jj).Value;
        if(ninfo.Dataset(ii).Attribute(jj).Nctype == 1); % nc_byte = 1;
        elseif(ninfo.Dataset(ii).Attribute(jj).Nctype == 2); value = char(ninfo.Dataset(ii).Attribute(jj).Value); % nc_char = 2;
        elseif(ninfo.Dataset(ii).Attribute(jj).Nctype == 3); % nc_short = 3;
        elseif(ninfo.Dataset(ii).Attribute(jj).Nctype == 4); value = int32(ninfo.Dataset(ii).Attribute(jj).Value); % nc_int = 4
        elseif(ninfo.Dataset(ii).Attribute(jj).Nctype == 5); value = single(ninfo.Dataset(ii).Attribute(jj).Value);
        elseif(ninfo.Dataset(ii).Attribute(jj).Nctype == 6); value = double(ninfo.Dataset(ii).Attribute(jj).Value);
        end
        if(~isempty(value))
          nc_attput ( fout, varstruct.Name,ninfo.Dataset(ii).Attribute(jj).Name, value);
        end
        ninfo.Dataset(ii).Attribute(jj).Name
    end
 
    clear vardata %varstruct
end


varst.Name = 'P_FORCE';
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time',  'force_layers'};
nc_addvar ( fout, varst )
nc_varput ( fout, varst.Name, dfld.pres, [0 0],[ifaNt ifaNz] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'Z');
nc_attput ( fout, varst.Name, 'description', 'Pressure at the force layer');
nc_attput ( fout, varst.Name, 'units', 'hPa');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);


varst.Name = 'TAU_LARGESCALE'; 
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time',  'force_layers'};
nc_addvar ( fout, varst )
tau_ls = 3600+zeros(ifaNt,ifaNz);
nc_varput ( fout, varst.Name, tau_ls, [0 0],[ifaNt ifaNz] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'Z');
nc_attput ( fout, varst.Name, 'description', 'largescale timescale');
nc_attput ( fout, varst.Name, 'units', 'K');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);


varst.Name = 'TAU_LARGESCALE_TEND'; 
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time',  'force_layers'};
nc_addvar ( fout, varst )
tau_ls = 0 + zeros(ifaNt,ifaNz);
nc_varput ( fout, varst.Name, tau_ls, [0 0],[ifaNt ifaNz] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'Z');
nc_attput ( fout, varst.Name, 'description', 'tendency largescale timescale ');
nc_attput ( fout, varst.Name, 'units', 'K');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);


varst.Name = 'TH_LARGESCALE'; 
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time',  'force_layers'};
nc_addvar ( fout, varst )
nc_varput ( fout, varst.Name, dfld.theta(1:ifaNt,1:ifaNz), [0 0],[ifaNt ifaNz] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'Z');
nc_attput ( fout, varst.Name, 'description', 'theta');
nc_attput ( fout, varst.Name, 'units', 'K');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);


varst.Name = 'TMK_LARGESCALE'; 
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time',  'force_layers'};
nc_addvar ( fout, varst )
nc_varput ( fout, varst.Name, dfld.T(1:ifaNt,1:ifaNz)+273.15, [0 0],[ifaNt ifaNz] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'Z');
nc_attput ( fout, varst.Name, 'description', 'temperature');
nc_attput ( fout, varst.Name, 'units', 'K');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);


varst.Name = 'TH_LARGESCALE_TEND'; 
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time',  'force_layers'};
nc_addvar ( fout, varst );
% tadv = (dadv.hT(1:ifaNt,1:ifaNz)+dadv.vT(1:ifaNt,1:ifaNz));
% thetaadv = tadv.*(1000./dadv.pres(1:ifaNt,1:ifaNz)).^(r_d/cp);
% nc_varput ( fout, varst.Name, thetaadv , [0 0],[ifaNt ifaNz] );
tadv = (dadv.hT+dadv.vT);
%www= -dfld.omega*100/3600./(dfld.pres*1e2/r_d./(dfld.T+273.15))./9.81;
tmkv=(dfld.T+273.15).*(1+0.608*dfld.qv*1e-3);
www= -dfld.omega*100/3600./(dfld.pres*1e2/r_d./tmkv)./9.81;
thetaadv = tadv.*(1000./dadv.pres).^(r_d/cp);
thetaadv2 = (g/cp.*www + tadv)  .*(1000./dadv.pres).^(r_d/cp);


thetaadv2 = (dadv.hT)  .*(1000./dadv.pres).^(r_d/cp);
nc_varput ( fout, varst.Name, -thetaadv2(1:ifaNt,1:ifaNz) , [0 0],[ifaNt ifaNz] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'Z');
nc_attput ( fout, varst.Name, 'description', 'horizontal and vertical temperature tendency');
nc_attput ( fout, varst.Name, 'units', 'K s-1');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);

  

varst.Name = 'U_LARGESCALE'; 
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time',  'force_layers'};
nc_addvar ( fout, varst )
nc_varput ( fout, varst.Name, dfld.u(1:ifaNt,1:ifaNz), [0 0],[ifaNt ifaNz] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'Z');
nc_attput ( fout, varst.Name, 'description', 'U');
nc_attput ( fout, varst.Name, 'units', 'm s_1');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);
 

 

varst.Name = 'V_LARGESCALE'; 
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time',  'force_layers'};
nc_addvar ( fout, varst )
nc_varput ( fout, varst.Name, dfld.v(1:ifaNt,1:ifaNz), [0 0],[ifaNt ifaNz] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'Z');
nc_attput ( fout, varst.Name, 'description', 'V');
nc_attput ( fout, varst.Name, 'units', 'm s_1');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);
 
 


varst.Name = 'W_SUBS'; 
%www= -dfld.omega*100/3600./(dfld.pres*1e2/r_d./(dfld.T+273.15))./9.81;
tmkv=(dfld.T+273.15).*(1+0.608*dfld.qv*1e-3);
www= -dfld.omega*100/3600./(dfld.pres*1e2/r_d./tmkv)./9.81;
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time',  'force_layers'};
nc_addvar ( fout, varst )
nc_varput ( fout, varst.Name, www(1:ifaNt,1:ifaNz), [0 0],[ifaNt ifaNz] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'Z');
nc_attput ( fout, varst.Name, 'description', 'large-scale vertical velocity');
nc_attput ( fout, varst.Name, 'units', 'm s-1');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);
 


% varst.Name = 'W_SUBS_TEND'; 
% varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time',  'force_layers'};
% nc_addvar ( fout, varst )
% nc_varput ( fout, varst.Name,  dadv.hv(1:ifaNt,:)+dadv.vv(1:ifaNt,:) , [0 0],[ifaNt ifaNz] );
% nc_attput ( fout, varst.Name, 'FieldType', 104);
% nc_attput ( fout, varst.Name, 'MemoryOrder', 'Z');
% nc_attput ( fout, varst.Name, 'description', 'tendency large-scale vertical velocity');
% nc_attput ( fout, varst.Name, 'units', 'm s-2');
% nc_attput ( fout, varst.Name, '_FillValue', -99999999);


dfld.qv=dfld.qv*1e3;
dadv.hq=dadv.hq/1e3;
dadv.vq=dadv.vq/1e3;

varst.Name = 'QV_LARGESCALE'; 
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time',  'force_layers'};
nc_addvar ( fout, varst )
nc_varput ( fout, varst.Name, dfld.qv(1:ifaNt,1:ifaNz)/1e3, [0 0],[ifaNt ifaNz] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'Z');
nc_attput ( fout, varst.Name, 'description', 'Qv');
nc_attput ( fout, varst.Name, 'units', 'kg kg-1');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);
 


varst.Name = 'QV_LARGESCALE_TEND'; 
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time',  'force_layers'};
nc_addvar ( fout, varst )
nc_varput ( fout, varst.Name, (-dadv.hq(1:ifaNt,1:ifaNz)), [0 0],[ifaNt ifaNz] );
%nc_varput ( fout, varst.Name, (-dadv.hq(1:ifaNt,1:ifaNz)-dadv.vq(1:ifaNt,1:ifaNz)), [0 0],[ifaNt ifaNz] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'Z');
nc_attput ( fout, varst.Name, 'description', 'large scale qv tendency');
nc_attput ( fout, varst.Name, 'units', 'kg kg-1 s-1');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);




varst.Name = 'HFX_FORCE';
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time'};
nc_addvar ( fout, varst )
nc_varput ( fout, varst.Name, hfx, [0],[ifaNt] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'XY');
nc_attput ( fout, varst.Name, 'description', 'SCM ideal surface sensible heat flux');
nc_attput ( fout, varst.Name, 'units', 'W m-2');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);


varst.Name = 'LH_FORCE';
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time'};
nc_addvar ( fout, varst )
nc_varput ( fout, varst.Name, lhf, [0],[ifaNt] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'T');
nc_attput ( fout, varst.Name, 'description', 'SCM ideal surface latent heat flux');
nc_attput ( fout, varst.Name, 'units', 'W m-2');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);


varst.Name = 'TSK_FORCE';
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time'};
nc_addvar ( fout, varst )
nc_varput ( fout, varst.Name, sst, [0],[ifaNt] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'T');
nc_attput ( fout, varst.Name, 'description', 'SCM ideal surface skin temperature');
nc_attput ( fout, varst.Name, 'units', 'K');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);

varst.Name = 'HFX_FORCE_TEND';
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time'};
nc_addvar ( fout, varst )
nc_varput ( fout, varst.Name, hfx*0.0, [0],[ifaNt] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'T');
nc_attput ( fout, varst.Name, 'description', 'SCM ideal surface sensible heat flux tendency ');
nc_attput ( fout, varst.Name, 'units', 'W m-2 s-1');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);


varst.Name = 'LH_FORCE_TEND';
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time'};
nc_addvar ( fout, varst )
nc_varput ( fout, varst.Name, lhf*0.0, [0],[ifaNt] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'T');
nc_attput ( fout, varst.Name, 'description', 'SCM ideal surface latent heat flux tendency ');
nc_attput ( fout, varst.Name, 'units', 'W m-2 s-1');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);



varst.Name = 'TSK_FORCE_TEND';
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time'};
nc_addvar ( fout, varst )
nc_varput ( fout, varst.Name, sst*0.0, [0],[ifaNt] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'T');
nc_attput ( fout, varst.Name, 'description', 'SCM ideal surface skin temperature tendency');
nc_attput ( fout, varst.Name, 'units', 'K s-1');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);

varst.Name = 'RAIN_B';
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time'};
nc_addvar ( fout, varst )
nc_varput ( fout, varst.Name, dfld.rain(1:ifaNt), [0],[ifaNt] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'T');
nc_attput ( fout, varst.Name, 'description', 'Budget derived rainfall ');
nc_attput ( fout, varst.Name, 'units', 'W m-2 s-1');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);

if(1==2)
varst.Name = 'RAIN_T';
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Time'};
nc_addvar ( fout, varst )
rout = BlockMean(dfld.rain_trmm(1:end-6),2,1);
nc_varput ( fout, varst.Name, rout, [0],[ifaNt] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'T');
nc_attput ( fout, varst.Name, 'description', 'RRMM rainfall ');
nc_attput ( fout, varst.Name, 'units', 'W m-2 s-1');
nc_attput ( fout, varst.Name, '_FillValue', -99999999);
end


nc_dump( fout )
newinfo=nc_info(fout);

