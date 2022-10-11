%clear;  
%close all;

%addpath('./mexcdf/mexnc')
%addpath('./mexcdf/snctools/')
%addpath('./mexcdf/netcdf_toolbox/')

% 
% iloaddata=1;
% if(iloaddata==1)
%     [dfld,ifaTimes,ifaNz,ifaNt ]=read_fields_ifa;
%     [dadv,ifaTimes,ifaNz,ifaNt ]=read_advect_ifa;
% end
%addpath('../');
%[dfld,ifaTimes,ifaNz_tot,ifaNt_tot ]=read_dynamo_v1_aug2013('./lsf_v1_agu_2013/');

%ifaNt=301;
 
fname = './wrfinput_d01';


ninfo = nc_info(fname);
sstNx = 64; sstNy = 64;
sstNx = 128; sstNy = 128;
sstNx = 95; sstNy = 95;
sstNx = 256; sstNy = 256;

fout = ['dyna_force_sst_aug13_west.nc'] ;

nc_create_empty ( fout );
for ii=1:length(ninfo.Dimension)
    if(strcmp(ninfo.Dimension(ii).Name,'Time'));  % Special for wrf output
        nc_add_dimension ( fout, ninfo.Dimension(ii).Name, 0 );
    elseif(strcmp(ninfo.Dimension(ii).Name,'force_layers')); 
        nc_add_dimension ( fout, ninfo.Dimension(ii).Name, ifaNz );  
    elseif(strcmp(ninfo.Dimension(ii).Name,'west_east'));         
        nc_add_dimension ( fout, ninfo.Dimension(ii).Name, sstNx );  
    elseif(strcmp(ninfo.Dimension(ii).Name,'south_north')); 
        nc_add_dimension ( fout, ninfo.Dimension(ii).Name, sstNy );  
        
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


for ii=1:1 %length(ninfo.Dataset);
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
        %nc_varput ( fout, varstruct.Name, reshape(vardata,[1 size(vardata)]),zeros(1,length(varstruct.Dimension)),[1 size(vardata)] );
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
        ninfo.Dataset(ii).Attribute(jj).Name;
    end
 
    clear vardata %varstruct
end


 
 
 
%---------------------------

sst=repmat(302.497,[801,sstNy,sstNx]);


%load ../toga/misc.ifa
%sst=repmat(misc(:,6)+273.15,[1,sstNy,sstNx]);

varst.Name = 'SST'; 
varst.Nctype=5; varst.Unlimited=1; varst.Dimension={'Times',  'south_north', 'west_east'};
nc_addvar ( fout, varst )

nc_varput ( fout, varst.Name, sst, [0 0 0],[ifaNt sstNy sstNx] );
nc_attput ( fout, varst.Name, 'FieldType', 104);
nc_attput ( fout, varst.Name, 'MemoryOrder', 'XY ');
nc_attput ( fout, varst.Name, 'description', 'SEA SURFACE TEMPERATURE');
nc_attput ( fout, varst.Name, 'units', 'K');
nc_attput ( fout, varst.Name, '_FillValue', '-999');


 



nc_dump( fout )
newinfo=nc_info(fout);

