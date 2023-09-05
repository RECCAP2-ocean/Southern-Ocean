dir='my_path/';

prefix=[dir 'RECCAP/'];

models={'ORCA025-GEOMAR','CESM-ETHZ','FESOM_REcoM_LR','OCIM','NorESM-OC1.2','MRI-ESM2-1','CNRM-ESM2-1',...
    'ORCA1-LIM3-PISCES','planktom12','MPIOM-HAMOCC','EC-Earth3','ROMS-SouthernOcean-ETHZ','CCSM-WHOI','OCIM-CTL'};
mod_id={'ORCA025_GEOMAR','CESM_ETHZ','FESOM_REcoM_LR','OCIM_v2021','NorESM_OC12','MRI_ESM2_1','CNRM_ESM2_1',...
    'ORCA1_LIM3_PISCES','planktom12','MPIOM_HAMOCC','EC_Earth3','ROMS_SouthernOcean_ETHZ','CCSM_WHOI','OCIM_v2014'};
ver={'v20210804','v20211122','v20211119','v20210511','v20211125','v20220502','v20211208',...
    'v20211215','v20220404','v20220110','v20220323','v20220630','v20211125','v20210607'};


%%
%#######################################################
% save ancillary
%#######################################################

for m=1:numel(models)

    tmp=ncread([prefix models{m} '/volume_' models{m} '_1_gr_' ver{m} '.nc'],'volume');
    eval(['volume_',mod_id{m},'=tmp;']);
    tmp=ncread([prefix models{m} '/mask_vol_' models{m} '_1_gr_' ver{m} '.nc'],'mask_vol');
    eval(['mask_vol_',mod_id{m},'=tmp;']);
    tmp=ncread([prefix models{m} '/mask_vol_' models{m} '_1_gr_' ver{m} '.nc'],'depth');
    eval(['depth_',mod_id{m},'=tmp;']);
    tmp=ncread([prefix models{m} '/mask_vol_' models{m} '_1_gr_' ver{m} '.nc'],'lat');
    eval(['lat_',mod_id{m},'=tmp;']);
    tmp=ncread([prefix models{m} '/mask_vol_' models{m} '_1_gr_' ver{m} '.nc'],'lon');
    eval(['lon_',mod_id{m},'=tmp;']);
    tmp=ncread([prefix models{m} '/area_' models{m} '_1_gr_' ver{m} '.nc'],'area');
    eval(['area_',mod_id{m},'=tmp;']);
    
    if strcmp(mod_id{m},'CCSM_WHOI')
        mask_vol_CCSM_WHOI(mask_vol_CCSM_WHOI==1)=2;
        mask_vol_CCSM_WHOI(mask_vol_CCSM_WHOI==0)=1;
        mask_vol_CCSM_WHOI(mask_vol_CCSM_WHOI==2)=0;
        eval(['mask_sfc_',mod_id{m},'=squeeze(mask_vol_',mod_id{m},'(:,:,1));']);
        
    else
        tmp=ncread([prefix models{m} '/mask_sfc_' models{m} '_1_gr_' ver{m} '.nc'],'mask_sfc');
        eval(['mask_sfc_',mod_id{m},'=tmp;']);
    end
    
    
    eval(['area_',mod_id{m},'=double(area_',mod_id{m},');']);
    eval(['volume_',mod_id{m},'=double(volume_',mod_id{m},');']);
    eval(['mask_sfc_',mod_id{m},'=double(mask_sfc_',mod_id{m},');']);
    if strcmp(mod_id{m},'CESM-ETHZ')
        eval(['mask_vol_',mod_id{m},'=-double(mask_vol_',mod_id{m},');']);
    else
        eval(['mask_vol_',mod_id{m},'=-double(mask_vol_',mod_id{m},');']);
    end
    eval(['depth_',mod_id{m},'=double(depth_',mod_id{m},');']);
    
    % compute dz
    load([prefix models{m} '/depth_' models{m} '_1_gr_1980-2018_' ver{m} '.mat']);
    
    eval(['dz_',mod_id{m},'=NaN(size(depth_',mod_id{m},'));']);
    
    eval(['dz_',mod_id{m},'(1)=depth_',mod_id{m},'(1)*2;']);
    eval(['depthb=dz_',mod_id{m},'(1);']);
    eval(['dummy=numel(depth_',mod_id{m},');']);
    
    for z=2:dummy
        
        eval(['dz_',mod_id{m},'(z)=(depth_',mod_id{m},'(z)- depthb)*2;'])
        eval(['depthb=depthb+dz_',mod_id{m},'(z);']);
    end
    
    save([prefix models{m} '/dz_' models{m} '_1_gr_1980-2018_' ver{m} '.mat'],['dz_' mod_id{m}]);
    
    % save in matfiles
    save([prefix models{m} '/lat_' models{m} '_1_gr_1980-2018_' ver{m} '.mat'],['lat_' mod_id{m}]);
    save([prefix models{m} '/lon_' models{m} '_1_gr_1980-2018_' ver{m} '.mat'],['lon_' mod_id{m}]);
    save([prefix models{m} '/depth_' models{m} '_1_gr_1980-2018_' ver{m} '.mat'],['depth_' mod_id{m}]);
    save([prefix models{m} '/volume_' models{m} '_1_gr_1980-2018_' ver{m} '.mat'],['volume_' mod_id{m}]);
    save([prefix models{m} '/mask_vol_' models{m} '_1_gr_1980-2018_' ver{m} '.mat'],['mask_vol_' mod_id{m}]);
    save([prefix models{m} '/area_' models{m} '_1_gr_1980-2018_' ver{m} '.mat'],['area_' mod_id{m}]);
    save([prefix models{m} '/mask_sfc_' models{m} '_1_gr_1980-2018_' ver{m} '.mat'],['mask_sfc_' mod_id{m}]);
end

% compute 1x1 grid spacing
lat_out = [-89.5:1:89.5];
lon_out = [-179.5:1:179.5];

nm_orca=1.853315; % nautical mile in NEMO (computed in this way):
%nmax(sum(dx'))/360/60 % divided by (360 degrees * 60 minutes)
%nmax(sum(dx))/pi/2 % Earth radius in m

%radius=nm_orca*360*60/2/pi; % the corresponding Earth radius is 6.3712289

for i=1:numel(lon_out)
    for j=1:numel(lat_out)
        [dist_dx,phaseangle] = sw_dist( [lat_out(j) lat_out(j)],[lon_out(i)-0.5 lon_out(i)+0.5],'km');
        [dist_dy,phaseangle] = sw_dist( [lat_out(j)-0.5 lat_out(j)+0.5],[lon_out(i) lon_out(i)],'km');
        dx(j,i)=dist_dx*1e3*nm_orca/1.8520;  % scale by the nautical mile used in NEMO
        dy(j,i)=dist_dy*1e3*nm_orca/1.8520;  % scale by the nautical mile used in NEMO
    end
end

save([prefix '1x1_dx_dy.mat'],'dx','dy')


%%
%#######################################################
% save C_ant all years
%#######################################################
var='dissic'; % in mol/m3
years=[1980:2018];
% y=15 --> 1994
% y=28 --> 2007

for m=1:numel(models)

    models(m)
    if ~strcmp(mod_id{m},'CCSM_WHOI') % CCSM has monthly resolution
        tmpA=ncread([prefix models{m} '/' var '_' models{m} '_A_1_gr_1980-2018_' ver{m} '.nc'],var);
        tmpB=ncread([prefix models{m} '/' var '_' models{m} '_B_1_gr_1980-2018_' ver{m} '.nc'],var);
        
        
        if strcmp(mod_id{m},'ROMS_SouthernOcean_ETHZ') ||  strcmp(mod_id{m},'OCIM_v2014')  %ROMS doesn't have C
            
            tmpC=NaN(size(tmpB)); 
        else
            tmpC=ncread([prefix models{m} '/' var '_' models{m} '_C_1_gr_1980-2018_' ver{m} '.nc'],var);
            
        end
        
        if strcmp(mod_id{m},'OCIM_v2021') ||  strcmp(mod_id{m},'OCIM_v2014') 
            tmpD=NaN(size(tmpB));
        else
            tmpD=ncread([prefix models{m} '/' var '_' models{m} '_D_1_gr_1980-2018_' ver{m} '.nc'],var);         
        end
        
    else %CCSM-WHOI
        
        tmpA=ncread([prefix models{m} '/' var '_' models{m} '_A_1_gr_1980-2017_' ver{m} '.nc'],var);
        tmpA=reshape(tmpA,[360 180 size(tmpA,3) 12 size(tmpA,4)/12]);
        tmpA=squeeze(nanmean(tmpA,4));        
        
        tmpB=ncread([prefix models{m} '/' var '_' models{m} '_B_1_gr_1980-2017_' ver{m} '.nc'],var);
        tmpB=reshape(tmpB,[360 180 size(tmpB,3) 12 size(tmpB,4)/12]);
        tmpB=squeeze(nanmean(tmpB,4));
        
        tmpC=ncread([prefix models{m} '/' var '_' models{m} '_C_1_gr_1980-2017_' ver{m} '.nc'],var);
        tmpC=reshape(tmpC,[360 180 size(tmpC,3) 12 size(tmpC,4)/12]);
        tmpC=squeeze(nanmean(tmpC,4));
        
        tmpD=ncread([prefix models{m} '/' var '_' models{m} '_D_1_gr_1980-2017_' ver{m} '.nc'],var);
        tmpD=reshape(tmpD,[360 180 size(tmpD,3) 12 size(tmpD,4)/12]);
        tmpD=squeeze(nanmean(tmpD,4));
        
    end  
       
        for y=15:28
        years(y)
        varA=squeeze(tmpA(:,:,:,y));
        varB=squeeze(tmpB(:,:,:,y));
        varC=squeeze(tmpC(:,:,:,y));
        varD=squeeze(tmpD(:,:,:,y));
        
        if  strcmp(mod_id{m},'OCIM_v2014')              
            varC=varA; %they are equivalent        
        end
        
        eval(['cant_',mod_id{m},'_AminusD=varA-varD;']);
        eval(['cant_',mod_id{m},'_CminusB=varC-varB;']);
        
      
        save([prefix models{m} '/Cant_years/cant_' mod_id{m} '_' num2str(years(y)) '.mat'],...
            ['cant_' mod_id{m} '_AminusD'],['cant_' mod_id{m} '_CminusB']);
        end
    
end


%%
%#######################################################
% save 1994-2007 of surface variables
%#######################################################
vecm=[169:336]; % 1994:2007 (for monthly data)
vecy=[15:28]; % 1994:2007 (for yearly data)

for m=1:numel(models)

    models(m)
    
    if strcmp(mod_id{m},'OCIM_v2021') ||  strcmp(mod_id{m},'MPIOM_HAMOCC') || strcmp(mod_id{m},'CNRM_ESM2_1')
         exper='B';
    else
        exper='A';
    end
    
    % annual sigma and salinity
    
    tmp1=ncread([prefix models{m} '/tos_' models{m} '_' exper '_1_gr_1980-2018_' ver{m} '.nc'],'tos');
    tmp2=ncread([prefix models{m} '/sos_' models{m} '_' exper '_1_gr_1980-2018_' ver{m} '.nc'],'sos');
    tmp1=squeeze(tmp1);
    tmp2=squeeze(tmp2);
    sal=squeeze(nanmean(tmp2(:,:,vecm),3));
    temp=squeeze(nanmean(tmp1(:,:,vecm),3));
    
    eval(['sigma0_',mod_id{m},'=sw_dens0(sal,temp);'])
    
    save([prefix models{m} '/sigma0_' models{m} '_clim_' ver{m} '_94_07.mat'],['sigma0_',mod_id{m}])
    
    eval(['sos_',mod_id{m},'=sal;'])
    
    save([prefix models{m} '/sos_' models{m} '_clim_' ver{m} '_94_07.mat'],['sos_',mod_id{m}])
    
    if ~strcmp(mod_id{m},'OCIM_v2021')
        
        % September sea ice
        tmp=ncread([prefix models{m} '/fice_' models{m} '_' exper '_1_gr_1980-2018_' ver{m} '.nc'],'fice');
        tmp=squeeze(tmp);
        tmp=reshape(tmp,[size(tmp,1) size(tmp,2) 12 size(tmp,3)/12]);
        tmp=squeeze(tmp(:,:,9,:));
        fice=squeeze(nanmean(tmp(:,:,vecy),3));
        eval(['fice_sep_',mod_id{m},'=fice;'])
        save([prefix models{m} '/fice_sep_' models{m} '_clim_' ver{m} '_94_07.mat'],['fice_sep_',mod_id{m}])
        
        %September MLD (user-defined and fixed threshold)
        tmp=ncread([prefix models{m} '/mld_' models{m} '_' exper '_1_gr_1980-2018_' ver{m} '.nc'],'mld');
        tmp=squeeze(tmp);
        tmp=reshape(tmp,[size(tmp,1) size(tmp,2) 12 size(tmp,3)/12]);
        tmp=squeeze(tmp(:,:,9,:));
        mld=squeeze(nanmean(tmp(:,:,vecy),3));
        
        eval(['mld_sep_fthr_',mod_id{m},'=mld;'])
        
        save([prefix models{m} '/mld_sep_fix_thr_' models{m} '_clim_' ver{m} '_94_07.mat'],['mld_sep_fthr_',mod_id{m}])
    end
    
    % Cant flux
    
    if strcmp(mod_id{m},'MPIOM_HAMOCC') || strcmp(mod_id{m},'OCIM_v2021') || strcmp(mod_id{m},'CNRM_ESM2_1')
        exp_const='B';
        exp_var='C';
    else
        exp_const='D';
        exp_var='A';
    end
    tmp=ncread([prefix models{m} '/fgco2_' models{m} '_' exp_var '_1_gr_1980-2018_' ver{m} '.nc'],'fgco2');
    tmp=squeeze(tmp);
    fgco2_var=squeeze(nanmean(tmp(:,:,vecm),3)); %1994-2007
    tmp=ncread([prefix models{m} '/fgco2_' models{m} '_' exp_const '_1_gr_1980-2018_' ver{m} '.nc'],'fgco2');
    tmp=squeeze(tmp);
    fgco2_const=squeeze(nanmean(tmp(:,:,vecm),3)); %1994-2007
        
    eval(['cantfl_',mod_id{m},'=(fgco2_var-fgco2_const);'])
    
    if strcmp(mod_id{m},'CCSM_WHOI')
        eval(['cantfl_',mod_id{m},'=-cantfl_',mod_id{m},';'])
        
    end
    
    save([prefix models{m} '/cantflux_' models{m} '_clim_' ver{m} '_94_07.mat'],['cantfl_',mod_id{m}])
    
    % MLD annual diagnostic
    if ~strcmp(mod_id{m},'OCIM_v2021')
        if strcmp(mod_id{m},'MRI_ESM2_1') || strcmp(mod_id{m},'MPIOM_HAMOCC') || strcmp(mod_id{m},'ROMS_SouthernOcean_ETHZ')
        tmp=ncread([prefix 'MLD_variab_thr/mld_' exper '_' models{m} '.nc'],'mld');
        eval(['mld_yy_vthr_',mod_id{m},'=-squeeze(nanmean(tmp,3));'])
        
        else
       tmp=ncread([prefix 'MLD_variab_thr/mld_' exper '_' models{m} '.nc'],'mld');
        eval(['mld_yy_vthr_',mod_id{m},'=-tmp;'])
        
    end
      
        save([prefix models{m} '/mld_yy_var_thr_' models{m} '_clim_' ver{m} '_94_07.mat'],['mld_yy_vthr_',mod_id{m}])
    end
end


%%
%#######################################################
% save vertical integrals models
%#######################################################

year='2007';
for m=1:numel(models)
load([prefix models{m} '/Cant_years/cant_' mod_id{m} '_' year '.mat']);
end

for m=1:numel(models)

    models{m}
    load([prefix models{m} '/lat_' models{m} '_1_gr_1980-2018_' ver{m} '.mat']);
    load([prefix models{m} '/lon_' models{m} '_1_gr_1980-2018_' ver{m} '.mat']);
    load([prefix models{m} '/volume_' models{m} '_1_gr_1980-2018_' ver{m} '.mat']);
    load([prefix models{m} '/mask_vol_' models{m} '_1_gr_1980-2018_' ver{m} '.mat']);
    eval(['mask_vol=mask_vol_',mod_id{m},';'])
    
    load([prefix models{m} '/depth_' models{m} '_1_gr_1980-2018_' ver{m} '.mat']);
    %put mask_vol to NaN below 3000 m depth
    eval(['depth_mask=repmat(depth_',mod_id{m} ',[1 360 180]);'])
    depth_mask=permute(depth_mask,[2 3 1]);
    mask_vol(mask_vol==0)=NaN;
    mask_vol(depth_mask>3000)=NaN;
     
    load([prefix models{m} '/dz_' models{m} '_1_gr_1980-2018_' ver{m} '.mat']);
    eval(['dz_rep=repmat(dz_',mod_id{m},',[1 360 180]);'])
    dz_rep=permute(dz_rep,[2 3 1]);
    dz_sum=squeeze(nansum(dz_rep.*mask_vol,3));
    
    if strcmp(mod_id{m},'MPIOM_HAMOCC') || strcmp(mod_id{m},'OCIM_v2021') || strcmp(mod_id{m},'CNRM_ESM2_1')
        
        eval(['cant_',mod_id{m},'_vint_' year '_CminusB=squeeze(nansum(cant_',mod_id{m},'_CminusB.*dz_rep.*mask_vol,3));'])
        eval(['cant_',mod_id{m},'_vint_' year '_CminusB(squeeze(mask_vol(:,:,1))==0)=NaN;']);
        save([prefix models{m} '/cant_vert_int_' models{m} '_1_gr_' year '_' ver{m} '_CminusB.mat'],['cant_',mod_id{m},'_vint_' year '_CminusB'])
        
    else
        
        eval(['cant_',mod_id{m},'_vint_' year '_AminusD=squeeze(nansum(cant_',mod_id{m},'_AminusD.*dz_rep.*mask_vol,3));'])
        eval(['cant_',mod_id{m},'_vint_' year '_AminusD=cant_',mod_id{m},'_vint_' year '_AminusD.*squeeze(mask_vol(:,:,1));']);        
        save([prefix models{m} '/cant_vert_int_' models{m} '_1_gr_' year '_' ver{m} '_AminusD.mat'],['cant_',mod_id{m},'_vint_' year '_AminusD'])
        
        eval(['cant_',mod_id{m},'_vint_' year '_CminusB=squeeze(nansum(cant_',mod_id{m},'_CminusB.*dz_rep.*mask_vol,3));'])
        eval(['cant_',mod_id{m},'_vint_' year '_CminusB(squeeze(mask_vol(:,:,1))==0)=NaN;']);
        save([prefix models{m} '/cant_vert_int_' models{m} '_1_gr_' year '_' ver{m} '_CminusB.mat'],['cant_',mod_id{m},'_vint_' year '_CminusB'])
    end
end

%%
%#######################################################
% save vertical integral eMLRC* 
%#######################################################

% Gruber data set for 1994-2007 (use vertically-integrated)
cant_grub=ncread([prefix 'Gruber2019_Cant/inv_dcant_emlr_cstar_gruber_94-07_vs1.nc'],'DCANT_INV01'); %umol/kg
cant_grub_inv=permute(cant_grub,[2 1]);
save([prefix 'Gruber2019_Cant_inventory.mat'],'cant_grub_inv')


%%
%#######################################################
% save zonal integrals and averages Observations
%#######################################################
% note: permute dimensions to be depth, lat, lon

gamma_grub=ncread([prefix 'Gruber2019_Cant/dcant_emlr_cstar_gruber_94-07_vs1.nc'],'GAMMA_DENS'); %kg/m3
gamma_grub=permute(gamma_grub,[3 2 1])+1000;
cant_grub=ncread([prefix 'Gruber2019_Cant/dcant_emlr_cstar_gruber_94-07_vs1.nc'],'DCANT_01'); % in umol/kg
cant_grub=permute(cant_grub,[3 2 1]);
lon_grub=ncread([prefix 'Gruber2019_Cant/dcant_emlr_cstar_gruber_94-07_vs1.nc'],'LONGITUDE');
lat_grub=ncread([prefix 'Gruber2019_Cant/dcant_emlr_cstar_gruber_94-07_vs1.nc'],'LATITUDE');
dep_grub=ncread([prefix 'Gruber2019_Cant/dcant_emlr_cstar_gruber_94-07_vs1.nc'],'DEPTH'); % same as Sabine data set

[latgrub,depgrub]=meshgrid(lat_grub,dep_grub);

load([prefix 'Gruber2019_Cant_inventory.mat'])

% sigma from WOA 
load([prefix 'WOA/sigma0_WOA_clim_all_lev.mat']);
sigma0_woa_clim=permute(sigma0_woa_clim,[3 2 1]);

load([prefix '1x1_dx_dy.mat'])
mask_vol=ones(size(sigma0_woa_clim));
mask_vol(isnan(sigma0_woa_clim))=NaN;

dx_rep=repmat(dx,[1 1 size(mask_vol,1)]);
dx_rep=permute(dx_rep,[3 1 2]);
lonsum_glob=squeeze(nansum(dx_rep.*mask_vol,3));

sigma_woa_clim_glob=squeeze(nansum(sigma0_woa_clim.*dx_rep.*mask_vol,3))./lonsum_glob;
depth_woa=double(depth);

[latwoa,depwoa]=meshgrid(lat_grub,depth_woa);

sigma_woa_regr=griddata(latwoa,depwoa,sigma_woa_clim_glob,latgrub,depgrub); % interpolate sigma on gruber grid

save([prefix 'sigma0_woa_zonmean_clim.mat'],'sigma_woa_regr','depth_woa')

% zonal integrals and averages

mask_vol=ones(size(cant_grub));
mask_vol(isnan(cant_grub))=NaN;

dx_rep=repmat(dx,[1 1 size(mask_vol,1)]);
dx_rep=permute(dx_rep,[3 1 2]);
lonsum_glob=squeeze(nansum(dx_rep.*mask_vol,3));

cant_gruber_zonave=squeeze(nansum(cant_grub.*dx_rep.*mask_vol.*sigma_woa_regr,3))./lonsum_glob/1e6; % units: umol/kg * kg/m3 / 1e6--> mol/m3 
cant_gruber_inv_zonave=squeeze(nansum(cant_grub_inv.*dx.*squeeze(mask_vol(1,:,:)),2))./squeeze(lonsum_glob(1,:))'; % units: mol/m2 
cant_gruber_zonint=squeeze(nansum(cant_grub.*dx_rep.*mask_vol.*sigma_woa_regr,3))/1e6; % units: umol/kg * m * kg/m3 / 1e6--> mol/m2 
cant_gruber_inv_zonint=squeeze(nansum(cant_grub_inv.*dx.*squeeze(mask_vol(1,:,:)),2)); % units: mol/m1

save([prefix 'cant_zonint_gruber_glob.mat'],'cant_gruber_zonave','cant_gruber_inv_zonave',...
    'cant_gruber_zonint','cant_gruber_inv_zonint','dep_grub','lat_grub')

%%
%#######################################################
% save zonal integrals and averages models
%#######################################################

for m=1:numel(models)
    models{m}
    
    % cant interior    
    load([prefix models{m} '/depth_' models{m} '_1_gr_1980-2018_' ver{m} '.mat']);
    load([prefix models{m} '/mask_vol_' models{m} '_1_gr_1980-2018_' ver{m} '.mat']);
    load([prefix models{m} '/dz_' models{m} '_1_gr_1980-2018_' ver{m} '.mat']);
 
    eval(['depth=depth_' mod_id{m} ';'])
    eval(['mask_vol=mask_vol_' mod_id{m} ';'])      
    eval(['dz=dz_' mod_id{m} ';'])
    
    %put mask_vol to NaN below 3000 m depth
    load([prefix models{m} '/depth_' models{m} '_1_gr_1980-2018_' ver{m} '.mat']);    
    eval(['depth_mask=repmat(depth_',mod_id{m} ',[1 360 180]);'])
    depth_mask=permute(depth_mask,[2 3 1]);
    mask_vol(mask_vol==0)=NaN;
    mask_vol_all=mask_vol;
    mask_vol(depth_mask>3000)=NaN;
  
    load([prefix '1x1_dx_dy.mat']);
    dx=dx';
    dx_rep=repmat(dx,[1 1 numel(depth)]);
    dz_rep=repmat(dz,[1 size(dx)]);
    dz_rep=permute(dz_rep,[2 3 1]);
    lonsum_glob=squeeze(nansum(dx_rep.*mask_vol_all));
    lonsum_1d_glob=squeeze(lonsum_glob(:,1))';
   
    if strcmp(mod_id{m},'MPIOM_HAMOCC') || strcmp(mod_id{m},'OCIM_v2021') || strcmp(mod_id{m},'CNRM_ESM2_1')
        load([prefix models{m} '/Cant_years/cant_' mod_id{m} '_1994.mat']);
        eval(['cant_1994=cant_',mod_id{m},'_CminusB;']);
        load([prefix models{m} '/Cant_years/cant_' mod_id{m} '_2007.mat']);
        eval(['cant_2007=cant_',mod_id{m},'_CminusB;']);
        
    else
        load([prefix models{m} '/Cant_years/cant_' mod_id{m} '_1994.mat']);
        eval(['cant_1994=cant_',mod_id{m},'_AminusD;']);
        load([prefix models{m} '/Cant_years/cant_' mod_id{m} '_2007.mat']);
        eval(['cant_2007=cant_',mod_id{m},'_AminusD;']);
    end
    
    delta_cant_3D=cant_2007-cant_1994;
    delta_cant_3D=double(delta_cant_3D);
    
    eval(['cant_',mod_id{m},'_zonave=squeeze(nansum(delta_cant_3D.*dx_rep.*mask_vol))./lonsum_glob;']) % mol/m3
    eval(['cant_',mod_id{m},'_vinv_zonave=squeeze(nansum(nansum(delta_cant_3D.*dx_rep.*dz_rep.*mask_vol,1),3))./lonsum_1d_glob;']) % mol/m2
    eval(['cant_',mod_id{m},'_zonint=squeeze(nansum(delta_cant_3D.*dx_rep.*mask_vol));']) % mol/m2
    eval(['cant_',mod_id{m},'_vinv_zonint=squeeze(nansum(nansum(delta_cant_3D.*dx_rep.*dz_rep.*mask_vol,1),3));']) % mol/m
    eval(['cant_',mod_id{m},'_vinvall_zonint=squeeze(nansum(nansum(delta_cant_3D.*dx_rep.*dz_rep.*mask_vol_all,1),3));']) % mol/m
    
    % regrid on eMLRC* grid
    [latmod,depmod]=meshgrid(lat_grub,depth);
    latmod=latmod';
    depmod=depmod';

    eval(['cant_',mod_id{m},'_zonave_regrid=griddata(latmod,depmod,cant_',mod_id{m},'_zonave,latgrub,depgrub);'])
     eval(['cant_',mod_id{m},'_zonint_regrid=griddata(latmod,depmod,cant_',mod_id{m},'_zonint,latgrub,depgrub);'])

    save([prefix models{m} '/zonmeans/delta_cant_sections_' models{m} '_1994_2007_' ver{m} '.mat'],...
        ['cant_',mod_id{m},'_zonave_regrid'],['cant_',mod_id{m},'_vinv_zonave'],['cant_',mod_id{m},'_zonint_regrid'],...
        ['cant_',mod_id{m},'_vinv_zonint'],['cant_',mod_id{m},'_vinvall_zonint'])  

    % sigma
   if ~strcmp(mod_id{m},'OCIM_v2021')
       if strcmp(mod_id{m},'MPIOM_HAMOCC')==1 || strcmp(mod_id{m},'CNRM_ESM2_1')
           exper='B'; 
       else
           exper='A';
       end

    load([prefix models{m} '/sigma0_' mod_id{m} '_' exper '_1994_2007_all_lev.mat'])
    sigma0_CNRM_ESM2_1_B=sigma0_CNRM_ESM2_1_A;
    eval(['sig0_',mod_id{m},'_zon=squeeze(nansum(sigma0_',mod_id{m},'_' exper,'.*dx_rep.*mask_vol_all))./lonsum_glob;'])
    
    eval(['sig0_',mod_id{m},'_zon=double(sig0_',mod_id{m},'_zon);'])
    eval(['sig0_',mod_id{m},'_zon_regrid=griddata(latmod,depmod,sig0_',mod_id{m},'_zon,latgrub,depgrub);'])
    
    
    save([prefix models{m} '/zonmeans/sigma0_zonave_' models{m} '_1994_2007_' ver{m} '.mat'],...
        ['sig0_',mod_id{m},'_zon'],['sig0_',mod_id{m},'_zon_regrid'])
    
   end
   
   % cantflux
   load([prefix models{m} '/cantflux_' models{m} '_clim_' ver{m} '_94_07.mat']); %,['cantfl_',mod_id{m}])  
   
   eval(['cantfl_',mod_id{m},'_zonave=squeeze(nansum(cantfl_',mod_id{m},'.*dx.*mask_vol(:,:,1)))./lonsum_1d_glob;'])
   eval(['cantfl_',mod_id{m},'_zonint=squeeze(nansum(cantfl_',mod_id{m},'.*dx.*mask_vol(:,:,1)));'])
   
   save([prefix models{m} '/zonmeans/cantflux_sections_' models{m} '_1994_2007_' ver{m} '_glob.mat'],...
       ['cantfl_',mod_id{m},'_zonave'],['cantfl_',mod_id{m},'_zonint'])
    
end
%%
%#######################################################
% save volume integrals and area averages Models
%#######################################################

% load RECCAP regional mask
tmp=ncread([prefix 'RECCAP2_region_masks_all_v20210412.nc'],'southern');
sou_mask=double(tmp);
spss_stss_mask=double(tmp);

sou_mask(sou_mask==0)=NaN;
spss_stss_mask(spss_stss_mask==0)=NaN;

spss_stss_mask(spss_stss_mask==1)=1;
spss_stss_mask(spss_stss_mask==2)=1;
spss_stss_mask(spss_stss_mask==3)=NaN;

sou_mask(~isnan(sou_mask))=1;

% models
for m=1:numel(models)
    
    mod_id{m}
    
    load([prefix models{m} '/area_' models{m} '_1_gr_1980-2018_' ver{m} '.mat']);
    load([prefix models{m} '/mask_sfc_' models{m} '_1_gr_1980-2018_' ver{m} '.mat']);
    load([prefix models{m} '/volume_' models{m} '_1_gr_1980-2018_' ver{m} '.mat']);
    load([prefix models{m} '/mask_vol_' models{m} '_1_gr_1980-2018_' ver{m} '.mat']);
    
    
    eval(['mask_sfc=mask_sfc_',mod_id{m},';'])
    
    if strcmp(mod_id{m},'planktom12')
        eval(['area=area_',mod_id{m},'*1e-6;']) % area is 6 order of magnitude higher than it should be
    end
    if strcmp(mod_id{m},'CESM_ETHZ')
        eval(['area=-area_',mod_id{m},';'])
    end
    if ~strcmp(mod_id{m},'planktom12') & ~strcmp(mod_id{m},'CESM_ETHZ')
        eval(['area=area_',mod_id{m},';'])
    end
    
    area_SO=squeeze(nansum(nansum(area.*mask_sfc.*sou_mask)));
    area_spss_stss=squeeze(nansum(nansum(area.*mask_sfc.*spss_stss_mask)));
    
    % volume SO
    eval(['mask_vol=mask_vol_',mod_id{m},';'])
    eval(['volume=volume_',mod_id{m},';'])
    sou_rep=repmat(sou_mask,[1 1 size(volume,3)]);
    
    %put mask_vol to NaN below 3000 m depth
    load([prefix models{m} '/depth_' models{m} '_1_gr_1980-2018_' ver{m} '.mat']);    
    eval(['depth_mask=repmat(depth_',mod_id{m} ',[1 360 180]);'])
    depth_mask=permute(depth_mask,[2 3 1]);
    mask_vol(mask_vol==0)=NaN;
     
    mask_vol(depth_mask>3000)=NaN;

    
    eval(['volume_',mod_id{m},'_SO=squeeze(nansum(nansum(nansum(volume.*mask_vol.*sou_rep))));'])
    
    % cantflux (area integral)
    load([prefix models{m} '/cantflux_' models{m} '_clim_' ver{m} '_94_07.mat']);  
    eval(['cantfl_int_',mod_id{m},'_SO=squeeze(nansum(nansum(cantfl_',mod_id{m},'.*area.*mask_sfc.*sou_mask)));'])
  
    % SSS (area average)
    load([prefix models{m} '/sos_' models{m} '_clim_' ver{m} '_94_07.mat']);  
    eval(['sos_',mod_id{m},'_SO=squeeze(nansum(nansum(sos_',mod_id{m},'.*area.*mask_sfc.*spss_stss_mask)))/area_spss_stss;'])
    
    % cant (volume integral)
    
    if strcmp(mod_id{m},'MPIOM_HAMOCC') || strcmp(mod_id{m},'OCIM_v2021') || strcmp(mod_id{m},'CNRM_ESM2_1')
        
        load([prefix models{m} '/Cant_years/cant_' mod_id{m} '_1994.mat']);
        eval(['var_1994=cant_',mod_id{m},'_CminusB;'])
        load([prefix models{m} '/Cant_years/cant_' mod_id{m} '_2007.mat']);
        eval(['var_2007=cant_',mod_id{m},'_CminusB;'])
        delta_cant=var_2007-var_1994;
        
    else
        load([prefix models{m} '/Cant_years/cant_' mod_id{m} '_1994.mat']);
        eval(['var_1994=cant_',mod_id{m},'_AminusD;'])
        load([prefix models{m} '/Cant_years/cant_' mod_id{m} '_2007.mat']);
        eval(['var_2007=cant_',mod_id{m},'_AminusD;'])
        delta_cant=var_2007-var_1994;
    end
    
    % compute integral and convert from mol/m3 to PgC/year: mol/m3 * g/mol
    % * m3 /1e15 /13 (years) (divide by 13 only the delta Cant, not the
    % Cant in 2007)
    eval(['delta_cant_',mod_id{m},'_int=squeeze(nansum(nansum(nansum(delta_cant*12.01.*volume.*mask_vol.*sou_rep))))/1e15/13;']) 
    eval(['cant_2007_',mod_id{m},'_int=squeeze(nansum(nansum(nansum(var_2007*12.01.*volume.*mask_vol.*sou_rep))))/1e15;']) 
    eval(['cant_1994_',mod_id{m},'_int=squeeze(nansum(nansum(nansum(var_1994*12.01.*volume.*mask_vol.*sou_rep))))/1e15;'])

    save([prefix models{m} '/fig12_integrals_sou_' models{m} '_1_gr_' ver{m} '.mat'],['delta_cant_',mod_id{m},'_int'],...
      ['sos_',mod_id{m},'_SO'],['cant_1994_',mod_id{m},'_int'],['cant_2007_',mod_id{m},'_int'],['cantfl_int_',mod_id{m},'_SO'])
  
   
end

%%
%#######################################################
% save volume integrals and area averages Observations
%#######################################################
load([prefix '1x1_dx_dy.mat']);

tmp=ncread([prefix 'RECCAP2_region_masks_all_v20210412.nc'],'southern');
sou_mask=double(tmp);
spss_stss_mask=double(tmp);

sou_mask(sou_mask==0)=NaN;
spss_stss_mask(spss_stss_mask==0)=NaN;

spss_stss_mask(spss_stss_mask==1)=1;
spss_stss_mask(spss_stss_mask==2)=1;
spss_stss_mask(spss_stss_mask==3)=NaN;

sou_mask(~isnan(sou_mask))=1;
sou_mask=sou_mask';
spss_stss_mask=spss_stss_mask';

area=dx.*dy.*sou_mask;

cant_sabine=ncread([prefixo 'Sabine2005_Cant/glodap_Anth_CO2.nc'],'anth_co2'); %umol/kg
cant_sabine=permute(cant_sabine,[3 2 1]);
depthedges_sabine=ncread([prefixo 'Sabine2005_Cant/glodap_Anth_CO2.nc'],'DEPTHedges');
depth_sabine=ncread([prefixo 'Sabine2005_Cant/glodap_Anth_CO2.nc'],'DEPTH');
lat_sabine=ncread([prefixo 'Sabine2005_Cant/glodap_Anth_CO2.nc'],'LATITUDE');
lon_sabine=ncread([prefixo 'Sabine2005_Cant/glodap_Anth_CO2.nc'],'LONGITUDE360_720');

cant_sab_resh=NaN([33 180 360]);
cant_sab_resh(:,:,1:180)=cant_sabine(:,:,181:360); % last long value is = the first --> disregard lon 361
cant_sab_resh(:,:,181:360)=cant_sabine(:,:,1:180); 
cant_sab_resh(29:end,:,:)=NaN; % only upper 3000 m
cant_sab_resh(cant_sab_resh<0)=0;

% convert to mol/m3 (use density from Gruber data set)
gamma_grub=ncread([prefixo 'Gruber2019_Cant/dcant_emlr_cstar_gruber_94-07_vs1.nc'],'GAMMA_DENS'); %kg/m3
gamma_grub=permute(gamma_grub,[3 2 1]);

% compute dz for Sabine et al data set
new_bord=[depth_sabine(1:end-1) depth_sabine(2:end)];
new_bord=mean(new_bord,2);
new_bord=[0 ; new_bord ; depth_sabine(end)+(depth_sabine(end)- new_bord(end))];
dz_sabine=diff(new_bord);

dz_sab_3d=repmat(dz_sabine,[1 180 360]);

area_rep=repmat(area,[1 1 size(dz_sab_3d,1)]);
area_rep=permute(area_rep,[3 1 2]);
volume=area_rep.*dz_sab_3d;


% Delta Cant 1994-2007 (Gruber 2019) in umol/kg
gamma_grub=ncread([prefixo 'Gruber2019_Cant/dcant_emlr_cstar_gruber_94-07_vs1.nc'],'GAMMA_DENS'); %kg/m3
gamma_grub=permute(gamma_grub,[3 2 1]);
cant_grub=ncread([prefixo 'Gruber2019_Cant/dcant_emlr_cstar_gruber_94-07_vs1.nc'],'DCANT_01'); %umol/kg
cant_grub=permute(cant_grub,[3 2 1]);
cant_grub(29:end,:,:)=NaN; % only upper 3000 m 
cant_2007=cant_sab_resh+cant_grub; %umol/kg

volume_cant=volume;
volume_cant(isnan(cant_2007))=NaN;

%compute integral and convert from umol/kg to PgC and PgC/yr (divide by 13
%years delta Cant) --> umol/kg * kg/m3 * m3 * g/mol /1e15 / 1e6
delta_cant_obs_int=nansum(nansum(nansum(nansum(cant_grub.*(gamma_grub+1000).*volume))))/1e15/1e6*12.01/13
cant_2007_obs_int=nansum(nansum(nansum(nansum(cant_2007.*(gamma_grub+1000).*volume))))/1e15/1e6*12.01;
cant_1994_obs_int=nansum(nansum(nansum(nansum(cant_sab_resh.*(gamma_grub+1000).*volume))))/1e15/1e6*12.01

% SSS SOMFFN
prefix=[dir 'Elements/RECCAP/Observations/SOMFFN_v20211121/'];
tmp=ncread([prefix 'sss_MPI_SOMFFN_1982-2019_v20211121.nc'],'sss');
tmp=double(tmp);
tmp=reshape(tmp,[12 38 size(tmp,2) size(tmp,3)]);
tmp=squeeze(nanmean(tmp,1));

sos_somffn=squeeze(nanmean(tmp(13:26,:,:))); % 1994-2007
sos_somffn(sos_somffn<0)=NaN;

tmp=ncread([prefix 'area_MPI_SOMFFN_1982-2019_v20211121.nc'],'area');
area_somffn=double(tmp);
area_somffn=area_somffn';
area_somffn(isnan(sos_somffn))=NaN;
areasum_somffn=squeeze(nansum(nansum(area_somffn.*spss_stss_mask)));

sos_somffn_mean=nansum(nansum(sos_somffn.*area_somffn.*spss_stss_mask))./areasum_somffn;

% save
save([prefixo 'SO_area_ave_observations_MLD_Cant.mat'],'cant_2007_obs_int','cant_1994_obs_int',...
    'delta_cant_obs_int','sos_somffn_mean')
