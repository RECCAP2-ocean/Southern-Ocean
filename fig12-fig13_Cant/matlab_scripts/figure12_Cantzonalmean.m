%% 
prefix='../ready_to_plot_netcdfs/';
filename='figure12_Cantzonalmean.nc';

cant_gruber=ncread([prefix filename],'cant_gruber');
cant_mmm_all=ncread([prefix filename],'cant_mmm_all');
cant_OCIM_v2021=ncread([prefix filename],'cant_OCIM_v2021');
sigma_woa_regr=ncread([prefix filename],'sigma_woa_regr');
sig_mmm=ncread([prefix filename],'sig_mmm');
cant_gruber_inv=ncread([prefix filename],'cant_gruber_inv');
cant_OCIM_v2021_vinv=ncread([prefix filename],'cant_OCIM_v2021_vinv');
cant_vinv_mmm_high=ncread([prefix filename],'cant_vinv_mmm_high');
cant_vinv_mmm_low=ncread([prefix filename],'cant_vinv_mmm_low');
cantfl_OCIM_v2021=ncread([prefix filename],'cantfl_OCIM_v2021');
cantfl_mmm_high=ncread([prefix filename],'cantfl_mmm_high');
cantfl_mmm_low=ncread([prefix filename],'cantfl_mmm_low');
min_ice_mean=ncread([prefix filename],'min_ice_mean');
max_ice_mean=ncread([prefix filename],'max_ice_mean');
max_spss_mean=ncread([prefix filename],'max_spss_mean');
max_stss_mean=ncread([prefix filename],'max_stss_mean');
lat=ncread([prefix filename],'lat');
depth=ncread([prefix filename],'depth');

%%         

subpl_pos_up={[0.11 0.75 0.28 0.2],[0.4 0.75 0.28 0.2],[0.69 0.75 0.28 0.2]}                  
subpl_pos_bot={[0.11 0.53 0.28 0.21],[0.4 0.53 0.28 0.21],[0.69 0.53 0.28 0.21]}      
load('Colormaps_Spectral_RdBu_r.mat')
%%
% Conversions:
% For Delta Cant, divide  by 13 years (2007-1994) to obtain accumulation rates (mol/m2/yr for Cant sections, mol/m/yr for Cant inventories).
% For Cant air-sea fluxes, multiply seconds_per_year to convert to mol/m/yr.

figure(1)
clf

%%%% Cant inventories and Cant air-sea fluxes %%%%

subplot(3,2,1,'Position',subpl_pos_up{1})

var_cantv=cant_gruber_inv/13*1e-7;
plot(lat,var_cantv,'Color','k','LineWidth',1.2)

hold on
ll=line([lat(min_ice_mean) lat(min_ice_mean)],[3.05 3.5])
set(ll,'Color',[0 0.8 0.933],'LineWidth',0.7)
hold on
ll=line([lat(max_ice_mean) lat(max_ice_mean)],[03.05 3.5])
set(ll,'Color',[0 0.8 0.933],'LineWidth',0.7)
hold on
ll1=line([lat(max_spss_mean) lat(max_spss_mean)],[03.05 3.5])
set(ll1,'Color',[0 0.8 0.933],'LineWidth',0.7)
hold on
ll2=line([lat(max_stss_mean) lat(max_stss_mean)],[03.05 3.5])
set(ll2,'Color',[0 0.8 0.933],'LineWidth',0.7)

t=text(-69,2.8,'ICE')
set(t,'FontSize',9.0,'FontWeight','bold','Color',[0 0.8 0.933])
t=text(-60.5,2.8,'SPSS')
set(t,'FontSize',9.0,'FontWeight','bold','Color',[0 0.8 0.933])
t=text(-50,2.8,'STSS')
set(t,'FontSize',9.0,'FontWeight','bold','Color',[0 0.8 0.933])

ylabel({'10^7','mol m^-^1 yr^-^1'},'Color','k')

set(gca,'FontSize',11.0,'XLim',[-80 -20],'YLim',[0 3.5],'XTickLabel',[],'YTick',[0:1.5:3.5],'Xgrid','on','YGrid','on')

title('a) eMLR(C*)','FontSize',12.0)

t=text(-52,0.6,'\DeltaC_{ant}')
set(t,'FontSize',11.0,'FontWeight','bold')


%%%%%% 

subplot(3,2,2,'Position',subpl_pos_up{2})

var_plot=cantfl_OCIM_v2021*86400*365*1e-7; 
var_cantv=cant_OCIM_v2021_vinv/13*1e-7; 

plot(lat,var_cantv,'Color','k','LineWidth',1.2)
hold on 
plot(lat,var_plot,'Color',[0.5 0.5 0.5],'LineWidth',1.2)

hold on
ll=line([lat(min_ice_mean) lat(min_ice_mean)],[3.05 3.5])
set(ll,'Color',[0 0.8 0.933],'LineWidth',0.7)
hold on
ll=line([lat(max_ice_mean) lat(max_ice_mean)],[03.05 3.5])
set(ll,'Color',[0 0.8 0.933],'LineWidth',0.7)
hold on
ll1=line([lat(max_spss_mean) lat(max_spss_mean)],[03.05 3.5])
set(ll1,'Color',[0 0.8 0.933],'LineWidth',0.7)
hold on
ll2=line([lat(max_stss_mean) lat(max_stss_mean)],[03.05 3.5])
set(ll2,'Color',[0 0.8 0.933],'LineWidth',0.7)


set(gca,'FontSize',11.0,'XLim',[-80 -20],'YLim',[0 3.5],'XTickLabel',[],'YTickLabel',[],'YTick',[0:1.5:3.5],'Xgrid','on','YGrid','on')
title('b) OCIM v2021 ','FontSize',12.0)

t=text(-52,0.6,'\DeltaC_{ant}')
set(t,'FontSize',11.0,'FontWeight','bold')
t=text(-72,2,{'C_{ant}','flux'})
set(t,'FontSize',11.0,'FontWeight','bold','Color',[0.5 0.5 0.5],'HorizontalAlignment','center')

%
subplot(3,2,3,'Position',subpl_pos_up{3})
var_plot_high=cantfl_mmm_high*86400*365*1e-7;
var_plot_low=cantfl_mmm_low*86400*365*1e-7;
var_cantv_high=cant_vinv_mmm_high/13*1e-7;
var_cantv_low=cant_vinv_mmm_low/13*1e-7;

plot(lat,var_cantv_high,'Color','k','LineWidth',1.2)
hold on 
plot(lat,var_cantv_low,'Color','k','LineWidth',1.2,'LineStyle','--')
hold on 
plot(lat,var_plot_high,'Color',[0.5 0.5 0.5],'LineWidth',1.2)
hold on 
plot(lat,var_plot_low,'Color',[0.5 0.5 0.5],'LineWidth',1.2,'LineStyle','--')

hold on
ll=line([lat(min_ice_mean) lat(min_ice_mean)],[3.05 3.5])
set(ll,'Color',[0 0.8 0.933],'LineWidth',0.7)
hold on
ll=line([lat(max_ice_mean) lat(max_ice_mean)],[03.05 3.5])
set(ll,'Color',[0 0.8 0.933],'LineWidth',0.7)
hold on
ll1=line([lat(max_spss_mean) lat(max_spss_mean)],[03.05 3.5])
set(ll1,'Color',[0 0.8 0.933],'LineWidth',0.7)
hold on
ll2=line([lat(max_stss_mean) lat(max_stss_mean)],[03.05 3.5])
set(ll2,'Color',[0 0.8 0.933],'LineWidth',0.7)

set(gca,'FontSize',11.0,'XLim',[-80 -20],'YLim',[0 3.5],'XTickLabel',[],'YTickLabel',[],'YTick',[0:1.5:3.5],'Xgrid','on','YGrid','on')
title('c) GOBMs','FontSize',12.0)

t=text(-52,0.6,'\DeltaC_{ant}')
set(t,'FontSize',11.0,'FontWeight','bold')
t=text(-72,2,{'C_{ant}','flux'})
set(t,'FontSize',11.0,'FontWeight','bold','Color',[0.5 0.5 0.5],'HorizontalAlignment','center')

l=legend('GOBMs high','GOBMs low')
set(l,'FontSize',9.0)

%%%%% interior
subplot(3,2,4,'Position',subpl_pos_bot{1})
cant_grub_plot=cant_gruber/13*1e-4; 
pcolor(lat,depth/1000,cant_grub_plot)
shading flat
caxis([0 3])
colormap(flipud(Spectral))
hold on
contour(lat,depth/1000,sigma_woa_regr,[1026:0.3:1028],'LineWidth',0.5,'Color',[0 0 0])
hold on
contour(lat,depth/1000,sigma_woa_regr,[1026.9 1027.5],'LineWidth',1.5,'Color',[0 0 0])

set(gca,'Layer','top')
set(gca,'FontSize',11.0,'YDir','reverse','YLim',[0 3],'YTick',[1 2 3],'XLim',[-80 -20],'XTick',[-60 -40],'XTicklabel',{'60^oS','40^oS'})

ylabel('depth (km)')

t=text(-78,2.6,'d)')
set(t,'FontSize',13.0,'FontWeight','bold','BackgroundColor',[1 1 1])
xlabel('latitude')

% ocim

subplot(3,2,5,'Position',subpl_pos_bot{2})
var_plot=cant_OCIM_v2021/13*1e-4;

pcolor(lat,depth/1000,var_plot)
hold on
contour(lat,depth/1000,sigma_woa_regr,[1026:0.3:1028],'LineWidth',0.5,'Color',[0 0 0])
hold on
contour(lat,depth/1000,sigma_woa_regr,[1026.9 1027.5],'LineWidth',1.5,'Color',[0 0 0])
shading flat
caxis([0 3])

colormap(flipud(Spectral))
set(gca,'Layer','top')
set(gca,'FontSize',11.0,'YDir','reverse','YLim',[0 3],'YTick',[1 2 3],'YTicklabel',[],'XLim',[-80 -20],'XTick',[-60 -40],'XTicklabel',{'60^oS','40^oS'})

t=text(-78,2.6,'e)')
set(t,'FontSize',13.0,'FontWeight','bold','BackgroundColor',[1 1 1])
xlabel('latitude')

% mmm
subplot(3,2,6,'Position',subpl_pos_bot{3})
var_plot=cant_mmm_all/13;
pcolor(lat,depth/1000,var_plot*1e-4)
hold on
contour(lat,depth/1000,sig_mmm,[1026:0.3:1028],'LineWidth',0.5,'Color',[0 0 0])
hold on
contour(lat,depth/1000,sig_mmm,[1026.9 1027.5],'LineWidth',1.5,'Color',[0 0 0])
shading flat
caxis([0 3])
colormap(flipud(Spectral))
set(gca,'Layer','top')
set(gca,'FontSize',11.0,'YDir','reverse','YLim',[0 3],'YTick',[1 2 3],'YTicklabel',[],'XLim',[-80 -20],'XTick',[-60 -40],'XTicklabel',{'60^oS','40^oS'})

xlabel('latitude')

cb=colorbar;
set(cb,'Location','SouthOutside','Position',[0.35 0.37 0.35 0.03],'FontSize',12.0)
t=text(-53,5.1,{'\DeltaC_{ant} rate','(10^4 mol m^-^2 yr^-^1)'})
set(t,'FontSize',11.0,'FontWeight','bold','HorizontalAlignment','center')
t=text(-78,2.6,'f)')
set(t,'FontSize',13.0,'FontWeight','bold','BackgroundColor',[1 1 1])


%%%% Anomalies %%%%

cant_ocim_m_emlrc=(cant_OCIM_v2021-cant_gruber)/13;
cant_gobm_all_m_emlrc=(cant_mmm_all-cant_gruber)/13;
cant_gobm_all_m_ocim=(cant_mmm_all-cant_OCIM_v2021)/13;


figure(2)
clf

subplot(3,3,1,'Position',subpl_pos_bot{1})
pcolor(lat,depth/1000,cant_ocim_m_emlrc*1e-4)
hold on
contour(lat,depth/1000,sigma_woa_regr,[1026:0.3:1028],'LineWidth',0.5,'Color',[0 0 0])
hold on
contour(lat,depth/1000,sigma_woa_regr,[1026.9 1027.5],'LineWidth',1.5,'Color',[0 0 0])
shading flat
caxis([-0.3 0.3])
colormap(RdBu_r)
shading flat
set(gca,'Layer','top')
set(gca,'FontSize',11.0,'YDir','reverse','YLim',[0 3],'YTick',[1 2 3],'XLim',[-80 -20],'XTick',[-60 -40],'XTicklabel',{'60^oS','40^oS'})
ylabel('depth (km)')
xlabel('latitude')
title({'g) OCIM-v2021 -','eMLR(C*)'})

subplot(3,3,2,'Position',subpl_pos_bot{2})

pcolor(lat,depth/1000,cant_gobm_all_m_emlrc*1e-4)
hold on
contour(lat,depth/1000,sigma_woa_regr,[1026:0.3:1028],'LineWidth',0.5,'Color',[0 0 0])
hold on
contour(lat,depth/1000,sigma_woa_regr,[1026.9 1027.5],'LineWidth',1.5,'Color',[0 0 0])
shading flat
caxis([-0.3 0.3])
colormap(RdBu_r)
shading flat
set(gca,'Layer','top')
set(gca,'FontSize',11.0,'YDir','reverse','YLim',[0 3],'YTick',[1 2 3],'XLim',[-80 -20],'XTick',[-60 -40],'XTicklabel',{'60^oS','40^oS'},'YTicklabel',{})
title({'h) GOBMs -',' eMLR(C*)'})
xlabel('latitude')

%
subplot(3,3,3,'Position',subpl_pos_bot{3})
pcolor(lat,depth/1000,cant_gobm_all_m_ocim*1e-4)
hold on
contour(lat,depth/1000,sigma_woa_regr,[1026:0.3:1028],'LineWidth',0.5,'Color',[0 0 0])
hold on
contour(lat,depth/1000,sigma_woa_regr,[1026.9 1027.5],'LineWidth',1.5,'Color',[0 0 0])
shading flat
caxis([-0.3 0.3])
colormap(RdBu_r)
shading flat
set(gca,'Layer','top')
set(gca,'FontSize',11.0,'YDir','reverse','YLim',[0 3],'YTick',[1 2 3],'XLim',[-80 -20],'XTick',[-60 -40],...
    'XTicklabel',{'60^oS','40^oS'},'YTicklabel',{})
title({'i) GOBMs -',' OCIM-v2021'})
xlabel('latitude')

cb=colorbar;
set(cb,'Location','SouthOutside','Position',[0.35 0.37 0.35 0.03],'FontSize',12.0)
t=text(-53,5.1,{'\DeltaC_{ant} rate','(10^4 mol m^-^2 yr^-^1)'})
set(t,'FontSize',11.0,'FontWeight','bold','HorizontalAlignment','center')
