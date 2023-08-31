
%%
prefix='../ready_to_plot_netcdfs/';
filename='figure13_Cantscatter.nc';

delta_cant_SO=ncread([prefix filename],'delta_cant_SO');
cant_1994_SO=ncread([prefix filename],'cant_1994_SO');
cantflux_SO=ncread([prefix filename],'cantflux_SO');
sos_SPSS_STSS=ncread([prefix filename],'sos_SPSS_STSS');
delta_cant_gruber=ncread([prefix filename],'delta_cant_gruber');
cant_1994_sabine=ncread([prefix filename],'cant_1994_sabine');
sos_SPSS_STSS_somffn=ncread([prefix filename],'sos_SPSS_STSS_somffn');

%% Regression statistics (excluding ROMS-SouthernOcean-ETHZ and OCIM-v2021 products)

cant_1994_gobm=cant_1994_SO(1:11)';
cantflux_gobm=cantflux_SO(1:11)';
delta_cant_gobm=delta_cant_SO(1:11)';
sos_gobm=sos_SPSS_STSS(1:11)';

%necessary for the computation of regression statistics:
cant_regr=[cant_1994_gobm' ones(size(cant_1994_gobm))'];
cantfl_regr=[cantflux_gobm' ones(size(cantflux_gobm))'];
sos_regr=[sos_gobm' ones(size(sos_gobm))'];

[B,BINT,R,RINT,STATS] = regress(delta_cant_gobm',cant_regr);
cantdiff_lin=polyval(B,cant_1994_gobm); 
R2_cant1994=STATS(1);
p_cant1994=STATS(3);

[B,BINT,R,RINT,STATS] = regress(delta_cant_gobm',cantfl_regr);
cantfl_lin=polyval(B,cantflux_gobm); 
R2_cantfl=STATS(1);
p_cantfl=STATS(3);

[B,BINT,R,RINT,STATS] = regress(delta_cant_gobm',sos_regr);
sos_lin=polyval(B,sos_gobm); 
R2_sos=STATS(1);
p_sos=STATS(3);


% compute uncertainty envelopes for observation-based products (18.9% for eMLR(C*) and 20% for Sabine et al., (2004))
unc=delta_cant_gruber*0.189;
unc_up=delta_cant_gruber+unc;
unc_bot=delta_cant_gruber-unc;

unc94=cant_1994_sabine*0.2;
unc_up_94=cant_1994_sabine+unc94;
unc_bot_94=cant_1994_sabine-unc94;

%%
models_title={'CESM-ETHZ','MRI-ESM2-1','NorESM-OC1.2','NEMO-plankTOM12',...
    'CCSM-WHOI','CNRM-ESM2-1','EC-Earth3','FESOM-REcoM-LR','ORCA025-GEOMAR','ORCA1-LIM3-PISCES','MPIOM-HAMOCC', ...
    'ROMS-SouthernOcean-ETHZ','OCIM-v2021','Observation-based products'};

  
figure(1)

clf
cf=255;

colorscale={[1 0 0],[1 0 0],[1 0 0],[1 0 0],...
    [0 0 1],[0 0 1],[0 0 1],[0 0 1],[0 0 1],[0 0 1],[0 0 1],...
    [0.4 0.4 0.4],[0 0 0],[0 0 0]};

mark={'o','d','+','>',...
    '*','d','+','s','<','h','x',...
    'v','p'}

subplot(2,2,1,'Position',[0.1 0.6 0.28 0.32])

for m=1:13
h=scatter(cant_1994_SO(m),delta_cant_SO(m),'Marker',mark{m},'MarkerEdgeColor',colorscale{m},'LineWidth',1.5)
h.SizeData = 60; 
hold on
end

h=scatter(cant_1994_sabine,delta_cant_gruber,'Marker','o','LineWidth',2.5,'MarkerEdgeColor','k')
h.SizeData = 60; 

hold on
fill([3.5,37.5,37.5,3.5],[unc_up,unc_up,unc_bot,unc_bot],[0.5,0.5,0.5],'edgecolor','none')
grid on
set(gca, 'Layer', 'top');
hold off
alpha(0.2)


hold on
fill([unc_bot_94,unc_up_94,unc_up_94,unc_bot_94],[3.95/13 3.95/13 9.9/13 9.9/13],[0.5,0.7,0.1],'edgecolor','none')
grid on
set(gca, 'Layer', 'top');
hold off
alpha(0.1)
hold on


hold on 
for m=1:13
h=scatter(cant_1994_SO(m),delta_cant_SO(m),'Marker',mark{m},'MarkerEdgeColor',colorscale{m},'LineWidth',1.5)
h.SizeData = 60; 
hold on
end

h=scatter(cant_1994_sabine,delta_cant_gruber,'Marker','o','LineWidth',2.5,'MarkerEdgeColor','k')
h.SizeData = 60; 

hold on
plot(cant_1994_gobm,cantdiff_lin,'LineWidth',0.7,'Color','k')
hold on


set(gca,'FontSize',8.0,'XGrid','on','YGrid','on','YLim',[3.95/13 9.9/13],'XLim',[3.5 37.5])
xlabel('C_{ant} inventory in 1994 [PgC ]')
ylabel({'\DeltaC_{ant} [PgC year^-^1]'})

t=text(5,0.71,['R^2 = ' num2str(round(R2_cant1994*100)/100)]);
set(t,'FontSize',12.0)
t=text(5,0.66,['p = ' num2str(round(p_cant1994*100000)/100000)]);
set(t,'FontSize',12.0)


t=text(0,0.3,'a');
set(t,'FontSize',20.0,'FontWeight','bold')

%
subplot(2,2,2,'Position',[0.47 0.6 0.28 0.32])

for m=1:13
h=scatter(cantflux_SO(m),delta_cant_SO(m),'Marker',mark{m},'MarkerEdgeColor',colorscale{m},'LineWidth',1.5)
h.SizeData = 60; 

hold on
end


hold on
plot(cantflux_gobm,cantfl_lin,'LineWidth',0.7,'Color','k')

set(gca,'FontSize',8.0,'XGrid','on','YGrid','on','YLim',[3.95/13 9.9/13],'XLim',[0.55 0.99])
ylabel('\DeltaC_{ant} [PgC year^-^1]')
xlabel({'C_{ant} air-sea fluxes [PgC year^-^1]'})

t=text(0.57,0.71,['R^2 = ' num2str(round(R2_cantfl*100)/100)]);
set(t,'FontSize',12.0)
t=text(0.57,0.66,['p = ' num2str(round(p_cantfl*10000)/10000)]);
set(t,'FontSize',12.0)

t=text(0.5,0.3,'b');
set(t,'FontSize',20.0,'FontWeight','bold')


subplot(2,2,3,'Position',[0.1 0.1 0.28 0.32])

for m=1:13
h=scatter(sos_SPSS_STSS(m),delta_cant_SO(m),'Marker',mark{m},'MarkerEdgeColor',colorscale{m},'LineWidth',1.5)
h.SizeData = 60; 

hold on
end


h=scatter(sos_SPSS_STSS_somffn,delta_cant_gruber,'Marker','o','LineWidth',2.5,'MarkerEdgeColor','k')
h.SizeData = 60; 

hold on
plot(sos_gobm,sos_lin,'LineWidth',0.7,'Color','k')

set(gca,'FontSize',8.0,'XGrid','on','YGrid','on','YLim',[3.95/13 9.9/13],'XLim',[33.8 34.7])
ylabel('\DeltaC_{ant} [PgC year^-^1]')
xlabel({'SSS'})

t=text(33.83,0.71,['R^2 = ' num2str(round(R2_sos*100)/100)]);
set(t,'FontSize',12.0)
t=text(33.83,0.66,['p = ' num2str(round(p_sos*100000)/100000)]);
set(t,'FontSize',12.0)
t=text(33.6,0.3,'c');
set(t,'FontSize',20.0,'FontWeight','bold')

t=text(35.7,0.8,'GOBMs high')
set(t,'FontSize',12.0,'FontWeight','bold','Color','r')

t=text(35.7,0.54,'GOBMs low')
set(t,'FontSize',12.0,'FontWeight','bold','Color','b')

l1=legend(models_title) 
set(l1,'Position',[0.65 0.3 0.4 0.4],'FontSize',10.0)


%%
legend('boxoff')


