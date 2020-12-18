function compaction_rate_mean = HL_analytic(A, Tc)
%Herron and Langway (1980) analytic solution
%Max Stevens, University of Washington
%2016
% clear all;
%T - Temperature in C
%A - Accumulation, product needs to be m W.E., e.g. South Pole is 0.08 m WE; multiply I.E. by 0.917 to get W.E.
close all;
%tic
%% Variables to change
rho_0=0.350; % Surface density in g cm^-3

dz=1; % grid spacing, m
maxdepth=100; % maximum depth of grid, m

plotter='on';
%% Model
T = 273.15+Tc; % T in K
% T = 265;
% A = 0.08;
h=0:dz:maxdepth; % grid
rho_ice=917.0;

R = 8.314;       %gas constant 
rho_i = rho_ice/1000;   %ice density (Mg m^{-3}
rho_c = 0.550;   %critical density (stage 1 - stage 2)

dh=diff(h);
dh=[dh dh(end)];

% site specific rate-constants, eqns 6a and 6b from Herron+Langway
k_0 = 11  * exp(-10160/(R*T));  
k_1 = 575 * exp(-21400/(R*T));

% Given site conditions T, A(ccumulation w.e.) and surface density,
% can calculate the density-depth profile and age.

%Stage 1 Densification
%depth of critical density, eqn 8 from Herron and Langway
h0_55 = 1/(rho_i*k_0)*(log(rho_c/(rho_i-rho_c))-...
                       log(rho_0/(rho_i-rho_0)));


Z_0 = exp(rho_i*k_0*h + log(rho_0/(rho_i-rho_0))); 

%age of critical density, eq. 9
t0_55 = 1/(k_0*A)*log((rho_i-rho_0)/(rho_i-rho_c));
rho_h0 = (rho_i * Z_0)./(1+Z_0);
t_0 =   1/(k_0*A)*log((rho_i-rho_0)./(rho_i-rho_h0));


Z_1 = exp(rho_i*k_1*(h-h0_55)/sqrt(A) +  log(rho_c/(rho_i-rho_c)));

%combine Z for Z_0 less than critical density and Z_1 greater than critical
%density
Z = [Z_0(h<h0_55) Z_1(h>h0_55)]';

% determine rho
rho_h = (rho_i * Z)./(1+Z);

t_p = 1/(k_1*sqrt(A))*log((rho_i-rho_c)./(rho_i-rho_h))+ t0_55; %Equation 11

age = [t_0(h<h0_55) t_p(h>h0_55)']'; %Eq. 12, firn age

rho_hkg=rho_h*1000; %density in kg m^-3
HLDIP=cumsum(((rho_ice-rho_hkg)./rho_ice).*dh'); %Depth-integrated porosity
HLDIPtot=HLDIP(end); %total DIP
ind=find(rho_hkg>=815.0,1);
ind2=find(rho_hkg>=801.0,1);
BCOHL=h(ind); %bubble close-off depth (where rho=815)
LIDHL=h(ind2); % lock in depth using the rho_cod minus 14 kg/m^3 method
LIZ_th=BCOHL-LIDHL;
bco_vec=815*ones(length(h),1); %for plotting

dage=age(ind);

mass=dz*rho_hkg;
sigma_MPa=9.8*cumsum(mass)/10^6;

ind=find(rho_hkg<550);

del_L = diff(h);
thickness = (max(h) - h(1:end-1));
strain = del_L./thickness;
delage = diff(age)';
strain_rate = strain./delage;
compaction_rate_mean = mean(strain_rate);

% if strcmp(plotter,'on')==1
% fig1=figure(1);
% clf;
% hold on;
% box on;
% grid on;
% plot(rho_hkg,h,'b','linewidth',2)
% % plot(bco_vec,h,'k','linewidth',1)
% set(gca,'ydir','reverse','fontsize',16)
% xlabel('Density (kg m^{-3})','fontsize',18)
% ylabel('Depth (m)','fontsize',18)
% ylim([0 100])
% xlim([300 900])
% % title('Herron and Langway (1980) firn depth/density profile','fontsize',18)
% savename='HLgeneric.png';
% print(fig1,savename,'-dpng')
% 
% fig2=figure(2);
% hold on;
% box on;
% grid on;
% plot(rho_hkg,age,'r','linewidth',2)
% % plot(bco_vec,h,'k','linewidth',1)
% set(gca,'ydir','reverse','fontsize',16)
% xlabel('Density (kg m^{-3})','fontsize',18)
% ylabel('Age (years)','fontsize',18)
% title('Herron and Langway (1980) firn depth/age profile','fontsize',18)
% text(300,1000,sprintf('Delta age = %g',dage),'fontsize',18);
% 
% fig3=figure(3);
% hold on;
% box on;
% grid on;
% plot(age,h,'r','linewidth',2)
% % plot(bco_vec,h,'k','linewidth',1)
% set(gca,'ydir','reverse','fontsize',16)
% ylabel('Depth (m)','fontsize',18)
% xlabel('Age (years)','fontsize',18)
% title('Herron and Langway (1980) firn depth/age profile','fontsize',18)
% text(300,110,sprintf('Delta age = %g',dage),'fontsize',18);
% end

%toc
end

