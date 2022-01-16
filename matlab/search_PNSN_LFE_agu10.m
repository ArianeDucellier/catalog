% % Use times of known LFEs to cut windows out of nearCross correlate a reference time event recorded at an array with continuous data to find repeating events
%clear

home_dir=pwd;
% tic
% % Use data for stations in the ascii file containing staCode tab staChannel on each line
% %station_file='good_stations.txt';   % list of stations that have good, correlated signals for the reference event
% %[sta_list.staCode,sta_list.staChannel]=textread(station_file,'%s %s','headerlines',0,'delimiter','\t');
% Filt=[1.5 8]; % 2-8?
% LFE_list=load('~/UW/IRIS_2008/LFEs/2008.073.09/LFEs2.txt','ascii'); % use LFEs whose times are described in this file
% RefTime      = LFE_list(1,1:6)'; % column vector if reference time [YY;MM;DD;HH;MM;SS.SSS]
% Ref_P_window = LFE_list(2,1:2)';
% Ref_S_window = LFE_list(2,3:4)';
% LFE_list     = LFE_list(3:end,:); % other events
% 
% % read all LFE P and S waves in the LFE_list
% windowLFE.RefTime = RefTime;
% P_start           = LFE_list(2:end,1)';
% P_start=P_start(7:end);
% Ref_P_window = [-10 20];
% windowLFE.P       = [P_start + Ref_P_window(1) ;P_start + Ref_P_window(2) ];
% %windowLFE.S       = [P_start + Ref_S_window(1) ;P_start + Ref_S_window(2) ];
% 
% sta = {'2008.073.09.00.00.XU.H05.EHZ..sac' '0803130900.B001.EHZ'	'0803130900.B007.EHZ'	'0803130900.B013.EHZ'	...
%   '0803130900.B001.EH1' '0803130900.B001.EH2' '0803130900.B013.EH1' '0803130900.B013.EH2' '0803130900.HDW.EHZ' '0803130900.GNW.BHZ' '0803130900.SQM.BHZ' ...
%   '0803130900.GNW.BHE' '0803130900.SQM.BHE' '0803130900.GNW.BHN' '0803130900.SQM.BHN' '0803130900.FACB.BHE' ... 
%   '0803130900.FACB.BHN' '0803130900.FACB.BHZ' '2008.073.09.15.58.XU.BS11.EHE..sac' '2008.073.09.15.58.XU.BS11.EHN..sac' ...
%   '2008.073.09.15.58.XU.BS11.EHZ..sac'};
% for ksta=1:length(sta);D(ksta)=coralReadSAC(sta{ksta});end
% D=coralDemean(D);
% D=coralTaper(D,-10);
% D=coralFilter(D,Filt(2),'lowpass');
% D=coralFilter(D,Filt(1),'highpass');
% %Dkeep=D;D=Dkeep([1 2 4]); tshift = [11.3439 12.6496 11.4015];
% dt = [12.6496-11.3439, 11.4015-11.3439, 12.6496-11.4015]
% % B001-H05=1.31s; B013-H05=.06s; B001-B013=1.25s; S-P(H05)=6s;
% %ind=[2 4 6 7 11 18 19 20 22 26 30 42 45 46 49 52 53 55];
% 
% figure(1);clf;
% for ksta=1:length(D);
%   D0=D(ksta);
%   for kwin=1:length(P_start);
%     optC=struct('cutType','absTime','absStartTime',timeadd(windowLFE.RefTime,P_start(kwin) + Ref_P_window(1)),...
%       'absEndTime',timeadd(windowLFE.RefTime,P_start(kwin) + Ref_P_window(2)));
%     Dwin(kwin)=coralCut(D0,optC);
%   end
%   %Dwin=Dwin(ind);
%   [stack,t] = coralBeam(Dwin, struct('n',1, 'norm',1) );
%   plot(t,stack+ksta/3);
%   text(max(t)+.1,ksta/3,sprintf('%s %s', Dwin(1).staCode,Dwin(1).staChannel))
%   hold on
% end
% xlabel('time(s)')
% title(sprintf('stacks of %d LFEs per station at %.2f-%.2fHz',length(Dwin),Filt))
% orient landscape;
% print -deps PNSNStacks.eps
% 
% % axis([8 23 0 2.8]);
% axis([6 23 0 7.4]);
% for ksta=1:21; text(23,ksta/3,sprintf('%s %s', D(ksta).staCode,D(ksta).staChannel)); end
% print -deps PNSNstacksZoom.eps
% 

% 
% IND=[1,2,4:8,15:21];
% for ksta=1:length(IND);
%   D0=D(IND(ksta));
%   for kwin=1:length(P_start);
%     optC=struct('cutType','absTime','absStartTime',timeadd(windowLFE.RefTime,P_start(kwin) + Ref_P_window(1)),...
%       'absEndTime',timeadd(windowLFE.RefTime,P_start(kwin) + Ref_P_window(2)));
%     Dwin(kwin)=coralCut(D0,optC);
%   end
%   figure(3);clf;coralPlot(Dwin);%title(sprintf('%s %s', Dwin(1).staCode,Dwin(1).staChannel))
%   %Dwin=Dwin(ind);
%   %[stack,t] = coralBeam(Dwin, struct('n',1, 'norm',1) );
%   %plot(t,stack+ksta/3);
%   title(sprintf('%s %s at %.2f-%.2fHz', Dwin(1).staCode,Dwin(1).staChannel,Filt ))
%   xlabel(num2str(IND(ksta)));ylabel(num2str(ksta))
%   axis([8 23 -inf inf]);
%   orient tall
%   eval(sprintf('print -deps LFE_%s.%s_%.2f-%.2fHz.eps', Dwin(1).staCode,Dwin(1).staChannel,Filt ))
%   %pause
% end

% load station_data.mat
% load station_data_LFE1.mat
% load LFE1_loc.mat
%%
% -122.9317 -123.1314 -122.9108 -122.9456
figure(2);clf;args.file='puget_big';
% args.map_limits=[-123.5 -122.3 47.5 48.5];
args.map_limits=[-123.5 -122 47.5 48.5];
mapp(args);
set(gca,'dataaspectratio',[1, cosd(mean(args.map_limits(3:4))), 1])
hold on;
% lons=[stacks50.staLon]; %BH1-30, BS31-57, CL58-87, GC88-114, TB-115-144
% lats=[stacks50.staLat];

lons=[mean([stacks50(1:30).staLon]),mean([stacks50(31:57).staLon]),...
  mean([stacks50(58:87).staLon]),mean([stacks50(88:114).staLon]),...
  mean([stacks50(115:141).staLon])];
lats=[mean([stacks50(1:30).staLat]),mean([stacks50(31:57).staLat]),...
  mean([stacks50(58:87).staLat]),mean([stacks50(88:114).staLat]),...
  mean([stacks50(115:141).staLat])];

staNames={stacks50([1,31,58,88,115]).staCode};



for k=1:length(lons)
  if k==1
    plot(lons(k),lats(k),'sr','markerfacecolor','r'); %plot red square for BS11 location
  else
    plot(lons(k),lats(k),'^k','markerfacecolor','k');
  end
end
%text(lons,lats,staNames,'verticalalign','bottom','FontSize',16)
%plot([-122.4,-122.8],[48.2,47.0],'r*');
vp_vs=1.8; %vp/vs
DT=[-.2:.2:.2];
% grd.lon = [-123:.01:-122.5]; grd.lat=[47.7:.01:48.1]; grd.depth=[49];[TTgrid,LAT,LON,DEP] = makeTTgrid(lats, lons, grd);
%contour(grd.lon,grd.lat,(TTgrid(2).TIM-TTgrid(1).TIM)/vp_vs,1.3+DT,'b');
%contour(grd.lon,grd.lat,(TTgrid(3).TIM-TTgrid(1).TIM)/vp_vs,0.1+DT,'b');
%contour(grd.lon,grd.lat,(TTgrid(2).TIM-TTgrid(3).TIM)/vp_vs,1.2+DT,'b');
%contour(grd.lon,grd.lat,TTgrid(1).TIM*(1-1/vp_vs),6.28+DT,'k');
%contour(grd.lon,grd.lat,TTgrid(4).TIM*(1-1/vp_vs),6.18+DT,'k');  % 6.15 to 6.18

grd.lon = [-123.5:.01:-122]; grd.lat=[47.5:.01:48.5]; grd.depth=[20:1:70];[TTgrid,LAT,LON,DEP] = makeTTgrid(lats, lons, grd);

% Data:
% B001_S - BS04_S = -0.4 s   S-S 2 - 1
% W020_S - BS04_S = 0.84 s   S-S 3 - 1
% B001_S - B001_P = 4.23 s   S-P 2 - 2
% BS04_S - BS04_P = 4.47 s   S-P 1 - 1
% B001_P - BS04_P = -0.13s   P-P 2 - 1
% TTgrid.TIM = S time
% TTgrid.TIM/vp_vs  = P time
misfit = (((TTgrid(1).TIM-TTgrid(2).TIM) - 0.74).^2 + ((TTgrid(1).TIM-TTgrid(3).TIM) + 1.30).^2 + ((TTgrid(1).TIM-TTgrid(5).TIM) - 0).^2 + ((TTgrid(2).TIM-TTgrid(3).TIM) + 2.04).^2 + ((TTgrid(2).TIM-TTgrid(5).TIM) + 0.74).^2 + ((TTgrid(3).TIM-TTgrid(5).TIM) - 1.30).^2 +...
  ((TTgrid(1).TIM-TTgrid(5).TIM)/vp_vs - 0.13).^2 +...
  ((TTgrid(1).TIM*(1-1/vp_vs)-5.27).^2) + (TTgrid(5).TIM*(1-1/vp_vs)-5.40).^2); % + ...
  %0.1*(TTgrid(2).TIM*(1-1/vp_vs)-4.23).^2 ;
%misfit = (TTgrid(1).TIM*(1-1/vp_vs)-4.47).^2;
[minmisfit,minindx]=min(misfit(:));

pred_tim=      [TTgrid(1).TIM(minindx)-TTgrid(2).TIM(minindx);
                TTgrid(1).TIM(minindx)-TTgrid(3).TIM(minindx)
				TTgrid(1).TIM(minindx)-TTgrid(5).TIM(minindx)
				TTgrid(2).TIM(minindx)-TTgrid(3).TIM(minindx)
				TTgrid(2).TIM(minindx)-TTgrid(5).TIM(minindx)
				TTgrid(3).TIM(minindx)-TTgrid(5).TIM(minindx)
                (TTgrid(1).TIM(minindx)-TTgrid(5).TIM(minindx))/vp_vs;
                TTgrid(1).TIM(minindx)*(1-1/vp_vs);
                TTgrid(5).TIM(minindx)*(1-1/vp_vs)];





% pred_tim=    [TTgrid(2).TIM(minindx)-TTgrid(1).TIM(minindx) ;
%              (TTgrid(3).TIM(minindx)-TTgrid(1).TIM(minindx))  ; 
%              (TTgrid(2).TIM(minindx)-TTgrid(1).TIM(minindx))/vp_vs ;
%              TTgrid(1).TIM(minindx)*(1-1/vp_vs);
%              TTgrid(2).TIM(minindx)*(1-1/vp_vs) ];
 obs_tim= [.74; -1.3; 0; -2.04; -.74; 1.3; .13; 5.27; 5.40]; %[-.4; .84 ; -.13 ; 4.47 ; 4.23 ];         
 format bank; [obs_tim,pred_tim,obs_tim-pred_tim]
           
% datamisfit= -[TTgrid(2).TIM(minindx)-TTgrid(1).TIM(minindx) + 0.4  ;
%              (TTgrid(3).TIM(minindx)-TTgrid(1).TIM(minindx)) - 0.84 ; 
%              (TTgrid(2).TIM(minindx)-TTgrid(1).TIM(minindx))/vp_vs - .13;
%              TTgrid(1).TIM(minindx)*(1-1/vp_vs)-4.47 ;
%              TTgrid(2).TIM(minindx)*(1-1/vp_vs)-4.23];
% disp('Data misfit ')
% disp(datamisfit)

[LON(minindx),LAT(minindx),DEP(minindx)]
% plot(LON(minindx),LAT(minindx),'.k','markersize',30);
dep_indx=find(DEP(minindx)==grd.depth)
contour(grd.lon,grd.lat,squeeze(misfit(:,:,dep_indx)),minmisfit+[0:.2:1]);
grd2.lon = [-123.5:.01:-122]; grd2.lat=[47.5:.01:48.5]; grd2.depth=[20:1:70];[TTgrid,LAT,LON,DEP] = makeTTgrid(lats, lons, grd2);
contour(grd2.lon,grd2.lat,squeeze(misfit(:,:,dep_indx)),minmisfit+[0:.2:1]);

% load beam_locs.mat
% plot(ABHI(k,7),ABHI(k,8))

h=gca;
set(h,'FontSize',16);

xlabel('Longitude','FontSize',16);
ylabel('Latitude','FontSize',16);
del='{\delta}';
% title(sprintf('LFE from %sP times and two S-P times; depth: %2.1fkm;
% 2008/03/13; 09:00-10:00',del,DEP(minindx)),'FontSize',16);
% title(sprintf('LFE from 3 %sP times, 1 %sS time, and 2 S-P times; depth: %2.1fkm',del,del,DEP(minindx)),'FontSize',16);
title(sprintf('LFE from 3 %sP times, 1 %sS time, and 2 S-P times',del,del),'FontSize',16);
% orient tall;
%print -depsc LFE_locMap.eps

figure(4);clf;
% plot(LON(minindx),-DEP(minindx),'.k','markersize',30);
hold on
% plot([LON(1) LON(numel(LON))],[-DEP(minindx) -DEP(minindx)],'--k');
hold on
for k=1:length(lons)
  if k==1
    plot(lons(k),0,'sr','markerfacecolor','r'); %plot red square for BS11 location
  else
    plot(lons(k),0,'^k','markerfacecolor','k');
  end
end
lat_indx=find(LAT(minindx)==grd.lat);
contour(grd.lon,-grd.depth,squeeze(misfit(lat_indx,:,:))',minmisfit+[0:.2:1]);
contour(grd2.lon,-grd2.depth,squeeze(misfit(lat_indx,:,:))',minmisfit+[0:.2:1]);
h=gca;
set(h,'FontSize',16);
xlabel('Longitude','FontSize',16);
ylabel('Depth (km)','FontSize',16);
title('LFE hypocenter and plate interface','FontSize',16);
hold on
home=pwd;
% eval('cd ~/UW/fromWECH/');
[d,x,y,dist]=plate_slice(args.map_limits(1),LAT(minindx),args.map_limits(2),LAT(minindx));
plot(x,d,'b')
% disp(dist)
ylim([-max(dist) 0]) %resize plot so x and y axes are of same length (dist) in km
eval(sprintf('cd %s',home));
orient tall;
% print -depsc LFE_Lon-Dep.eps


% figure(5);clf;
% plot(LAT(minindx),-DEP(minindx),'.k','markersize',30);
% hold on
% plot(lats,0,'^k','markerfacecolor','k');
% lon_indx=find(LON(minindx)==grd.lon);
% contour(grd.lat,-grd.depth,squeeze(misfit(:,lon_indx,:))',minmisfit+[0:.2:1]);
% xlabel('Latitude');
% ylabel('Depth (km)');
% hold on
% %eval('cd ~/UW/fromWECH/');
% [d,x,y,dist]=plate_slice(LON(minindx),args.map_limits(3),LON(minindx),args.map_limits(4));
% plot(y,d,'b')
% % disp(dist)
% ylim([-max(dist) 0]) %resize plot so x and y axes are of same length (dist) in km
% eval(sprintf('cd %s',home));
% print -depsc LFE_Lat-Dep.eps

%%
figure(6);clf;
del='{\delta}';
% title(sprintf('LFE from %sP times and two S-P times; depth: %2.1fkm;  2008/03/13; 09:00-10:00',del,DEP(minindx)));

subplot(2,1,1)
args.file='puget_big';
args.map_limits=[-123.5 -122 47.5 48.5];
mapp(args);set(gca,'dataaspectratio',[1,cosd(48),1]);
% set(gca,'dataaspectratio',[1, cosd(mean(args.map_limits(3:4))), 1])
hold on;
% plot(LON(minindx),LAT(minindx),'.k','markersize',30);
% for k=1:length(lons)
%   if k==1
%     plot(lons(k),lats(k),'sr','markerfacecolor','r'); %plot red square for BS11 location
%   else
%     plot(lons(k),lats(k),'^k','markerfacecolor','k');
%   end
% end
plot(lons,lats,'^k','markerfacecolor','k');
% text(lons,lats,staNames,'verticalalign','bottom','FontSize',16)
contour(grd.lon,grd.lat,squeeze(misfit(:,:,dep_indx)),minmisfit+[0:.2:1])
contour(grd2.lon,grd2.lat,squeeze(misfit(:,:,dep_indx)),minmisfit+[0:.2:1])

% load beam_locs.mat
% plot(ABHI(k,7),ABHI(k,8))
% load slab
% contour (slablon, slablat, slabdepth, [-80: 10: -20]);
h=gca;
set(h,'FontSize',16);
xlabel('Longitude','FontSize',16);
ylabel('Latitude','FontSize',16);
title('LFEs Location - Map View','FontSize',16);
orient tall

subplot(2,1,2)
% plot(LON(minindx),-DEP(minindx),'.k','markersize',30);
hold on
plot(lons,0,'^k','markerfacecolor','k');
% plot(lons(1),0,'sr','markerfacecolor','r');
lat_indx=find(LAT(minindx)==grd.lat);
contour(grd.lon,-grd.depth,squeeze(misfit(lat_indx,:,:))',minmisfit+[0:.2:1]);
contour(grd2.lon,-grd2.depth,squeeze(misfit(lat_indx,:,:))',minmisfit+[0:.2:1]);
h=gca;
set(h,'FontSize',16);
xlabel('Longitude','FontSize',16);
ylabel('Depth (km)','FontSize',16);
title('LFEs Location - Cross Section','FontSize',16);
hold on
home=pwd;
% eval('cd ~/UW/fromWECH/');
[d,x,y,dist]=plate_slice(args.map_limits(1),LAT(minindx),args.map_limits(2),LAT(minindx));
plot(x,d,'b')
disp(max(dist))
ylim([-max(dist) 0]) %resize plot so x and y axes are of same length (dist) in km
xlim([-123.5 -122]);
grid on
set(gca,'dataaspectratio',[1,111.1*cosd(48),1]);
orient tall

% cd ~/Desktop/LFEs/agu_2009_figs/
% print('-depsc','map_xsection_LFEs.eps')
% cd(home_dir)
% eval(sprintf('cd %s',home));