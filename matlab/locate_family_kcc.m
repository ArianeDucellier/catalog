clear;
PPweights=[0 .1 .3 .5 .75 1 2 10];
PPweight=0.5;  
% first section finds location for S-P times and plots data and maps
% second section finds location for P-P times down weighted by PPweight and
% S-P times and makes plots for S-P times at the location for PP and SP
grd.lon = [-123.8:.01:-122.3];
grd.lat = [39.5:.01:40.8];
grd.depth = [20:1:70];
vp_vs = 1.8;

% Using S minus P times
stations = {'KHBB', 'ME08', 'ME12', 'ME28', 'ME37', 'ME39', 'ME41', 'ME57', 'WDC'};
lats = [40.65990, 40.222, 40.104, 40.327, 40.285, 40.188202, 39.884998, 39.9118, 40.57988];
lons = [-123.21966, -123.305, -122.498, -122.471, -123.653999, -123.594299, -123.361, -122.5676, -122.54113];
S_P = [15.65, 12.70, 2.20, 6.75, 17.90, 16.85, 15.35, 3.85, 7.65];
S_P = S_P(:);
[TTgrid, LAT, LON, DEP] = makeTTgrid(lats, lons, grd);
misfit=0;
for k=1:9;
  misfit=misfit + abs((vp_vs - 1) * TTgrid(k).TIM - S_P(k));
end
misfit_SP=misfit;
[minmisfit, minindx] = min(misfit(:));
LAT0 = LAT(minindx);
LON0 = LON(minindx);
DEP0 = DEP(minindx);
fprintf('%f %f %f\n', LON0, LAT0, DEP0);
% Plot cross section of misfit with fixed latitude
figure(1); clf;
k = find(grd.lat == LAT0);
contour(reshape(LON(k, : , :), 151, 51), reshape(DEP(k, :, :), 151, 51), reshape(misfit(k, :, :), 151, 51),'ShowText','on');
hold on;
plot(LON0, DEP0, '*');
xlabel('Longitude', 'FontSize', 16);
ylabel('Depth', 'FontSize', 16);
% Plot cross section of misfit with fixed longitude
figure(2); clf;
k = find(grd.lon == LON0);
contour(reshape(LAT(:, k, :), 131, 51), reshape(DEP(:, k, :), 131, 51), reshape(misfit(:, k, :), 131, 51),'ShowText','on');
hold on;
plot(LAT0, DEP0, '*');
xlabel('Latitude', 'FontSize', 16);
ylabel('Depth', 'FontSize', 16);
% Plot map of misfit with fixed depth
figure(3); clf;
k = find(grd.depth == DEP0);
contour(reshape(LON(:, : , k), 131, 151), reshape(LAT(:, :, k), 131, 151), reshape(misfit(:, :, k), 131, 151),'ShowText','on');
hold on;
plot(LON0, LAT0, '*');
xlabel('Longitude', 'FontSize', 16);
ylabel('Latitude', 'FontSize', 16);
%%
Dist = delaz(LAT0,LON0,lats,lons,0)*111.1; %epicentral distance (km)
k1 = find(grd.lat == LAT0);k2 = find(grd.lon == LON0);k3 = find(grd.depth == DEP0);
disp([k1,k2,k3]);
for k=1:9; S_P0(k,1)=(vp_vs - 1) * TTgrid(k).TIM(k1,k2,k3) ; end;  % get predicted S_P times
figure(10);clf; plot( Dist, S_P0,'bo',Dist,  S_P, 'kd'); xlabel('Range (km)'); ylabel('Travel Time (s)'); 
title('Predicted (o) and Observed (d)S-P times')
%%

% Using P-P times
stations = {'KBN', 'KHBB', 'LRB', 'ME08', 'ME12', 'ME28', 'ME37', 'ME39', 'ME41', 'ME57', 'WDC'};
lats = [39.89237, 40.65990, 40.14323, 40.222, 40.104, 40.327, 40.285, 40.188202, 39.884998, 39.9118, 40.57988];
lons = [-123.19503, -123.21966, -122.55772, -123.305, -122.498, -122.471, -123.653999, -123.594299, -123.361, -122.5676, -122.54113];
PP = [  1.95,  -7.15,   1.8 ,  -5.3 ,  -6.2 ,   4.6 ,   5.75,   2.4 , ...
        -4.5 ,  -4.  , -10.5 ,  -1.75,  -8.35,  -9.5 ,   4.05,   0.95, ...
         1.35,  -7.8 ,  -5.05,   7.8 ,   3.3 ,   1.  ,  13.15,  12.45, ...
        10.95,   2.25,   3.2 ,  -8.6 ,  -6.85,   5.3 ,   4.6 ,   3.6 , ...
        -5.1 ,  -3.3 ,  -2.3 ,  10.4 ,   9.25,  10.65,  -0.55,   0.85, ...
        13.1 ,  12.  ,  10.45,   1.75,   2.25,  -1.7 ,  -2.2 , -11.9 , ...
       -10.5 ,  -1.5 , -10.2 ,  -9.3 ,  -8.7 ,  -6.4 ,   0.9 ];
     PP=PP(:);
[TTgrid, LAT, LON, DEP] = makeTTgrid(lats, lons, grd);
misfit = 0;
k = 0;
for i = 1:10
    for j = (i+1):11
        k = k + 1;
        misfit = misfit + abs(TTgrid(j).TIM - TTgrid(i).TIM - PP(k));
    end
end
misfit_PP = misfit;
[minmisfit, minindx] = min(misfit(:));
LAT0 = LAT(minindx);
LON0 = LON(minindx);
DEP0 = DEP(minindx);
fprintf('%f %f %f\n', LON0, LAT0, DEP0);
% Plot cross section of misfit with fixed latitude
figure(4); clf;
k = find(grd.lat == LAT0);
contour(reshape(LON(k, : , :), 151, 51), reshape(DEP(k, :, :), 151, 51), reshape(misfit(k, :, :), 151, 51),'ShowText','on');
hold on;
plot(LON0, DEP0, '*');
xlabel('Longitude', 'FontSize', 16);
ylabel('Depth', 'FontSize', 16);
% Plot cross section of misfit with fixed longitude
figure(5); clf;
k = find(grd.lon == LON0);
contour(reshape(LAT(:, k, :), 131, 51), reshape(DEP(:, k, :), 131, 51), reshape(misfit(:, k, :), 131, 51),'ShowText','on');
hold on;
plot(LAT0, DEP0, '*');
xlabel('Latitude', 'FontSize', 16);
ylabel('Depth', 'FontSize', 16);
% Plot map of misfit with fixed depth
figure(6); clf;
k = find(grd.depth == DEP0);
contour(reshape(LON(:, : , k), 131, 151), reshape(LAT(:, :, k), 131, 151), reshape(misfit(:, :, k), 131, 151),'ShowText','on');
hold on;
plot(LON0, LAT0, '*');
xlabel('Longitude', 'FontSize', 16);
ylabel('Latitude', 'FontSize', 16);

%%
Dist = delaz(LAT0,LON0,lats,lons,0)*111.1; %epicentrl distance (km)
k1 = find(grd.lat == LAT0);k2 = find(grd.lon == LON0);k3 = find(grd.depth == DEP0);
disp([k1,k2,k3]);
k = 0;
for i = 1:10
    for j = (i+1):11
        k = k + 1;
        PP0(k,1) = TTgrid(j).TIM(k1,k2,k3) - TTgrid(i).TIM(k1,k2,k3);
        DIST(k,:)=[Dist(i),Dist(j)];
    end
end
DDIST=DIST(:,2)-DIST(:,1);
figure(11);clf; plot( DDIST, PP0,'bo',DDIST,  PP, 'kd'); xlabel('Range difference (km)'); ylabel('Travel Time (s)'); 
title('Predicted (o) and Observed (d) P-P Times')
figure(12); clf; histogram(PP-PP0,20); title('Histogram of PP residuals');xlabel('Time (s)');
format bank; disp([LON0,LAT0,DEP0, PPweight]); format short

misfit=misfit*PPweight+misfit_SP;  % Simultaneous solution for PP and SP data down weighting PP data
[minmisfit, minindx] = min(misfit(:));
LAT0 = LAT(minindx);
LON0 = LON(minindx);
DEP0 = DEP(minindx);
fprintf('%f %f %f\n', LON0, LAT0, DEP0);
%% results below is from running this code 8 times using the 
%%weithts in the 4th column to get the location in the first 3 columns
Results=[
       -122.92         40.30         21.00             0
       -122.81         40.27         20.00          0.10
       -122.75         40.25         20.00          0.30
       -122.74         40.25         21.00          0.50
       -122.71         40.24         27.00          0.75
       -122.69         40.23         30.00          1.00
       -122.65         40.23         39.00          2.00
       -122.56         40.21         54.00         10.00];
scal=[111.1*cosd(39),111.1,1];
figure(20);clf;
labls={'Longitude','Latitude','Depth'};
for k=1:3; subplot(3,1,k);plot(Results(:,k),'.-');
  limits=mean(Results(:,k))+[-20,20]/scal(k);
  ylim(limits);
  ylabel(labls{k});
end
xlabel('Left is S-P only, right is 99% P-P: weights: 0 .1 .3 .5 .75 1 2 10')

figure(21);clf; plot(Results(:,1),Results(:,2),'*');set(gca,'dataAspectRatio',[1,cosd(39),1])
hold on; plot(lons,lats,'^')
