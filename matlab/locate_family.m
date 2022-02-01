grd.lon = [-123.8:.01:-122.3];
grd.lat = [39.5:.01:40.8];
grd.depth = [20:1:70];
vp_vs = 1.8;

% Using S minus P times
stations = ['KHBB', 'ME08', 'ME12', 'ME28', 'ME37', 'ME39', 'ME41', 'ME57', 'WDC'];
lats = [40.65990, 40.222, 40.104, 40.327, 40.285, 40.188202, 39.884998, 39.9118, 40.57988];
lons = [-123.21966, -123.305, -122.498, -122.471, -123.653999, -123.594299, -123.361, -122.5676, -122.54113];
S_P = [14.00, 12.70, 5.50, 12.50, 6.45, 16.40, 17.10, 15.35, 3.50, 7.45];
[TTgrid, LAT, LON, DEP] = makeTTgrid(lats, lons, grd);
misfit = abs(TTgrid(1).TIM - S_P(1)) + abs(TTgrid(2).TIM - S_P(2)) + ...
         abs(TTgrid(3).TIM - S_P(3)) + abs(TTgrid(4).TIM - S_P(4)) + ...
         abs(TTgrid(5).TIM - S_P(5)) + abs(TTgrid(6).TIM - S_P(6)) + ...
         abs(TTgrid(7).TIM - S_P(7)) + abs(TTgrid(8).TIM - S_P(8)) + ...
         abs(TTgrid(9).TIM - S_P(9));
[minmisfit, minindx] = min(misfit(:));
LAT0 = LAT(minindx);
LON0 = LON(minindx);
DEP0 = DEP(minindx);
fprintf('%f %f %f', LON0, LAT0, DEP0);
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

% Using P-P times
stations = ['KBN', 'KHBB', 'LRB', 'ME08', 'ME12', 'ME28', 'ME37', 'ME39', 'ME41', 'ME57', 'WDC'];
lats = [39.89237, 40.65990, 40.14323, 40.222, 40.104, 40.327, 40.285, 40.188202, 39.884998, 39.9118, 40.57988];
lons = [-123.19503, -123.21966, -122.55772, -123.305, -122.498, -122.471, -123.653999, -123.594299, -123.361, -122.5676, -122.54113];
PP = [  1.95,  -7.15,   1.8 ,  -5.3 ,  -6.2 ,   4.6 ,   5.75,   2.4 , ...
        -4.5 ,  -4.  , -10.5 ,  -1.75,  -8.35,  -9.5 ,   4.05,   0.95, ...
         1.35,  -7.8 ,  -5.05,   7.8 ,   3.3 ,   1.  ,  13.15,  12.45, ...
        10.95,   2.25,   3.2 ,  -8.6 ,  -6.85,   5.3 ,   4.6 ,   3.6 , ...
        -5.1 ,  -3.3 ,  -2.3 ,  10.4 ,   9.25,  10.65,  -0.55,   0.85, ...
        13.1 ,  12.  ,  10.45,   1.75,   2.25,  -1.7 ,  -2.2 , -11.9 , ...
       -10.5 ,  -1.5 , -10.2 ,  -9.3 ,  -8.7 ,  -6.4 ,   0.9 ];
[TTgrid, LAT, LON, DEP] = makeTTgrid(lats, lons, grd);
misfit = 0;
k = 0;
for i = 1:10
    for j = (i+1):11
        k = k + 1;
        misfit = misfit + abs(TTgrid(j).TIM - TTgrid(i).TIM - PP(k));
    end
end
[minmisfit, minindx] = min(misfit(:));
LAT0 = LAT(minindx);
LON0 = LON(minindx);
DEP0 = DEP(minindx);
fprintf('%f %f %f', LON0, LAT0, DEP0);
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
