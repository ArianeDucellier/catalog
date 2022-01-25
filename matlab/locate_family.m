grd.lon = [-123.8:.01:-122.3];
grd.lat = [39.5:.01:40.8];
grd.depth = [20:1:70];
vp_vs = 1.8;
stations = ['KHBB', 'ME08', 'ME12', 'ME27', 'ME28', 'ME37', 'ME39', 'ME41', 'ME57', 'WDC'];
lats = [40.65990, 40.222, 40.104, 40.453, 40.327, 40.285, 40.188202, 39.884998, 39.9118, 40.57988];
lons = [-123.21966, -123.305, -122.498, -123.155, -122.471, -123.653999, -123.594299, -123.361, -122.5676, -122.54113];
S_P = [14.00, 12.70, 5.50, 12.50, 6.45, 16.40, 17.10, 15.35, 3.50, 7.45];
[TTgrid, LAT, LON, DEP] = makeTTgrid(lats, lons, grd);
misfit = abs(TTgrid(1).TIM - S_P(1)) + abs(TTgrid(2).TIM - S_P(2)) + ...
         abs(TTgrid(3).TIM - S_P(3)) + abs(TTgrid(4).TIM - S_P(4)) + ...
         abs(TTgrid(5).TIM - S_P(5)) + abs(TTgrid(6).TIM - S_P(6)) + ...
         abs(TTgrid(7).TIM - S_P(7)) + abs(TTgrid(8).TIM - S_P(8)) + ...
         abs(TTgrid(9).TIM - S_P(9)) + abs(TTgrid(10).TIM - S_P(10));
[minmisfit, minindx] = min(misfit(:));
