function [TTgrid,LAT,LON,DEP] = makeTTgrid(staLat, staLon, grd, model);
%  makeTTGrid      calculate travel time grid
% USAGE: [TTgrid,LAT,LON,DEP] = makeTTgrid(staLat, staLon, grd, model);
%
% Input:
%  staLat   vector of station latitudes (geographic degrees)
%  staLon   vector of station longitudes (geographic degrees)
%  grd      structure containing 3-D source grid 
%           grd.lon = vector of source longitudes (geographic degrees)
%           grd.lat = vector of source latitudes (geographic degrees)
%           grd.depth = vector of source depths (km positive down)
%
% Output:
%  TTgrid   vector of structures for each station containing containing fields 
%           TTgrid(ksta).TIM  3-D array of travel time (s)
%           TTgrid(ksta).XI   2-D array of azimuths (deg)
%           TTgrid(ksta).DELT 3-D array of horizontal distancd (km)
%  LAT      3-D array of latitudes (geographic degrees)
%  LON      3-D array of longitudes (geographic degrees)
%  DEP      3-D array of depths (km)
% 
%
%G(2).arrayName='SEQ';   G(2).array_loc=[47.99158525; -122.93713188];
%G(3).arrayName='SOOK';  G(3).array_loc=[48.38379033; -123.80361048];
%GROUP=G;
%grd.lat   = [  47.5 : 0.10  :  49.0];
%grd.lon   = [-125.0 : 0.10 : -122.0];
%grd.depth = [25:10:65];
% modified 7/20/07 to optionaly accept a new earth model
  
if nargin<4; model=''; end; % use default (PNSN P2 S-wave model) if no model is passed in

[LON,LAT,DEP]          = meshgrid(grd.lon,grd.lat,grd.depth);  % latitude and longitude of grid points
LAT2  = squeeze(LAT(:,:,1)); 
LON2  = squeeze(LON(:,:,1));
[ngrd1,ngrd2,ngrd3]= size(LON);                  % dimensions of the grids 
ngrd               = ngrd1*ngrd2*ngrd3;          % total number of grid points

for kdepth=1:ngrd3;   % loop over depths
  
  [Xray,Tray,Pray]=trace_ray(grd.depth(kdepth),model);  % calculate distance, time, ray parameter for 'model' and source depth.
  
  for ksta=1:length(staLat);  % loop over stations

    [DELTS,AZIMS] = delaz(staLat(ksta),staLon(ksta),LAT2(:),LON2(:),0);  % distance and azimuth from station to each grid point
    DELTS=DELTS*111.1;                           % station to source epicentral distance (km) for lat, lon each grid point
    tmp  = interp1(Xray,Pray,DELTS,'linear','extrap');% linearly interpolate horizontal slowness to each lat,lon grid point
    S    = reshape(tmp,  ngrd1,ngrd2);           % slowness of S wave (s/km) from each grid point
    tmp  = interp1(Xray,Tray,DELTS,'linear','extrap'); % linearly interpolate travel time curve to each lat,lon grid point
    TIM  = reshape(tmp,  ngrd1,ngrd2);           % travel time (s) to center of array from each grid point
    XI   = reshape(AZIMS,ngrd1,ngrd2);           % station to source azimuth (deg) at each lat,lon grid point
    DELT = reshape(DELTS,ngrd1,ngrd2);           % station to source azimuth (deg) at each grid point
    %    SX  = S.*sin(XI*pi/180);                % X-component of slowness from center of array to grid points
    %   SY  = S.*cos(XI*pi/180);                 % Y-component of slowness from center of array to grid points
    if kdepth==1;
      TTgrid(ksta).XI   = XI;
      TTgrid(ksta).DELT = DELT;
    end
    TTgrid(ksta).TIM(:,:,kdepth)=TIM;
    TTgrid(ksta).S(:,:,kdepth)  =S;
  end
end
return
cont_color={'b','k','r'};
for kpick=1:length(GROUP(1).picks);
  err23 = GROUP(2).picks(kpick)-GROUP(3).picks(kpick) - [GROUP(2).TIM-GROUP(3).TIM];
  err13 = GROUP(1).picks(kpick)-GROUP(3).picks(kpick) - [GROUP(1).TIM-GROUP(3).TIM];
  err_total = sqrt(err13.^2+err23.^2);
  %surf(LON,LAT,err_total);
  %view(2);
  %shading('interp');
  %caxis([0 5]);
  %colorbar
  %hold on;
  contour(LON,LAT,err_total,[.5 GROUP(1).picksQual(kpick)]); %,cont_color(pickQual(kpick,1)));
  [tmp,indmin]=min(err_total(:));
  text(LON(indmin),LAT(indmin), sprintf('%.0f', GROUP(1).picks(kpick)-GROUP(1).TIM(indmin)))
  hold on;
end
for k=1:3;
  plot(GROUP(k).array_loc(2),GROUP(k).array_loc(1),'^','MarkerFaceColor','k','markersize',15);
  text(GROUP(k).array_loc(2),GROUP(k).array_loc(1),GROUP(k).arrayName);
end
