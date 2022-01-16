function [X,T,P,vs,vr]=trace_ray(depth, model);
%   trace_ray    calculate travel time curve by tracing rays
% USAGE: [X,T,P,vs,vr]=trace_ray(depth, model);
%
%  Input: 
%  depth is earthquake depth (km)
%  model is an optional velocity model (Nx2 array) with
%  velocity (km/s) in first column and depth (km) in second
%  velocity is a layered model, each line of the model gives the
%  velocity of a layer and depth at the top of the layer
%  use default (PNSN Puget Sound (P2)) model if model='' or if 
%  there is only one input argument
% 
% Output:
%  X horizontal distance (km)
%  T travel time (s)
%  P ray parameter (dT/dX s/km)
%  vs velocity at source (km/s)
%  vr velocity at receiver (surface) (km/s)
%
% example:
% clear; figure(20);clf; depth=[25 35 45 55]; colr = {'k' 'm' 'r' 'b'} ;
% for k=1:4; [X,T,P]=trace_ray(depth(k));hold on;h(k)=plot(X,T-X/4.5033,['.-' colr{k}]);end;axis([0 300 3 13]); grid on;
% legend(h,'25 km','35 km','45 km','55 km');xlabel('Range (km)'); ylabel('Travel Time - Range/4.5033 (s)')

%clear;close all;  % clear matlab workspace and close all figures

h0=depth;   % source depth

if nargin<2 || length(model)==0;  % use default model (PNSN model P2 converted to a shear wave model)

  vel_P2=[
   5.40  0.0
   6.38  4.0
   6.59  9.0
   6.73 16.0
   6.86 20.0
   6.95 25.0
   7.80 51.0;
   8.00 81.0];
 
  model = [vel_P2(:,1)/sqrt(3) ,  vel_P2(:,2)];
  
end

h  = model(:,2);    % column vector of depths corresponding to the tops of layers (km)
vl = model(:,1);    % velocity within each layer (km/s)

% add a new layer boundary at the source depth without changing the model
if length(find(h0==h))==0;
  i=max(find(h0>=h));
  h = [h(1:i); h0; h(i+1:end)];
  vl= vl([1:i i:end]);
else
  i=find(h0==h);
end
i0=i;
v0= vl(i0);
vs=v0;
vr=vl(1); %velocity at surface (receiver).

%   take_off_layer{1}  defines rays going up
%   take_off_layer{k}  (k>1) defines rays going down and reflecting critically off the layer boundary that is k-1 interfaces below the source

k0=1;
take_off_layer{k0} = -[0.01:5:80 81:2:87 88 89 89.9]*pi/180; % column vector of desired take-off angles (radians) < 0 means up-going ray

for k=i0+1:length(vl)-1;   % define rays reflecting critically off each interface below the source
  k0=k0+1;
  pmin = 1/vl(k+1);
  pmax = 1/vl(k);
  prange = pmin+(pmax-pmin)*[0:.05:.95 .96:.01:.99 .993:.002:.999];  
  prange = pmin+(pmax-pmin)*[.9:-.1:.1 .01];    
  take_off_layer{k0} = asin(v0*prange);
end

% combine the take-off angles into on vector 
take_off=[];
for k=1:1; % end loop at 1 for upgoing rays only and at k0 for all rays
  take_off=[take_off take_off_layer{k}]; 
end

k=find(take_off>0);
if length(k)>0;
  take_off(k) = max(take_off(k),asin(v0/vl(end))+eps );   % force all rays to turn with the model.
end


%h=[0:5.5:6310]';       % column vector of depths in the earth (km)
%r=6371-h;              % vector of equivalent radii
%v=prem(r);             % column vector of P-wave velocities (km/s) at depths h
%vl=(v(1:end-1)+v(2:end))/2; % mean velocity within each layer
dh=diff(h);                  % thickness of each layer (km)

%take_off=[30:1:85]'*pi/180;% column vector of desired take-off angles (radians)
%take_off=[26:.1:32]'*pi/180% column vector of desired take-off angles (radians)
%figure(1);clf;hold on     % open a new plot window
for l=1:length(take_off)   % loop over each ray
  clear dz dx dt x t i interface; % clear vectors containing ray increments
  k=1;
  interface(k) = find(h0==h);
  i(k)         = take_off(l);      % start a new ray
  x(k)         = 0;
  t(k)         = 0;
  while ~(interface(k)==1 & i(k)<0) % stop if interface is 1 and ray is going up.
    k=k+1;
    if i(k-1)<0;   % if ray going up then layer=interface-1
      layer = interface(k-1)-1; 
      layer_next = layer-1;
      interface(k) = interface(k-1)-1;
    else         
      layer = interface(k-1); % if ray going down that layer=interface
      layer_next = layer+1;
      interface(k) = interface(k-1)+1;
    end
    dz = h(interface(k)) - h(interface(k-1));
    dx = dz*tan(i(k-1)); % determine horizontal distance traveled in this layer
    dt = abs(dz/(vl(layer)*cos(i(k-1)))); % calculte travel time through this layer
    x(k)=x(k-1)+dx;
    t(k)=t(k-1)+dt;
    if layer_next == 0; break;end
    sini  = vl(layer_next)*sin(i(k-1))/vl(layer);% solve Snell's Law for ray angle in next layer
    if sini>1;                  % ray must reflect, change sign or ray angle 
      i(k) = -i(k-1);
    else
      i(k) = asin(sini);        % calculate ray angle in next layer
    end
    if k>20; break;end
  end
  %plot(x,-h(interface),'.-');axis('equal');hold on
  T(l,1)=t(end);
  X(l,1)=x(end);
end
P=sin(abs(take_off(:)))/v0; % horizontal slowness
%figure(2);subplot(211);plot(X,T,'.-r');axis([0 600 0 100]); subplot(212);plot(X,P,'.-r');axis([0 600 0 .3])
return

