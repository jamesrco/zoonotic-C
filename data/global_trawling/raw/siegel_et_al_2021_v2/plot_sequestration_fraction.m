% load the sequestration fraction data and the model grid
%load fseq_OCIM2_48L.mat

% variables are as follows
% 1. fseq = the fraction of carbon remaining sequestered at each ocean grid
% cell as a function of time since injection; dimensions (m x n) where m is
% the number of ocean grid cells and n is the number of years since
% injection
% 2. MASK = 3-dimensional land-sea mask for the model, with 1 = ocean and 0
% = land;dimensions (ny x nx x nz) where ny = 91 (number of grid points in
% latitudinal direction), nx = 180 (number of grid points in the
% longitudinal direction), and nz = 48 (number of grid points in the
% vertical direction)
% 3. LAT = 3-d array of the latitudes (degrees north) of the model grid
% cells; same dimensions as MASK
% 4. LON = 3-d array of the longitudes (degrees east) of the model grid
% cells; same dimensions as MASK
% 5. DEPTH = 3-d array of the depths (meters) of the model grid cells; same
% dimensions as MASK
% 6. VOL = 3-d array of grid box volumes (m^3) of the model grid cells; same
% dimensions as MASK
% 7. AREA = 3-d array of grid box areas (m^2) of the model grid cells; same
% dimensions as MASK
% 8. time = 1-d array of time (in years) since injection of CO2; range from
% 0-1000 years

% Example 1: Find the near-bottom ocean grid cells and plot the fraction of
% CO2 remaining sequestered after 100 years
% first make a 3-d array of the fraction sequestered after 100 years
[ny,nx,nz] = size(MASK);
FSEQ_100yr = 0*MASK+NaN;
t_indx = find(time==100); % find the index of the 100'th year
FSEQ_100yr(MASK==1) = fseq(:,t_indx);
% now find the near-bottom ocean grid cells
fseq_bottom_100yr = 0*MASK(:,:,1)+NaN;
TOPO = sum(MASK,3); % number of grid cells in the vertical direction
for i = 1:ny
  for j = 1:nx
    if TOPO(i,j)~=0
      fseq_bottom_100yr(i,j) = FSEQ_100yr(i,j,TOPO(i,j));
    end
  end
end
% find the bottom depth
bottom_depth = sum(VOL.*MASK,3)./AREA(:,:,1); % depth of the ocean at each water column
% plot both the fraction remaining after 100 years and the bottom depth
figure(1)
subplot(2,1,1)
pcolor(LON(:,:,1),LAT(:,:,1),fseq_bottom_100yr)
set(gca,'clim',[0 1])
colormap(parula(10))
colorbar
set(gca,'color',[.8 .8 .8])
set(gcf,'color','w')
shading('flat')
title('Fraction of CO2 injected at seafloor remaining after 100 years')
subplot(2,1,2)
pcolor(LON(:,:,1),LAT(:,:,1),bottom_depth)
set(gca,'clim',[0 5500])
colormap(parula(11))
colorbar
set(gca,'color',[.8 .8 .8])
set(gcf,'color','w')
shading('flat')
title('Depth of seafloor (m)')

% Example 2: Find a particular point in the ocean and plot the fraction
% remaining sequestered over time
% first specify a location by latitude, longitude, and depth
lat = 40; % degrees north
lon = 220; % degrees east, from 0 - 360
depth = 1000; % meters
% the following code will find the nearest model grid cell(s)
iy = find(abs(LAT(:)-lat)==min(abs(LAT(:)-lat)));
ix = find(abs(LON(:)-lon)==min(abs(LON(:)-lon)));
iz = find(abs(DEPTH(:)-depth)==min(abs(DEPTH(:)-depth)));
indx = intersect(intersect(ix,iy),iz);
% check if that grid cell(s) is land or ocean; if in ocean, plot fseq
if sum(MASK(indx))==0
  fprintf('Location is not in the ocean. Please choose another location.\n')
else
  % if in ocean, find location in ocean grid coordinates and plot
  iy2 = find(abs(LAT(MASK==1)-lat)==min(abs(LAT(MASK==1)-lat)));
  ix2 = find(abs(LON(MASK==1)-lon)==min(abs(LON(MASK==1)-lon)));
  iz2 = find(abs(DEPTH(MASK==1)-depth)==min(abs(DEPTH(MASK==1)-depth)));
  indx2 = intersect(intersect(ix2,iy2),iz2);
  figure(2)
  plot(time,mean(fseq(indx2,:)))
  xlabel('Time since injection (years)')
  ylabel('Fraction of CO2 remaining sequestered')
  title('Sequestration fraction over time')
  set(gca,'ylim',[0 1])
end
