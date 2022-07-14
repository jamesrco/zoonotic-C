% gen_fracs_to_constrain_trawlCO2.m
% Created June 7, 2022 by Jamie Collins, jcollins@edf.org, based on the script 
% plot_sequestration_fraction.m provided by Siegel et al. at
%  https://doi.org/10.6084/m9.figshare.15228690.v2
% Note: Assumes user has loaded "fseq_OCIM2_48L.mat" into memory 

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

% Find the near-bottom ocean grid cells and plot the fraction of
% CO2 remaining sequestered after 1, 25, 50, and then 100 years ... plus
% 1000 y as a check

% make a 3-d array of the fraction sequestered after 1 year
[ny,nx,nz] = size(MASK);
FSEQ_1yr = 0*MASK+NaN;
t_indx = find(time==1); % find the index of the 1st year
FSEQ_1yr(MASK==1) = fseq(:,t_indx);
% now find the near-bottom ocean grid cells
fseq_bottom_1yr = 0*MASK(:,:,1)+NaN;
TOPO = sum(MASK,3); % number of grid cells in the vertical direction
for i = 1:ny
  for j = 1:nx
    if TOPO(i,j)~=0
      fseq_bottom_1yr(i,j) = FSEQ_1yr(i,j,TOPO(i,j));
    end
  end
end

% make a 3-d array of the fraction sequestered after 5 years
[ny,nx,nz] = size(MASK);
FSEQ_5yr = 0*MASK+NaN;
t_indx = find(time==5); % find the index of the 5'th year
FSEQ_5yr(MASK==1) = fseq(:,t_indx);
% now find the near-bottom ocean grid cells
fseq_bottom_5yr = 0*MASK(:,:,1)+NaN;
TOPO = sum(MASK,3); % number of grid cells in the vertical direction
for i = 1:ny
  for j = 1:nx
    if TOPO(i,j)~=0
      fseq_bottom_5yr(i,j) = FSEQ_5yr(i,j,TOPO(i,j));
    end
  end
end

% make a 3-d array of the fraction sequestered after 10 years
[ny,nx,nz] = size(MASK);
FSEQ_10yr = 0*MASK+NaN;
t_indx = find(time==10); % find the index of the 10'th year
FSEQ_10yr(MASK==1) = fseq(:,t_indx);
% now find the near-bottom ocean grid cells
fseq_bottom_10yr = 0*MASK(:,:,1)+NaN;
TOPO = sum(MASK,3); % number of grid cells in the vertical direction
for i = 1:ny
  for j = 1:nx
    if TOPO(i,j)~=0
      fseq_bottom_10yr(i,j) = FSEQ_10yr(i,j,TOPO(i,j));
    end
  end
end

% make a 3-d array of the fraction sequestered after 25 years
[ny,nx,nz] = size(MASK);
FSEQ_25yr = 0*MASK+NaN;
t_indx = find(time==25); % find the index of the 25'th year
FSEQ_25yr(MASK==1) = fseq(:,t_indx);
% now find the near-bottom ocean grid cells
fseq_bottom_25yr = 0*MASK(:,:,1)+NaN;
TOPO = sum(MASK,3); % number of grid cells in the vertical direction
for i = 1:ny
  for j = 1:nx
    if TOPO(i,j)~=0
      fseq_bottom_25yr(i,j) = FSEQ_25yr(i,j,TOPO(i,j));
    end
  end
end

% make a 3-d array of the fraction sequestered after 75 years
[ny,nx,nz] = size(MASK);
FSEQ_75yr = 0*MASK+NaN;
t_indx = find(time==75); % find the index of the 75'th year
FSEQ_75yr(MASK==1) = fseq(:,t_indx);
% now find the near-bottom ocean grid cells
fseq_bottom_75yr = 0*MASK(:,:,1)+NaN;
TOPO = sum(MASK,3); % number of grid cells in the vertical direction
for i = 1:ny
  for j = 1:nx
    if TOPO(i,j)~=0
      fseq_bottom_75yr(i,j) = FSEQ_75yr(i,j,TOPO(i,j));
    end
  end
end

% make a 3-d array of the fraction sequestered after 100 years
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

% make a 3-d array of the fraction sequestered after 1000 years
[ny,nx,nz] = size(MASK);
FSEQ_1000yr = 0*MASK+NaN;
t_indx = find(time==1000); % find the index of the 1000'th year
FSEQ_1000yr(MASK==1) = fseq(:,t_indx);
% now find the near-bottom ocean grid cells
fseq_bottom_1000yr = 0*MASK(:,:,1)+NaN;
TOPO = sum(MASK,3); % number of grid cells in the vertical direction
for i = 1:ny
  for j = 1:nx
    if TOPO(i,j)~=0
      fseq_bottom_1000yr(i,j) = FSEQ_1000yr(i,j,TOPO(i,j));
    end
  end
end

% calculate the sum of the emission fractions (1 - sequestration fraction)
% for every year from 1 to 100 years; store
SUM_FEMIT_1to100 = zeros(100,2);
SUM_FEMIT_1to100(:,1) = 1:100;

for i = 1:100
    [ny,nx,nz] = size(MASK);
    FSEQ_thisyr = 0*MASK+NaN;
    t_indx = find(time==i); % find the index of the i'th year
    FSEQ_thisyr(MASK==1) = fseq(:,t_indx);
    % now find the near-bottom ocean grid cells
    fseq_bottom_thisyr = 0*MASK(:,:,1)+NaN;
    TOPO = sum(MASK,3); % number of grid cells in the vertical direction
    for j = 1:ny
        for k = 1:nx
            if TOPO(j,k)~=0
                fseq_bottom_thisyr(j,k) = FSEQ_thisyr(j,k,TOPO(j,k));
            end
        end
    end
    SUM_FEMIT_thisyr = sum((1-fseq_bottom_thisyr()), 'all', 'omitnan');
    SUM_FEMIT_1to100(i,2) = SUM_FEMIT_thisyr;
end

% make a 3-d array of the fractions sequestered at the benthic depth
% after 1-200 (inclusive) years, plus years 300, 400, ... 900, 1000
fseq_bottom_multyears = zeros([size(MASK,1:2), 208]);
years = [1:200 300 400 500 600 700 800 900 1000];

for i = 1:length(years)
    [ny,nx,nz] = size(MASK);
    FSEQ_thisyr = 0*MASK+NaN;
    t_indx = find(time==years(i)); % find the index of the i'th year
    FSEQ_thisyr(MASK==1) = fseq(:,t_indx);
    % now find the near-bottom ocean grid cells
    fseq_bottom_thisyr = 0*MASK(:,:,1)+NaN;
    TOPO = sum(MASK,3); % number of grid cells in the vertical direction
    for j = 1:ny
        for k = 1:nx
            if TOPO(j,k)~=0
                fseq_bottom_thisyr(j,k) = FSEQ_thisyr(j,k,TOPO(j,k));
            end
        end
    end
    (:,:,i) = fseq_bottom_thisyr;
end

% find the bottom depth
bottom_depth = sum(VOL.*MASK,3)./AREA(:,:,1); % depth of the ocean at each water column

% export the individual year arrays as .csv; multiple year array of
% matrixes .mat file

writematrix(fseq_bottom_1yr,'/Users/jamesrco/Code/zoonotic-C/data/global_trawling/derived/benthic_seqfractions/fseq_bottom_1yr.csv')
writematrix(fseq_bottom_5yr,'/Users/jamesrco/Code/zoonotic-C/data/global_trawling/derived/benthic_seqfractions/fseq_bottom_5yr.csv')
writematrix(fseq_bottom_10yr,'/Users/jamesrco/Code/zoonotic-C/data/global_trawling/derived/benthic_seqfractions/fseq_bottom_10yr.csv')
writematrix(fseq_bottom_25yr,'/Users/jamesrco/Code/zoonotic-C/data/global_trawling/derived/benthic_seqfractions/fseq_bottom_25yr.csv')
writematrix(fseq_bottom_50yr,'/Users/jamesrco/Code/zoonotic-C/data/global_trawling/derived/benthic_seqfractions/fseq_bottom_50yr.csv')
writematrix(fseq_bottom_75yr,'/Users/jamesrco/Code/zoonotic-C/data/global_trawling/derived/benthic_seqfractions/fseq_bottom_75yr.csv')
writematrix(fseq_bottom_100yr,'/Users/jamesrco/Code/zoonotic-C/data/global_trawling/derived/benthic_seqfractions/fseq_bottom_100yr.csv')
writematrix(fseq_bottom_1000yr,'/Users/jamesrco/Code/zoonotic-C/data/global_trawling/derived/benthic_seqfractions/fseq_bottom_1000yr.csv')
writematrix(bottom_depth,'/Users/jamesrco/Code/zoonotic-C/data/global_trawling/derived/benthic_seqfractions/bottom_depth_m.csv')
writematrix(SUM_FEMIT_1to100,'/Users/jamesrco/Code/zoonotic-C/data/global_trawling/derived/benthic_seqfractions/sum_of_fseq_bottom_1to100years.csv')
writematrix(LAT(:,1,1),'/Users/jamesrco/Code/zoonotic-C/data/global_trawling/derived/benthic_seqfractions/lat_degN.csv')
writematrix(LON(1,:,1),'/Users/jamesrco/Code/zoonotic-C/data/global_trawling/derived/benthic_seqfractions/long_degE.csv')

writematrix(years,'/Users/jamesrco/Code/zoonotic-C/data/global_trawling/derived/benthic_seqfractions/benthic_years.csv')
save('/Users/jamesrco/Code/zoonotic-C/data/global_trawling/derived/benthic_seqfractions/fseq_bottom_multyears.mat','fseq_bottom_multyears')

% % plot both the fraction remaining after 100 years and the bottom depth
% figure(1)
% subplot(2,1,1)
% pcolor(LON(:,:,1),LAT(:,:,1),fseq_bottom_100yr)
% set(gca,'clim',[0 1])
% colormap(parula(10))
% colorbar
% set(gca,'color',[.8 .8 .8])
% set(gcf,'color','w')
% shading('flat')
% title('Fraction of CO2 injected at seafloor remaining after 100 years')
% subplot(2,1,2)
% pcolor(LON(:,:,1),LAT(:,:,1),bottom_depth)
% set(gca,'clim',[0 5500])
% colormap(parula(11))
% colorbar
% set(gca,'color',[.8 .8 .8])
% set(gcf,'color','w')
% shading('flat')
% title('Depth of seafloor (m)')

% % Example 2: Find a particular point in the ocean and plot the fraction
% % remaining sequestered over time
% % first specify a location by latitude, longitude, and depth
% lat = 40; % degrees north
% lon = 220; % degrees east, from 0 - 360
% depth = 1000; % meters
% % the following code will find the nearest model grid cell(s)
% iy = find(abs(LAT(:)-lat)==min(abs(LAT(:)-lat)));
% ix = find(abs(LON(:)-lon)==min(abs(LON(:)-lon)));
% iz = find(abs(DEPTH(:)-depth)==min(abs(DEPTH(:)-depth)));
% indx = intersect(intersect(ix,iy),iz);
% % check if that grid cell(s) is land or ocean; if in ocean, plot fseq
% if sum(MASK(indx))==0
%   fprintf('Location is not in the ocean. Please choose another location.\n')
% else
%   % if in ocean, find location in ocean grid coordinates and plot
%   iy2 = find(abs(LAT(MASK==1)-lat)==min(abs(LAT(MASK==1)-lat)));
%   ix2 = find(abs(LON(MASK==1)-lon)==min(abs(LON(MASK==1)-lon)));
%   iz2 = find(abs(DEPTH(MASK==1)-depth)==min(abs(DEPTH(MASK==1)-depth)));
%   indx2 = intersect(intersect(ix2,iy2),iz2);
%   figure(2)
%   plot(time,mean(fseq(indx2,:)))
%   xlabel('Time since injection (years)')
%   ylabel('Fraction of CO2 remaining sequestered')
%   title('Sequestration fraction over time')
%   set(gca,'ylim',[0 1])
% end
