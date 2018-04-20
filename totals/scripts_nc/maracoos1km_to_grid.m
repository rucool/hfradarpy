function [P] =  maracoos1km_to_grid(grid_file)

% codar_grid_dir = '/home/codaradm/data/grid_files/';%'/Volumes/lemus/MARACOOS/grid_hfrnet/';%
% codar_grid_file = '1km_sr-grid.txt';

try
  d = load(grid_file);
catch
  codar_grid_dir = [];
  d = load(grid_file);
end
  
% split out the lon/lat columns of the grid file
x = d(:,1);
y = d(:,2);

% compute dx and dy (lon/lat grid increments) from all possible values in
% the grid file
dx = unique(diff(sort(x)));
dx = dx(2);
dy = unique(diff(sort(y)));
dy = dy(2);

% calculate number of lon/lat coordinates when put on a 2D grid
nx = 1+round(diff([min(x) max(x)])/dx);
ny = 1+round(diff([min(y) max(y)])/dy);

% make vectors of grid lon/lat values that span the entire range
lon = linspace(min(x),max(x),nx);
lat = linspace(min(y),max(y),ny);

% expand into a 2-D set of lon,lat coordinates
[X,Y] = meshgrid(lon,lat);

% find the 1-D index into coordinates corresponding to each row of the
% codar grid file
P = NaN*ones(size(x));
for k=1:length(x)
  del=[(X(:)-x(k)) + sqrt(-1)*(Y(:)-y(k))];
  [ignore,P(k)]=min(del); % don't need the min, only it's location P(k)
end

% example: dummy data are the pointwise lon values
data = x;
DATA = NaN*ones(size(X));
DATA(P) = data;
end
