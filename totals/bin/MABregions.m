function [region,rotang]=MABregions(lon,lat)


% Region 1: north of 41.66, v is along, u is cross
% Region 2: south of 41.66 and east of -69.95, rotation angle 36.88 (that
% was based on the old grid I was using for 03-07 - may not be the best
% number but it worked well enough for us)
% Region 3: west of -69.95, north of the line lat = -.45*lon+7.1725 (a
% couple km north of the canyon), u is along, -v is cross
% Region 4: south of that line, same rotation angle as Region 2

% Regions go north (1) to south (4)
% Rotang is the angle clockwise from straight north-south
% We haven't looked at anything north of 44 degrees
% 37 degrees and south looks like it should be a new defined region (maybe 2)



region=NaN*ones(size(lon));
rotang=NaN*ones(size(lon));


ind=find(lat>=41.66);
region(ind)=1;
rotang(ind)=0;

ind=find(lon>=-69.95&isnan(region));
region(ind)=2;
rotang(ind)=36.88;

latline=-.45*lon+7.1725;
ind=find(lat>=latline&isnan(region));
region(ind)=3;
rotang(ind)=90;

ind=find(lat<latline&isnan(region));
region(ind)=4;
rotang(ind)=36.88;