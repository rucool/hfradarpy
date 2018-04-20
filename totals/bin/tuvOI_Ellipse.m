
function [u,v,xi] = tuvOI_Ellipse(rad_Speed, angle, rloc, xyloc, mdlvar, errvar, sx, sy, orient, option)
%
% TUVOI  Calculate total current vectors from radial currents using optimal interpolation method
%
% tuvOI  - calculates the u/v components of a total vector from 2 to n
%   radial vector components in a given model variance and error variance
%
% Usage:
%           [u,v,xi] = tuvOI(rad_Speed, angle, rloc, xyloc, mdlvar, errvar, sx, sy, option)
%
% Inputs:  
%   rad_Speed - radial current speed
%
%   rad_Speed - radial velocity falling into the serach radius of
%               each vector grid point.
%
%   angle     - bearing angle, in degrees
%
%   rloc      - radial grid locations within a search radius 
%               (:,1) = longitude (:,2) = latitude
%
%   xyloc     - grid location where uv vector components are estimated
%
%   errvar    - observational error variance, which can be constant 
%               or hourly standard deviation of radial velocity (HSTD, 
%               so called temporal uncertainty). For example, if the 
%               observational uncertainty in HF radar observations is 3-5 cm/s,
%               the corresponding error variance = 9-25 (cm/s)^2. 
%				As a standard error, HSTD^2/N can be used as the error variance.
%				N is the number of cross spectra within a given time period.
%				Although N is unknown, the bound of N is known. 1<= N <= 6.
%
%   mdlvar    - a prior model covariance of the surface currents.
%               In the pointwise approach, user can set the a prior model 
%               variance as a funcion of water depth, the length from the 
%               coast line, or constant.
%					
%   sx, sy    - decorrelation legnth scale defining the elipse along the
%               major (sy) and minor axis (sx) Added by Josh
%
%   orient    - orientation of the averaging elipse (degrees - positive 
%               clockwise) Added by Josh
%   option    - for correlation function; option == 1,
%               Gaussian, option == 2, 
%               exponential function
%
% Outputs:
%         u,v - the total currents (2 by 1 matrix)
%         xi - uncertainty normalized by the a prior model covariance. (2 by 2 matrix)
%			xi(1,1) : normalized uncertainty of u = <(u_hat - u)^2>/<u^2> (good :0, poor: 1)
%			xi(2,2) : normalized uncertainty of v = <(v_hat - v)^2>/<v^2> (good :0, poor: 1)
%			xi(1,2) : directional information of u and v = <(u_hat -u)(v_hat- v)>/sqrt(<u^2><v^2>)
%

%  Make sure that angle and rad_Speed are vectors
if ~isvector(rad_Speed)  ||  ~isvector(angle)
   error('%s: angle and radial velocity inputs *MUST* be vectors',mfilename);
end

if ~exist( 'raderr', 'var' ), raderr = []; end

% Make sure the vectors are column vectors.
angle = angle(:);
rad_Speed = rad_Speed(:);

% Form the angle (A) and radial_matrices
d2r = pi/180;
angle = angle*d2r;

nr = length(angle); 
if length(errvar) > 1.  R = diag(errvar); %for HSTD (socalled hourly temporal uncertainty)
else, R = eye(nr)*errvar; %%for constant errvariance
end
P_ = eye(2)*mdlvar;


[dx, dy] = lonlat2km(rloc(:,1), rloc(:,2), xyloc(:,1), xyloc(:,2));
[x1, x2] = meshgrid(rloc(:,1), rloc(:,1)); [y1, y2] = meshgrid(rloc(:,2), rloc(:,2)); 
[ang1, ang2] = meshgrid(angle, angle);
[drx_, dry_] = lonlat2km(x1, y1, x2, y2); 

% Add a rotation step given the orientation of the ellipse
% this step will rotate the distance and location to a coordinate system
% that matcheds the desired ellipse orientation
[dx_rot, dy_rot]=rot(dx, dy, orient);
[drx_rot, dry_rot]=rot(drx_, dry_, orient);

if option == 1,
	cmd = exp(-sqrt(dx_rot.^2/sx^2 + dy_rot.^2/sy^2));
	w_ = exp(-sqrt(drx_rot.^2/sx^2 + dry_rot.^2/sy^2))*mdlvar;
end

if option == 2,
	cmd = exp(-sqrt(dx_rot.^2/sx^2 + dy_rot.^2/sy^2));
	w_ = exp(-sqrt(drx_rot.^2/sx^2 + dry_rot.^2/sy^2))*mdlvar;
end

cmdu = cmd.*cos(angle); cmdv = cmd.*sin(angle); clear cmd
cmd(:,1) = cmdu; cmd(:,2) = cmdv;
cdd = w_.*(cos(ang1).*cos(ang2) + sin(ang1).*sin(ang2)); clear x1 x2 ang1 ang2 ang x y
cmdicdd = P_*cmd'*inv(cdd +R);  a = cmdicdd*rad_Speed; 
xi = inv(P_)*(P_ - cmdicdd*cmd*P_);



u = a(1);
v = a(2);