function [xr,yr]=rot(x,y,theta,ya)

% ROT rotates a vector counterclockwise,
%
% [xr,yr]=ROT(x,y,theta,ya) rotates a vector counterclockwise
% theta degrees OR rotates the coordinate system clockwise
% theta degrees. Example: rot(1,0,90) returns (0,1).
%
% If 4 arguments, then it uses xa and ya as orientation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ver. 1: 11/17/96 (RG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==4
    theta=-atan2(ya,theta)*180/pi;
end

costheta=cos(theta/180*pi);
sintheta=sin(theta/180*pi);

if length(theta)==1
    xr=x*costheta-y*sintheta;
    yr=x*sintheta+y*costheta;
else
    xr=x.*costheta-y.*sintheta;
    yr=x.*sintheta+y.*costheta;
end