function [solnu,solnv]=getDetidedCoeffs(time,u,v,period,varargin)

% Get detiding coefficients for a time-series of u and v currents.
%
% Outputs:
%     solnu: detiding coefficients for u (1st element is mean, following
%       elements alternate a and b coefficients for corresponding period)
%     solnv: detiding coefficients for v (1st element is mean, following
%       elements alternate a and b coefficients for corresponding period)
%
% Inputs:
%     time: 1-D array of time values (in matlab time or any unit increasing
%       by day
%     u: 1-D array of eastward velocity corresponding to times
%     v: 1-D array of northward velocity corresponding to times
%     period: 1-D array of period values for major tidal components
% 
% Options:
%     mincoverage: proportion (0-1) of time range that must have hourly u 
%       and v coverage in order to be detided (default: .25)
%     timeunits: unit that time increases in (default: 'days', other 
%       options: 'hours','minutes','seconds')
%     periodunits: units of period input (default: 'hours', other
%       options: 'days','minutes','seconds')
%

solnu=[];
solnv=[];

caller = [mfilename '.m'];

if nargin <4
    fprintf(2,...
        '%s E: not enough arguments.\n',...
        caller);
    return;
elseif ~isnumeric(time) | ~isnumeric(u) | ~isnumeric(v) | ~isnumeric(period)
    fprintf(2,...
        '%s E: time, u, v, and period must all be numeric.\n',...
        caller);
    return;
elseif size(time)~=size(u) | size(time)~=size(v)
    fprintf(2,'%s E: time, u, and v must have the same dimensions.\n',caller);
    return;
elseif length(size(time))~=2
    fprintf(2,'%s E: time, u, and v must be one-dimensional arrays.\n',caller);
    return;
end

if(size(time,2)~=1)
    time=time';
    u=u';
    v=v';
end


for x=1:2:length(varargin)
    name=varargin{x};
    value=varargin{x+1};
    switch lower(name)
        case 'mincoverage'
            if(~isnumeric(value))
                fprintf(2,'%s E: Value for %s must be a numeric\n',caller,name);
                return;
            end
            if(length(value)>1|value<0|value>1)
                fprintf(2,'%s E: Value for %s must be between 0 and 1\n',caller,name);
                return;
            end
            mincoverage=value;
        case 'timeunits'
            if(~ischar(value))
                fprintf(2,'%s E: Value for %s must be a string\n',caller,name);
                return;
            end
            switch lower(value)
                case 'days'
                case 'hours'
                    period=period*24;
                case 'minutes'
                    period=period*24*60;
                case 'seconds'
                    period=period*24*60*60;
                otherwise
                    fprintf(2,'%s E: Unknown option for %s: %s\n',caller,name,value);
                    return;
            end
        case 'periodunits'
            if(~ischar(value))
                fprintf(2,'%s E: Value for %s must be a string\n',caller,name);
                return;
            end
            switch lower(value)
                case 'days'
                    period=period*24;
                case 'hours'
                case 'minutes'
                    period=period/60;
                case 'seconds'
                    period=period/60/60;
                otherwise
                    fprintf(2,'%s E: Unknown option for %s: %s\n',caller,name,value);
                    return;
            end
        otherwise
            fprintf(2,...
                '%s E: Unknown option specified: %s\n',...
                caller,...
                name);
            return;
    end
end


mincoverage=.25;

omega=(2*24*pi)./period;

ind=find(~isnan(u));
totalnp=length(time);
range=time(end)-time(1)+1/24;

if(length(ind)>=mincoverage*totalnp)
    solnu=leastsquares(time(ind), omega, length(period), u(ind), length(ind));
    solnv=leastsquares(time(ind), omega, length(period), v(ind), length(ind));
else
    solnu=NaN*ones(length(period)*2+1,1);
    solnv=NaN*ones(length(period)*2+1,1);
    
end

