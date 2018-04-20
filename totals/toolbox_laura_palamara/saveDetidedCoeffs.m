function saveDetidedCoeffs(outfile,starttime,endtime,period,varargin)

% Inputs:
%     outfile: string, name of netcdf file to create with detiding
%       coefficients
%     starttime: earliest time to use for detiding (matlab time units)
%     endtime: latest time to use for detiding (matlab time units)
%     period: 1-D array of period values for major tidal components
% 
% Options:
%     mincoverage: proportion (0-1) of time range that must have hourly u 
%       and v coverage in order to be detided (default: .25)
%     timeformat: format of time in file name of surface current total
%       netcdfs (default: 'yyyy_mm_dd_HHMM')
%     periodunits: units of period input (default: 'hours', other
%       options: 'days','minutes','seconds')
%     maxerror: maximum error allowed for including surface current data 
%       (default: .6)
%     inputdir: string of directory with input surface current netcdf 
%       files (default: current working directory)
%     lonbounds: min and max longitude to read in (formatted [minlon
%       maxlon]) (default: read in all)
%     latbounds: min and max latitude to read in (formatted [minlat
%       maxlat]) (default: read in all)
%     overwrite: logical value indicating whether to replace data in netcdf
%       file (default: false)
%

caller = [mfilename '.m'];

if nargin <4
    fprintf(2,...
        '%s E: not enough arguments.\n',...
        caller);
    return;
elseif ~isnumeric(starttime) | ~isnumeric(endtime) | ~isnumeric(period)
    fprintf(2,...
        '%s E: starttime,endtime, and period must all be numeric.\n',...
        caller);
    return;
elseif length(starttime)~=1 | length(endtime)~=1
    fprintf(2,'%s E: starttime and endtime must be numeric scalars.\n',caller);
    return;
elseif endtime<=starttime
    fprintf(2,'%s E: endtime must be later than starttime.\n',caller);
    return;
elseif ~ischar(outfile)
    fprintf(2,'%s E: outfile must be a string.\n',caller);
    return;
end


inputdir=[];
filetimeformat='yyyy_mm_dd_HHMM';
maxerr=.6;
periodunits='hours';
mincoverage=.25;
lonbounds=[-360 360];
latbounds=[-90 90];
replacevals=false;

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
        case 'timeformat'
            if(~ischar(value))
                fprintf(2,'%s E: Value for %s must be a string\n',caller,name);
                return;
            end
            filetimeformat=value;
        case 'periodunits'
            if(~ischar(value))
                fprintf(2,'%s E: Value for %s must be a string\n',caller,name);
                return;
            end
            periodunits=value;
        case 'inputdir'
            if(~ischar(value))
                fprintf(2,'%s E: Value for %s must be a a string\n',caller,name);
                return;
            end
            if(~exist(value,'dir'))
                fprintf(2,'%s E: Input directory %s not found\n',caller,value);
                return;
            end
            if(value(end)~='/')
                value=[value '/'];
            end
            inputdir=value;
        case 'maxerror'
            if(~isnumeric(value))
                fprintf(2,'%s E: Value for %s must be a numeric scalar\n',caller,name);
                return;
            end
            if(length(value)~=1)
                fprintf(2,'%s E: Value for %s must be a numeric scalar\n',caller,name);
                return;
            end
            maxerr=value;
        case 'lonbounds'
            if(~isnumeric(value))
                fprintf(2,'%s E: Value for %s must be numeric\n',caller,name);
                return;
            end
            if(length(value)~=2)
                fprintf(2,'%s E: Value for %s must be 2 elements (min and max)\n',caller,name);
                return;
            end
            lonbounds=value;
        case 'latbounds'
            if(~isnumeric(value))
                fprintf(2,'%s E: Value for %s must be numeric\n',caller,name);
                return;
            end
            if(length(value)~=2)
                fprintf(2,'%s E: Value for %s must be 2 elements (min and max)\n',caller,name);
                return;
            end
            latbounds=value;
        case 'overwrite'
            if(~islogical(value))
                fprintf(2,'%s E: Value for %s must be logical\n',caller,name);
                return;
            end
            replacevals=value;
        otherwise
            fprintf(2,...
                '%s E: Unknown option specified: %s\n',...
                caller,...
                name);
            return;
    end
end


timeformatnums=datestr(datenum(1111,11,11,11,11,11),filetimeformat);
inputname=dir([inputdir '*.nc']);
if(isempty(inputname))
    fprintf(2,'%s E: No netcdf files in %s\n',caller,inputdir);
end
inputname=inputname(1).name;
inputnamenums=inputname;
inputnamenums(inputname>='0'&inputname<='9')='1';
timestart=strfind(inputnamenums,timeformatnums);
if(length(timestart)~=1)
    fprintf(2,'%s E: netcdf files do not match time format %s\n',caller,inputname,filetimeformat);
    return;
end
fileformatstart=inputname(1:timestart-1);
fileformatend=inputname(timestart+length(filetimeformat):end);

allfiles=dir([inputdir fileformatstart '*' fileformatend]);
allfiles=struct2cell(allfiles);
allfiles=char(allfiles(1,:));
indf=find(datenum(allfiles(:,timestart:timestart+length(filetimeformat)-1),filetimeformat)>=starttime&datenum(allfiles(:,timestart:timestart+length(filetimeformat)-1),filetimeformat)<=endtime);

if(exist(outfile,'file')&~replacevals)
    fprintf(2,'%s E: %s already exists. to replace, set ''overwrite'' to true\n',caller,outfile);
    return;
end
if(exist(outfile,'file'))
    delete(outfile);
end

ncid=netcdf.open([inputdir allfiles(indf(1),:)],'NC_NOWRITE');
tID=netcdf.inqVarID(ncid,'time');
tunits1=netcdf.getAtt(ncid,tID,'units');
lonID=netcdf.inqVarID(ncid,'lon');
latID=netcdf.inqVarID(ncid,'lat');
uID=netcdf.inqVarID(ncid,'u');
ufill=netcdf.getAtt(ncid,uID,'_FillValue');
vID=netcdf.inqVarID(ncid,'v');
vfill=netcdf.getAtt(ncid,vID,'_FillValue');
uerrID=netcdf.inqVarID(ncid,'u_err');
uerrfill=netcdf.getAtt(ncid,uerrID,'_FillValue');
verrID=netcdf.inqVarID(ncid,'v_err');
verrfill=netcdf.getAtt(ncid,verrID,'_FillValue');
lon=netcdf.getVar(ncid,lonID);
lat=netcdf.getVar(ncid,latID);
netcdf.close(ncid)

timeoptions={'days','hours','minutes','seconds'};
indt=find(strncmpi(timeoptions,tunits1,4));
tunits=timeoptions{indt};


outFile = [inputdir 'detidingcoefficients/Archive/' outfile];
coeffsncid=netcdf.create(outFile,'NETCDF4');
coefflondimID=netcdf.defDim(coeffsncid,'lon',length(lon));
coefflatdimID=netcdf.defDim(coeffsncid,'lat',length(lat));
coeffdimID=netcdf.defDim(coeffsncid,'soln',length(period)*2+1);
coeffperdimID=netcdf.defDim(coeffsncid,'period',length(period));
coefflonID=netcdf.defVar(coeffsncid,'lon','float',coefflondimID);
netcdf.putAtt(coeffsncid,coefflonID,'units','degrees_east');
coefflatID=netcdf.defVar(coeffsncid,'lat','float',coefflatdimID);
netcdf.putAtt(coeffsncid,coefflatID,'units','degrees_north');
coeffuID=netcdf.defVar(coeffsncid,'solnu','double',[coeffdimID coefflondimID coefflatdimID]);
coeffvID=netcdf.defVar(coeffsncid,'solnv','double',[coeffdimID coefflondimID coefflatdimID]);
netcdf.putAtt(coeffsncid,coeffuID,'description',['least-squares coefficients for u velocity, calculated ' datestr(starttime) ' through ' datestr(endtime) ' with time in ' tunits1]);
netcdf.putAtt(coeffsncid,coeffvID,'description',['least-squares coefficients for v velocity, calculated ' datestr(starttime) ' through ' datestr(endtime) ' with time in ' tunits1]);
netcdf.putAtt(coeffsncid,coeffuID,'missing_value',-999);
netcdf.defVarFill(coeffsncid,coeffuID,false,-999);
netcdf.putAtt(coeffsncid,coeffvID,'missing_value',-999);
netcdf.defVarFill(coeffsncid,coeffvID,false,-999);
coeffnpID=netcdf.defVar(coeffsncid,'np','float',[coefflondimID coefflatdimID]);
netcdf.putAtt(coeffsncid,coeffnpID,'description','number of points that went in to get detiding coefficients');
netcdf.putAtt(coeffsncid,coeffnpID,'totalpointspossible',length(indf))
netcdf.putAtt(coeffsncid,-1,'date_range',[datestr(starttime,'mmm dd yyyy HH:MM') ' to ' datestr(endtime,'mmm dd yyyy HH:MM')])
netcdf.defVarFill(coeffsncid,coeffnpID,false,single(-999));
coeffperID=netcdf.defVar(coeffsncid,'period','double',coeffperdimID);
netcdf.putAtt(coeffsncid,coeffperID,'description','periods of major tidal components included in detiding');
netcdf.putAtt(coeffsncid,coeffperID,'units',periodunits);
netcdf.endDef(coeffsncid);
netcdf.putVar(coeffsncid,coefflonID,lon);
netcdf.putVar(coeffsncid,coefflatID,lat);
netcdf.putVar(coeffsncid,coeffperID,period);

indlon=find(lon>=min(lonbounds)&lon<=max(lonbounds));
if(isempty(indlon))
    indlon=find(lon>=min(lonbounds+360)&lon<=max(lonbounds+360));
    if(isempty(indlon))
        indlon=find(lon>=min(lonbounds-360)&lon<=max(lonbounds-360));
    end
end
indlat=find(lat>=min(latbounds)&lat<=max(latbounds));
if(isempty(indlon)|isempty(indlat))
    netcdf.close(coeffsncid)
    fprintf(2,'%s E: no gridpoints found between latitudes %f and %f and longitudes %f and %f\n',caller,min(latbounds),max(latbounds),min(lonbounds),max(lonbounds));
    return;
end

lon=lon(indlon);
lat=lat(indlat);


for i=0:length(lon)-1
    time=nan(length(indf),1);
    u=nan(length(indf),length(lat));
    v=u;
    for ti=1:length(indf)
        
        ncid=netcdf.open([inputdir allfiles(indf(ti),:)],'NC_NOWRITE');
        codartime=netcdf.getVar(ncid,tID);
        if(~isempty(codartime))
            time(ti)=codartime;
            uerr=netcdf.getVar(ncid,uerrID,[i+min(indlon)-1 min(indlat)-1 0],[1 length(lat) 1]);
            verr=netcdf.getVar(ncid,verrID,[i+min(indlon)-1 min(indlat)-1 0],[1 length(lat) 1]);
            u(ti,:)=double(netcdf.getVar(ncid,uID,[i+min(indlon)-1 min(indlat)-1 0],[1 length(lat) 1]));
            v(ti,:)=double(netcdf.getVar(ncid,vID,[i+min(indlon)-1 min(indlat)-1 0],[1 length(lat) 1]));
            ind=find(abs(uerr)>maxerr|abs(verr)>maxerr);
            u(ti,ind)=NaN;
            v(ti,ind)=NaN;
        end
        netcdf.close(ncid)
    end
    u(u==ufill)=NaN;
    v(v==vfill)=NaN;
    for j=0:length(lat)-1
        [solnu,solnv]=getDetidedCoeffs(time,u(:,j+1),v(:,j+1),period,'periodunits',periodunits,'timeunits',tunits,'mincoverage',mincoverage);
        
        ind=find(~isnan(u(:,j+1)));
        netcdf.putVar(coeffsncid,coeffnpID,[i+min(indlon)-1 j+min(indlat)-1],[1 1],length(ind));
        ind=find(isnan(solnu));
        solnu(ind)=-999;
        ind=find(isnan(solnv));
        solnv(ind)=-999;
        netcdf.putVar(coeffsncid,coeffuID,[0 i+min(indlon)-1 j+min(indlat)-1],[length(period)*2+1 1 1],solnu);
        netcdf.putVar(coeffsncid,coeffvID,[0 i+min(indlon)-1 j+min(indlat)-1],[length(period)*2+1 1 1],solnv);
        clear ind solnu solnv
    end
end


netcdf.sync(coeffsncid)
netcdf.close(coeffsncid)

copyfile(outFile, [inputdir 'detidingcoefficients/detidingcoefficients_current.nc'])

