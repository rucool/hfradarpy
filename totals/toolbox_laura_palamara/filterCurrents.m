function [lon,lat,originalu,originalv,lowpassu,lowpassv,highpassu,highpassv]=filterCurrents(inputfile,varargin)

% Get tidal and detided current velocities.
%
% Outputs:
%     lon: longitude vector for grid
%     lat: latitude vector for grid
%     originalu: matrix of eastward velocity used in filtering (from inputfile)
%     originalv: matrix of northward velocity used in filtering (from inputfile)
%     lowpassu: matrix of low-pass filtered eastward velocity
%     lowpassv: matrix of low-pass filtered northward velocity
%     highpassu: matrix of high-pass filtered eastward velocity
%     highpassv: matrix of high-pass filtered northward velocity
%
% Inputs:
%     inputfile: file name of netcdf that includes the currents to filter
%
% Options:
%     starttime: earliest time to use in the filtering (default: 1 week
%       before inputfile)
%     endtime: latest time to use in the filtering (default: 1 week after
%       inputfile)
%     timestep: interval between each time to read in data from, in hours
%       (default: 1)
%     timeformat: format of time in file name of surface current total
%       netcdfs (default: 'yyyy_mm_dd_HHMM')
%     currenttype: variable name to read in and filter, formatted with 'X'
%       replacing 'u' or 'v' (i.e. to read in 'u_detided' and 'v_detided',
%       option would be 'X_detided'; default: 'X_detided')
%     filterlength: time length for low-pass filter, in hours (default: 30)
%     tail: 2-element vector indicating how much data needs to be available
%       on the before and after tail end of the file for the filtered data
%       to be considered "good" (default: [3 12] indicates that 3 out of
%       the 12 elements immediately before and 3 out of the 12 elements
%       immediately after the inputfile must have data
%     addlowpass: logical value indicating whether to add low-pass
%       filtered velocity to the original surface current netcdf file
%       (inputfile) (default: false)
%     addhighpass: logical value indicating whether to add high-pass
%       filtered velocity to the original surface current netcdf file
%       (inputfile) (default: false)
%     overwrite: logical value indicating whether to replace low-pass
%       and high-pass filtered velocity if they already exist in the
%       netcdf (addlowpass and/or addhighpass must be set to true;
%       default: false)
%     lonbounds: min and max longitude to read in (formatted [minlon
%       maxlon]) (default: read in all)
%     latbounds: min and max latitude to read in (formatted [minlat
%       maxlat]) (default: read in all)
%

lon=[];
lat=[];
originalu=[];
originalv=[];
lowpassu=[];
lowpassv=[];
highpassu=[];
highpassv=[];

caller = [mfilename '.m'];

if nargin <1
    fprintf(2,...
        '%s E: not enough arguments.\n',...
        caller);
    return;
elseif ~ischar(inputfile)
    fprintf(2,...
        '%s E: inputfile must be a string.\n',...
        caller);
    return;
elseif ~exist(inputfile,'file')
    fprintf(2,'%s E: inputfile %s not found.\n',caller,inputfile);
    return;
end


filetimeformat='yyyy_mm_dd_HHMM';
inputvariable='X_detided';
tail=[3 12];
timestep=1/24;
filttime=30;
addlowpass=false;
addhighpass=false;
replacevals=false;
lonbounds=[-360 360];
latbounds=[-90 90];

for x=1:2:length(varargin)
    name=varargin{x};
    value=varargin{x+1};
    switch lower(name)
        case 'timeformat'
            if(~ischar(value))
                fprintf(2,'%s E: Value for %s must be a string\n',caller,name);
                return;
            end
            filetimeformat=value;
    end
end


[inputdir,inputname,inputext]=fileparts(inputfile);
inputdir=[inputdir '/'];

timeformatnums=datestr(datenum(1111,11,11,11,11,11),filetimeformat);
inputnamenums=inputname;
inputnamenums(inputname>='0'&inputname<='9')='1';
timestart=strfind(inputnamenums,timeformatnums);
if(length(timestart)~=1)
    fprintf(2,'%s E: file %s does not match time format %s\n',caller,inputfile,filetimeformat);
    return;
end
fileformatstart=inputname(1:timestart-1);
fileformatend=[inputname(timestart+length(filetimeformat):end) inputext];

currenttime=datenum(inputname(timestart:timestart+length(filetimeformat)-1),filetimeformat);
starttime=currenttime-7;
endtime=currenttime+7;

inputinfo=ncinfo(inputfile);


for x=1:2:length(varargin)
    name=varargin{x};
    value=varargin{x+1};
    switch lower(name)
        case 'starttime'
            if(~isnumeric(value)|length(value)~=1)
                fprintf(2,'%s E: Value for %s must be a numeric scalar\n',caller,name);
                return;
            end
        case 'endtime'
            if(~isnumeric(value)|length(value)~=1)
                fprintf(2,'%s E: Value for %s must be a numeric scalar\n',caller,name);
                return;
            end
        case 'timestep'
            if(~isnumeric(value)|length(value)~=1)
                fprintf(2,'%s E: Value for %s must be a numeric scalar\n',caller,name);
                return;
            end
            timestep=value/24;
        case 'filterlength'
            if(~isnumeric(value)|length(value)~=1)
                fprintf(2,'%s E: Value for %s must be a numeric scalar\n',caller,name);
                return;
            end
            filttime=value;
        case 'tail'
            if(~isnumeric(value))
                fprintf(2,'%s E: Value for %s must be numeric\n',caller,name)
                return;
            end
            if(length(value)~=2)
                fprintf(2,'%s E: Value for %s must have two elements\n',caller,name);
                return;
            end
            if(value(2)<value(1))
                fprintf(2,'%s E: 2nd element in %s must be greater than the first\n',caller,name);
                return;
            end
            tail=value;
        case 'timeformat'
        case 'currenttype'
            if(~ischar(value))
                fprintf(2,'%s E: Value for %s must be a string\n',caller,name);
                return;
            end
            if(value(1)~='X')
                fprintf(2,'%s E: Value for %s must begin with ''X'' (will be replaced with ''u'' and ''v'' in variable name)\n',caller,name);
                return;
            end
            inputvariable=value;
        case 'addlowpass'
            if(~islogical(value))
                fprintf(2,'%s E: Value for %s must be logical\n',caller,name);
                return;
            end
            addlowpass=value;
        case 'addhighpass'
            if(~islogical(value))
                fprintf(2,'%s E: Value for %s must be logical\n',caller,name);
                return;
            end
            addhighpass=value;
        case 'overwrite'
            if(~islogical(value))
                fprintf(2,'%s E: Value for %s must be logical\n',caller,name);
                return;
            end
            if(value&~addlowpass&~addhighpass)
                fprintf(1,'%s Warning: Either ''addlowpass'' or ''addhighpass'' must be set to true for %s to be used; no data will be added to file\n',caller,name);
            end
            replacevals=value;
            if(addlowpass&~replacevals)
                fprintf(1,'%s Warning: Low-pass filtered velocity already exists in %s and ''replace'' is set to false; low-pass filter will not be added\n',caller,inputfile);
            end
            if(addhighpass&~replacevals)
                fprintf(1,'%s Warning: High-pass filtered velocity already exists in %s and ''replace'' is set to false; high-pass filter will not be added\n',caller,inputfile);
            end
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
        otherwise
            fprintf(2,...
                '%s E: Unknown option specified: %s\n',...
                caller,...
                name);
            return;
    end
end

if endtime<=starttime
    fprintf(2,'%s E: endtime must be later than starttime.\n',caller);
    return;
end

lowexist=false;
highexist=false;
varexist=false;

for k=1:length(inputinfo.Variables)
    if(strcmp(inputinfo.Variables(k).Name,['u' inputvariable(2:end) '_lowpass']))
        lowexist=true;
    elseif(strcmp(inputinfo.Variables(k).Name,['u' inputvariable(2:end) '_highpass']))
        highexist=true;
    elseif(strcmp(inputinfo.Variables(k).Name,['u' inputvariable(2:end)]))
        varexist=true;
    end
end

if(~varexist)
    fprintf(2,'%s E: Variables ''u %s'' and ''v %s'' not found in %s\n',caller,inputvariable(2:end),inputvariable(2:end),inputfile);
    return;
end


times=starttime:timestep:endtime;
tind=find(times==currenttime);
if(isempty(tind))
    d=abs(times-currenttime);
    if(min(d<.01))
        tind=find(d==min(d));
    end
end

% open file with raw data and get ID #s for variables
rawnc=netcdf.open(inputfile,'NC_NOWRITE');
rawtid=netcdf.inqVarID(rawnc,'time');
tunits=netcdf.getAtt(rawnc,rawtid,'units');
rawlonid=netcdf.inqVarID(rawnc,'lon');
rawlatid=netcdf.inqVarID(rawnc,'lat');
rawuid=netcdf.inqVarID(rawnc,['u' inputvariable(2:end)]);
ufill=netcdf.getAtt(rawnc,rawuid,'_FillValue');
rawvid=netcdf.inqVarID(rawnc,['v' inputvariable(2:end)]);
vfill=netcdf.getAtt(rawnc,rawvid,'_FillValue');

% read in time from raw hourly file
t=netcdf.getVar(rawnc,rawtid);

if(isempty(t))
    fprintf(2,'%s E: %s has no data\n',caller,inputfile);
    netcdf.close(rawnc)
    return;
end

% read in lon & lat
lon=netcdf.getVar(rawnc,rawlonid);
lat=netcdf.getVar(rawnc,rawlatid);
udescription=netcdf.getAtt(rawnc,rawuid,'longname');
vdescription=netcdf.getAtt(rawnc,rawvid,'longname');
netcdf.close(rawnc)

indlon=find(lon>=min(lonbounds)&lon<=max(lonbounds));
if(isempty(indlon))
    indlon=find(lon>=min(lonbounds+360)&lon<=max(lonbounds+360));
    if(isempty(indlon))
        indlon=find(lon>=min(lonbounds-360)&lon<=max(lonbounds-360));
    end
end
indlat=find(lat>=min(latbounds)&lat<=max(latbounds));
if(isempty(indlon)|isempty(indlat))
    fprintf(2,'%s E: no gridpoints found between latitudes %f and %f and longitudes %f and %f\n',caller,min(latbounds),max(latbounds),min(lonbounds),max(lonbounds));
    return;
end

sizelon=length(lon);
sizelat=length(lat);

lon=lon(indlon);
lat=lat(indlat);


timeoptions={'days','hours','minutes','seconds'};
indt=find(strncmpi(timeoptions,tunits,4));
tunits=timeoptions{indt};

u=nan(length(lon),length(lat),length(times));
v=u;
filtu=u;
filtv=u;

for ti=1:length(times)
    t=times(ti);
    fname=dir([inputdir '*' datestr(t,filetimeformat) '*.nc']);
    if(isempty(fname))
        [yt,mt,dt,Ht,Mt,St]=datevec(t);
        if(timestep>=1/24/60)
            if(St>=30)
                Mt=Mt+1;
            end
            t=datenum(yt,mt,dt,Ht,Mt,0);
            fname=dir([inputdir '*' datestr(t,filetimeformat) '*.nc']);
            if(isempty(fname))
                if(timestep>=1/24/6)
                    Mt=round(Mt/10)*10;
                    t=datenum(yt,mt,dt,Ht,Mt,0);
                    fname=dir([inputdir '*' datestr(t,filetimeformat) '*.nc']);
                end
            end
        end
    end
    if(length(fname)==1)
        ncid=netcdf.open([inputdir fname.name],'NC_NOWRITE');
        tID=netcdf.inqVarID(ncid,'time');
        codartime=netcdf.getVar(ncid,tID);
        if(~isempty(codartime))
            inputinfo=ncinfo([inputdir fname.name]);
            for k=1:length(inputinfo.Variables)
                if(strcmp(inputinfo.Variables(k).Name,['u' inputvariable(2:end)]))
                    uid=netcdf.inqVarID(ncid,['u' inputvariable(2:end)]);
                    ufill=netcdf.getAtt(ncid,uid,'_FillValue');
                    vid=netcdf.inqVarID(ncid,['v' inputvariable(2:end)]);
                    vfill=netcdf.getAtt(ncid,vid,'_FillValue');
                    
                    ut=netcdf.getVar(ncid,uid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1]);
                    ut(ut==ufill)=nan;
                    vt=netcdf.getVar(ncid,vid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1]);
                    vt(vt==vfill)=nan;
                    
                    u(:,:,ti)=ut;
                    v(:,:,ti)=vt;
                end
            end
        else times(ti)=nan;
        end
        netcdf.close(ncid)
    else times(ti)=nan;
    end
end

timestep=mode(diff(times(~isnan(times))));
switch tunits
    case 'hours'
        timestep=timestep/24;
    case 'minutes'
        timestep=timestep/24/60;
    case 'seconds'
        timestep=timestep/24/60/60;
end

% low-pass filter of detided data
[x,y]=cheb1(1,.25,(0.5/(filttime/(24*timestep))));

for i=1:length(lon)
    for j=1:length(lat)
        
        uloc=u(i,j,:);
        vloc=v(i,j,:);
        
        ind=find(~isnan(uloc));
        if(length(find(ind>tind&ind<tind+tail(2)))>tail(1)&length(find(ind<tind&ind>tind-tail(2)))>tail(1))
            detuf=reshape(uloc(ind),[length(ind),1]);
            detvf=reshape(vloc(ind),[length(ind),1]);
            filters=filtfilt(x,y,[detuf detvf]);
            filtu(i,j,ind)=filters(:,1);
            filtv(i,j,ind)=filters(:,2);
        end
    end
end

originalu=u(:,:,tind);
originalv=v(:,:,tind);
lowpassu=filtu(:,:,tind);
lowpassv=filtv(:,:,tind);
highpassu=originalu-lowpassu;
highpassv=originalv-lowpassv;


if(addlowpass|addhighpass)
    outnc=netcdf.open(inputfile,'NC_WRITE');
    outlondimid=netcdf.inqDimID(outnc,'lon');
    outlatdimid=netcdf.inqDimID(outnc,'lat');
    outtdimid=netcdf.inqDimID(outnc,'time');
    
    if(addlowpass)
        if(~lowexist)
            netcdf.reDef(outnc)
            % define low-pass filtered u variable
            outufiltid=netcdf.defVar(outnc,['u' inputvariable(2:end) '_lowpass'],'float',[outlondimid outlatdimid outtdimid]);
            netcdf.putAtt(outnc,outufiltid,'longname', 'Eastward Low-Pass Filtered Velocity');
            netcdf.putAtt(outnc,outufiltid,'units','cm/s');
            netcdf.putAtt(outnc,outufiltid,'description',[int2str(filttime) '-hour low-pass filter of ' udescription]);
            netcdf.putAtt(outnc,outufiltid,'_FillValue',single(-999.0));
            netcdf.defVarFill(outnc,outufiltid,false,single(-999));
            
            % define low-pass filtered v variable
            outvfiltid=netcdf.defVar(outnc,['v' inputvariable(2:end) '_lowpass'],'float',[outlondimid outlatdimid outtdimid]);
            netcdf.putAtt(outnc,outvfiltid,'longname', 'Northward Low-Pass Filtered Velocity');
            netcdf.putAtt(outnc,outvfiltid,'units','cm/s');
            netcdf.putAtt(outnc,outvfiltid,'description',[int2str(filttime) '-hour low-pass filter of ' vdescription]);
            netcdf.putAtt(outnc,outvfiltid,'_FillValue',single(-999.0));
            netcdf.defVarFill(outnc,outvfiltid,false,single(-999));
            netcdf.endDef(outnc)
            
        end
        
        if(lowexist&replacevals)
            outufiltid=netcdf.inqVarID(outnc,['u' inputvariable(2:end) '_lowpass']);
            outvfiltid=netcdf.inqVarID(outnc,['v' inputvariable(2:end) '_lowpass']);
            netcdf.reDef(outnc)
            netcdf.delAtt(outnc,outufiltid,'description');
            netcdf.putAtt(outnc,outufiltid,'description',[int2str(filttime) '-hour low-pass filter of ' udescription]);
            netcdf.delAtt(outnc,outvfiltid,'description');
            netcdf.putAtt(outnc,outvfiltid,'description',[int2str(filttime) '-hour low-pass filter of ' vdescription]);
            netcdf.endDef(outnc)
        end
        
        if(~lowexist|replacevals)
            lowpassu(isnan(lowpassu))=-999;
            lowpassv(isnan(lowpassv))=-999;
            netcdf.putVar(outnc,outufiltid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1],single(lowpassu));
            netcdf.putVar(outnc,outvfiltid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1],single(lowpassv));
            lowpassu(lowpassu==-999)=nan;
            lowpassv(lowpassv==-999)=nan;
        end
    end
    
    if(addhighpass)
        if(~highexist)
            netcdf.reDef(outnc)
            % define high-pass filtered u variable
            outuhighfiltid=netcdf.defVar(outnc,['u' inputvariable(2:end) '_highpass'],'float',[outlondimid outlatdimid outtdimid]);
            netcdf.putAtt(outnc,outuhighfiltid,'longname', 'Eastward High-Pass Filtered Velocity');
            netcdf.putAtt(outnc,outuhighfiltid,'units','cm/s');
            netcdf.putAtt(outnc,outuhighfiltid,'description',[int2str(filttime) '-hour high-pass filter of ' udescription]);
            netcdf.putAtt(outnc,outuhighfiltid,'_FillValue',single(-999.0));
            netcdf.defVarFill(outnc,outuhighfiltid,false,single(-999));
            
            % define high-pass filtered v variable
            outvhighfiltid=netcdf.defVar(outnc,['v' inputvariable(2:end) '_highpass'],'float',[outlondimid outlatdimid outtdimid]);
            netcdf.putAtt(outnc,outvhighfiltid,'longname', 'Northward High-Pass Filtered Velocity');
            netcdf.putAtt(outnc,outvhighfiltid,'units','cm/s');
            netcdf.putAtt(outnc,outvhighfiltid,'description',[int2str(filttime) '-hour high-pass filter of ' vdescription]);
            netcdf.putAtt(outnc,outvhighfiltid,'_FillValue',single(-999.0));
            netcdf.defVarFill(outnc,outvhighfiltid,false,single(-999));
            netcdf.endDef(outnc)
            
        end
        
        if(highexist&replacevals)
            outuhighfiltid=netcdf.inqVarID(outnc,['u' inputvariable(2:end) '_highpass']);
            outvhighfiltid=netcdf.inqVarID(outnc,['v' inputvariable(2:end) '_highpass']);
            netcdf.reDef(outnc)
            netcdf.delAtt(outnc,outuhighfiltid,'description');
            netcdf.putAtt(outnc,outuhighfiltid,'description',[int2str(filttime) '-hour high-pass filter of ' udescription]);
            netcdf.delAtt(outnc,outvhighfiltid,'description');
            netcdf.putAtt(outnc,outvhighfiltid,'description',[int2str(filttime) '-hour high-pass filter of ' vdescription]);
            netcdf.endDef(outnc)
        end
        
        if(~highexist|replacevals)
            highpassu(isnan(highpassu))=-999;
            highpassv(isnan(highpassv))=-999;
            netcdf.putVar(outnc,outuhighfiltid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1],single(highpassu));
            netcdf.putVar(outnc,outvhighfiltid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1],single(highpassv));
            highpassu(highpassu==-999)=nan;
            highpassv(highpassv==-999)=nan;
        end
    end
    
    netcdf.close(outnc)
end

