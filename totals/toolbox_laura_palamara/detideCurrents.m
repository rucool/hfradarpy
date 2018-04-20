function [lon,lat,rawu,rawv,tidalu,tidalv,detidedu,detidedv]=detideCurrents(inputfile,coefficientfile,varargin)

% Get tidal and detided current velocities.
%
% Outputs:
%     lon: longitude vector for grid
%     lat: latitude vector for grid
%     rawu: matrix of eastward velocity totals (from inputfile)
%     rawv: matrix of northward velocity totals (from inputfile)
%     tidalu: matrix of tidal eastward velocity estimated from
%       coefficients in coefficientfile
%     tidalv: matrix of tidal northward velocity estimated from
%       coefficients in coefficientfile
%     detidedu: matrix of estimated detided eastward velocity
%     detidedv: matrix of estimated detided northward velocity
%
% Inputs:
%     inputfile: file name of netcdf that includes the raw current totals
%     coefficientfile: file name of netcdf that includes detiding
%       coefficients (where the time format matches the time format in
%       inputfile)
%
% Options:
%     mintotalcoverage: amount of files that had valid data in the detiding
%       process (default: 336 = 2 weeks of hourly data)
%     maxerror: maximum allowable error value for raw u or v (default: .6)
%     addtide: logical value indicating whether to add tidal velocity
%       component to the original surface current netcdf file (inputfile)
%       (default: false)
%     adddetided: logical value indicating whether to add detided velocity
%       component to the original surface current netcdf file (inputfile)
%       (default: false)
%     overwrite: logical value indicating whether to replace tidal velocity
%       and detided velocity if they already exist in the netcdf (addtide
%       and/or adddetided must be set to true; default: false)
%     lonbounds: min and max longitude to read in (formatted [minlon
%       maxlon]) (default: read in all)
%     latbounds: min and max latitude to read in (formatted [minlat
%       maxlat]) (default: read in all)
%

lon=[];
lat=[];
rawu=[];
rawv=[];
tidalu=[];
tidalv=[];
detidedu=[];
detidedv=[];

caller = [mfilename '.m'];

if nargin <2
    fprintf(2,...
        '%s E: not enough arguments.\n',...
        caller);
    return;
elseif ~ischar(inputfile) | ~ischar(coefficientfile)
    fprintf(2,...
        '%s E: inputfile and coefficientfile must both be strings.\n',...
        caller);
    return;
elseif ~exist(inputfile,'file')
    fprintf(2,'%s E: inputfile %s not found.\n',caller,inputfile);
    return;
elseif ~exist(coefficientfile,'file')
    fprintf(2,'%s E: coefficientfile %s not found.\n',caller,coefficientfile);
    return;
end


maxerr=.6;
minp=24*14;
addtide=false;
adddetided=false;
replacevals=false;
lonbounds=[-360 360];
latbounds=[-90 90];

inputinfo=ncinfo(inputfile);
tideexist=false;
detideexist=false;

for k=1:length(inputinfo.Variables)
    if(strcmp(inputinfo.Variables(k).Name,'u_tide'))
        tideexist=true;
    elseif(strcmp(inputinfo.Variables(k).Name,'u_detided'))
        detideexist=true;
    end
end

for x=1:2:length(varargin)
    name=varargin{x};
    value=varargin{x+1};
    switch lower(name)
        case 'mintotalcoverage'
            if(~isnumeric(value)|length(value)>1)
                fprintf(2,'%s E: Value for %s must be a numeric scalar\n',caller,name);
                return;
            end
            if(value<1)
                fprintf(2,'%s E: Value for %s must be greater than 1\n',caller,name);
                return;
            end
            minp=value;
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
        case 'addtide'
            if(~islogical(value))
                fprintf(2,'%s E: Value for %s must be logical\n',caller,name);
                return;
            end
            addtide=value;
        case 'adddetided'
            if(~islogical(value))
                fprintf(2,'%s E: Value for %s must be logical\n',caller,name);
                return;
            end
            adddetided=value;
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
            if(value&~addtide&~adddetided)
                fprintf(1,'%s Warning: Either ''addtide'' or ''adddetided'' must be set to true for %s to be used; no data will be added to file\n',caller,name);
            end
            replacevals=value;
            if(addtide&~replacevals)
                fprintf(1,'%s Warning: Tidal component already exists in %s and ''replace'' is set to false; tidal component will not be added\n',caller,inputfile);
            end
            if(adddetided&~replacevals)
                fprintf(1,'%s Warning: Detided component already exists in %s and ''replace'' is set to false; detided component will not be added\n',caller,inputfile);
            end
        otherwise
            fprintf(2,...
                '%s E: Unknown option specified: %s\n',...
                caller,...
                name);
            return;
    end
end


% open file with raw data and get ID #s for variables
rawnc=netcdf.open(inputfile,'NC_NOWRITE');
rawtid=netcdf.inqVarID(rawnc,'time');
rawlonid=netcdf.inqVarID(rawnc,'lon');
rawlatid=netcdf.inqVarID(rawnc,'lat');
rawuid=netcdf.inqVarID(rawnc,'u');
ufill=netcdf.getAtt(rawnc,rawuid,'_FillValue');
rawvid=netcdf.inqVarID(rawnc,'v');
vfill=netcdf.getAtt(rawnc,rawvid,'_FillValue');
rawuerrid=netcdf.inqVarID(rawnc,'u_err');
rawverrid=netcdf.inqVarID(rawnc,'v_err');

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

indlon=find(lon>=min(lonbounds)&lon<=max(lonbounds));
if(isempty(indlon))
    indlon=find(lon>=min(lonbounds+360)&lon<=max(lonbounds+360));
    if(isempty(indlon))
        indlon=find(lon>=min(lonbounds-360)&lon<=max(lonbounds-360));
    end
end
indlat=find(lat>=min(latbounds)&lat<=max(latbounds));
if(isempty(indlon)|isempty(indlat))
    netcdf.close(rawnc)
    fprintf(2,'%s E: no gridpoints found between latitudes %f and %f and longitudes %f and %f\n',caller,min(latbounds),max(latbounds),min(lonbounds),max(lonbounds));
    return;
end

sizelon=length(lon);
sizelat=length(lat);

lon=lon(indlon);
lat=lat(indlat);

% pull in raw u and v (and errors)
rawu=netcdf.getVar(rawnc,rawuid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1]);
rawv=netcdf.getVar(rawnc,rawvid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1]);
uerr=netcdf.getVar(rawnc,rawuerrid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1]);
verr=netcdf.getVar(rawnc,rawverrid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1]);
netcdf.close(rawnc)

% find gridpoints where the error for either u or v is >.6 and change
% velocities to NaNs
inderr=find(uerr>maxerr|verr>maxerr);
clear uerr verr
rawu(inderr)=NaN;
rawv(inderr)=NaN;
rawu(rawu==ufill)=NaN;
rawv(rawv==vfill)=NaN;

% open file with detiding coefficients and find ID #s for variables
coeffnc=netcdf.open(coefficientfile,'NC_NOWRITE');
coeffuid=netcdf.inqVarID(coeffnc,'solnu');
coeffvid=netcdf.inqVarID(coeffnc,'solnv');
coeffnpid=netcdf.inqVarID(coeffnc,'np');
coeffperid=netcdf.inqVarID(coeffnc,'period');

period=netcdf.getVar(coeffnc,coeffperid);
periodunits=netcdf.getAtt(coeffnc,coeffperid,'units');
switch periodunits
    case 'days'
        period=period*24;
    case 'minutes'
        period=period/60;
    case 'seconds'
        period=period/60/60;
end

omega=(2*24*pi)./period;
nc=length(period);

% read in detiding coefficients and # points in timeseries used for
% detiding
coeffsnp=netcdf.getVar(coeffnc,coeffnpid,[min(indlon)-1 min(indlat)-1],[length(indlon) length(indlat)]);
solnuall=netcdf.getVar(coeffnc,coeffuid,[0 min(indlon)-1 min(indlat)-1],[length(period)*2+1 length(indlon) length(indlat)]);
solnvall=netcdf.getVar(coeffnc,coeffvid,[0 min(indlon)-1 min(indlat)-1],[length(period)*2+1 length(indlon) length(indlat)]);
detidingrange=netcdf.getAtt(coeffnc,-1,'date_range');
netcdf.close(coeffnc)

% create default NaN grid and assign to tidal and detided velocity
% grids
tideu=nan(size(rawu));
tidev=tideu;
detu=tideu;
detv=tideu;

% loop through grid to detide
for i=1:length(lon)
    % only detide points with at least 2 weeks' data used to train the
    % detiding coefficients
    inddet=find(coeffsnp(i,:)>=minp);
    
    for ind=1:length(inddet)
        j=inddet(ind);
        
        if(~isnan(rawu(i,j)))
            % get coefficients for that point
            solnu=solnuall(:,i,j);
            solnv=solnvall(:,i,j);
            
            
            soln=solnu';
            a=zeros(1,nc);
            b=zeros(1,nc);
            for ic=1:1:nc
                a(ic)=soln(2*ic);
                b(ic)=soln(2*ic+1);
            end
            
            % get tidal u (subtract from raw to get detided)
            pred=zeros(1,1);
            for jc=1:1:nc
                pred=pred+(a(jc)*cos(omega(jc)*t))+(b(jc)*sin(omega(jc)*t));
            end
            
            tideu(i,j)=pred;
            detu(i,j)=rawu(i,j)-pred;
            
            
            soln=solnv';
            a=zeros(1,nc);
            b=zeros(1,nc);
            for ic=1:1:nc
                a(ic)=soln(2*ic);
                b(ic)=soln(2*ic+1);
            end
            
            % get tidal v (subtract from raw to get detided)
            pred=zeros(1,1);
            for jc=1:1:nc
                pred=pred+(a(jc)*cos(omega(jc)*t))+(b(jc)*sin(omega(jc)*t));
            end
            
            tidev(i,j)=pred;
            detv(i,j)=rawv(i,j)-pred;
            
        end
    end
end

if(addtide|adddetided)
    outnc=netcdf.open(inputfile,'NC_WRITE');
    outlondimid=netcdf.inqDimID(outnc,'lon');
    outlatdimid=netcdf.inqDimID(outnc,'lat');
    outtdimid=netcdf.inqDimID(outnc,'time');
    perstr=[];
    for k=1:length(period)
        perstr=[perstr num2str(period(k)) ' '];
    end
    
    if(addtide)
        if(~tideexist)
            netcdf.reDef(outnc)
            % define tidal u variable
            oututideid=netcdf.defVar(outnc,'u_tide','float',[outlondimid outlatdimid outtdimid]);
            netcdf.putAtt(outnc,oututideid,'longname', 'Eastward Tidal Velocity');
            netcdf.putAtt(outnc,oututideid,'units','cm/s');
            netcdf.putAtt(outnc,oututideid,'description',['tidal u velocity based on least-squares fit of tidal periods ' perstr 'hours']);
            netcdf.putAtt(outnc,oututideid,'detiding_date_range',detidingrange);
            netcdf.putAtt(outnc,oututideid,'_FillValue',single(-999.0));
            netcdf.defVarFill(outnc,oututideid,false,single(-999));
            
            % define tidal v variable
            outvtideid=netcdf.defVar(outnc,'v_tide','float',[outlondimid outlatdimid outtdimid]);
            netcdf.putAtt(outnc,outvtideid,'longname', 'Northward Tidal Velocity');
            netcdf.putAtt(outnc,outvtideid,'units','cm/s');
            netcdf.putAtt(outnc,outvtideid,'description',['tidal v velocity based on least-squares fit of tidal periods ' perstr 'hours']);
            netcdf.putAtt(outnc,outvtideid,'detiding_date_range',detidingrange);
            netcdf.putAtt(outnc,outvtideid,'_FillValue',single(-999.0));
            netcdf.defVarFill(outnc,outvtideid,false,single(-999));
            netcdf.endDef(outnc)
        end
        
        if(tideexist&replacevals)
            oututideid=netcdf.inqVarID(outnc,'u_tide');
            outvtideid=netcdf.inqVarID(outnc,'v_tide');
            netcdf.reDef(outnc)
            netcdf.delAtt(outnc,oututideid,'detiding_date_range');
            netcdf.putAtt(outnc,oututideid,'detiding_date_range',detidingrange);
            netcdf.delAtt(outnc,outvtideid,'detiding_date_range');
            netcdf.putAtt(outnc,outvtideid,'detiding_date_range',detidingrange);
            netcdf.endDef(outnc)
        end
        
        if(~tideexist|replacevals)
            tideu(isnan(tideu))=-999;
            tidev(isnan(tidev))=-999;
            netcdf.putVar(outnc,oututideid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1],single(tideu));
            netcdf.putVar(outnc,outvtideid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1],single(tidev));
            tideu(tideu==-999)=nan;
            tidev(tidev==-999)=nan;
        end
    end
    
    if(adddetided)
        if(~detideexist)
            netcdf.reDef(outnc)
            % define detided u variable
            outudetid=netcdf.defVar(outnc,'u_detided','float',[outlondimid outlatdimid outtdimid]);
            netcdf.putAtt(outnc,outudetid,'longname', 'Eastward Detided Velocity');
            netcdf.putAtt(outnc,outudetid,'units','cm/s');
            netcdf.putAtt(outnc,outudetid,'description','raw u velocity minus tidal u velocity');
            netcdf.putAtt(outnc,outudetid,'detiding_date_range',detidingrange);
            netcdf.putAtt(outnc,outudetid,'_FillValue',single(-999.0));
            netcdf.defVarFill(outnc,outudetid,false,single(-999));
            
            % define detided v variable
            outvdetid=netcdf.defVar(outnc,'v_detided','float',[outlondimid outlatdimid outtdimid]);
            netcdf.putAtt(outnc,outvdetid,'longname', 'Northward Detided Velocity');
            netcdf.putAtt(outnc,outvdetid,'units','cm/s');
            netcdf.putAtt(outnc,outvdetid,'description','raw v velocity minus tidal v velocity');
            netcdf.putAtt(outnc,outvdetid,'detiding_date_range',detidingrange);
            netcdf.putAtt(outnc,outvdetid,'_FillValue',single(-999.0));
            netcdf.defVarFill(outnc,outvdetid,false,single(-999));
            netcdf.endDef(outnc)
        end
        
        if(detideexist&replacevals)
            outudetid=netcdf.inqVarID(outnc,'u_detided');
            outvdetid=netcdf.inqVarID(outnc,'v_detided');
            netcdf.reDef(outnc)
            netcdf.delAtt(outnc,outudetid,'detiding_date_range');
            netcdf.putAtt(outnc,outudetid,'detiding_date_range',detidingrange);
            netcdf.delAtt(outnc,outvdetid,'detiding_date_range');
            netcdf.putAtt(outnc,outvdetid,'detiding_date_range',detidingrange);
            netcdf.endDef(outnc)
        end
        
        if(~detideexist|replacevals)
            detu(isnan(detu))=-999;
            detv(isnan(detv))=-999;
            netcdf.putVar(outnc,outudetid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1],single(detu));
            netcdf.putVar(outnc,outvdetid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1],single(detv));
            detu(detu==-999)=nan;
            detv(detv==-999)=nan;
        end
    end
    
    netcdf.close(outnc)
end

