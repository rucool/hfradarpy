function [lon,lat,meandata,trendeddata]=trendData(inputfile,inputvariable,threshold,varargin)

% Get tidal and detided current velocities.
%
% Outputs:
%     lon: longitude vector for grid
%     lat: latitude vector for grid
%     meandata: matrix of absolute mean of data specified over time period
%       specified
%     trendeddata: matrix of "trend" of data specified over time period
%       specified, calculated by averaging together +1s for "high" (>high
%       threshold), -1s for "low" (<low threshold), and 0s for "neutral"
%       (between thresholds), instead of the instantaneous values
%
% Inputs:
%     inputfile: file name of netcdf that includes the data
%     inputvariable: string of variable name to read in from netcdf
%     threshold: 1- or 2- element vector with 1st element as the "low"
%       threshold and the 2nd element as the "high" threshold; set both
%       elements the same to have all values defined as high or low except
%       those exactly equal to the threshold; 1-element vectors must be
%       positive and assume a high threshold of that number and a low
%       threshold of the negative
%
% Options:
%     timeformat: format of time in file name of surface current total
%       netcdfs (default: 'yyyy_mm_dd_HHMM')
%     averagingtime: amount of time over which to get the mean and trend
%       (default: 7 days)
%     averagingformat: which direction to average through time; valid
%       options:
%           'previous': averages from the time of the given file BACK
%             through the amount of time provided (default)
%           'centered': averages from half the time provided before the
%             given file through half the time provided after the given
%             file
%           'future': averages from the time of the given file FORWARD
%             through the amount of time provided
%     addmean: logical value indicating whether to add time-averaged data
%       to the original surface current netcdf file (inputfile) (default:
%       false)
%     addtrend: logical value indicating whether to add trend of data
%       through time to the original surface current netcdf file
%       (inputfile) (default:false)
%     overwrite: logical value indicating whether to replace mean
%       and trend if they already exist in the netcdf (addmean
%       and/or addtrend must be set to true; default: false)
%     lonbounds: min and max longitude to read in (formatted [minlon
%       maxlon]) (default: read in all)
%     latbounds: min and max latitude to read in (formatted [minlat
%       maxlat]) (default: read in all)
%

lon=[];
lat=[];
meandata=[];
trendeddata=[];

caller = [mfilename '.m'];

if nargin <3
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
elseif ~ischar(inputvariable)
    fprintf(2,'%s E: inputvariable must be a string.\n',caller);
    return;
elseif ~isnumeric(threshold)
    fprintf(2,'%s E: threshold must be numeric.\n',caller);
    return;
elseif length(threshold)<1|length(threshold)>2
    fprintf(2,'%s E: threshold must be 1 or 2 elements.\n',caller);
    return;
end

inputinfo=ncinfo(inputfile);
varexist=false;
for k=1:length(inputinfo.Variables)
    if(strcmp(inputinfo.Variables(k).Name,inputvariable))
        varexist=true;
    end
end
if(~varexist)
    fprintf(2,'%s E: variable %s does not exist in file %s.\n',caller,inputvariable,inputfile);
    return;
end

if(length(threshold)==1)
    if(threshold<0)
        fprintf(2,'%s E: Value for threshold must be positive if only one element\n',caller);
        return;
    end
    threshold=[-threshold threshold];
    fprintf(1,['%s Warning: Only one threshold value provided; assuming values less than ' num2str(threshold(1)) ' are low and values greater than ' num2str(threshold(2)) ' are high.\n'],caller);
end
if(length(threshold)~=2)
    fprintf(2,'%s E: Value for threshold must be either 1 or 2 elements\n',caller);
    return;
end
if(threshold(2)<threshold(1))
    fprintf(2,'%s E: 2nd element of threshold must be greater than or equal to the first\n',caller);
    return;
end
if(threshold(1)==threshold(2))
    fprintf(1,'%s Warning: 1st and 2nd elements of threshold are equal; only values exactly equal to %f will be considered "neutral".\n',caller,threshold(1));
end


filetimeformat='yyyy_mm_dd_HHMM';
avgtime=7;
avgformat='previous';
addavg=false;
addtrend=false;
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

avgexist=false;
trendexist=false;


validformats={'previous','centered','future'};


for x=1:2:length(varargin)
    name=varargin{x};
    value=varargin{x+1};
    switch lower(name)
        case 'timeformat'
        case 'averagingtime'
            if(~isnumeric(value)|length(value)~=1)
                fprintf(2,'%s E: Value for %s must be a numeric scalar\n',caller,name);
                return;
            end
            avgtime=value;
        case 'averageingformat'
            if(~ischar(value))
                fprintf(2,'%s E: Value for %s must be a string\n',caller,name);
                return;
            end
            if(~ismember(value,validformats))
                fprintf(2,'%s E: Value for %s is not recognized as a valid format\n',caller,name);
                return;
            end
            avgformat=lower(value);
        case 'addmean'
            if(~islogical(value))
                fprintf(2,'%s E: Value for %s must be logical\n',caller,name);
                return;
            end
            addavg=value;
        case 'addtrend'
            if(~islogical(value))
                fprintf(2,'%s E: Value for %s must be logical\n',caller,name);
                return;
            end
            addtrend=value;
        case 'overwrite'
            if(~islogical(value))
                fprintf(2,'%s E: Value for %s must be logical\n',caller,name);
                return;
            end
            if(value&~addavg&~addtrend)
                fprintf(1,'%s Warning: Either ''addavg'' or ''addtrend'' must be set to true for %s to be used; no data will be added to file\n',caller,name);
            end
            replacevals=value;
            if(addavg&~replacevals)
                fprintf(1,'%s Warning: Mean %s already exists in %s and ''replace'' is set to false; mean will not be added\n',caller,inputvariable,inputfile);
            end
            if(addtrend&~replacevals)
                fprintf(1,'%s Warning: Trended %s already exists in %s and ''replace'' is set to false; trend will not be added\n',caller,inputvariable,inputfile);
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

switch avgformat
    case 'previous'
        endtime=currenttime;
        starttime=currenttime-avgtime;
    case 'centered'
        endtime=currenttime+avgtime/2;
        starttime=currenttime-avgtime/2;
    case 'future'
        endtime=currenttime+avgtime;
        starttime=currenttime;
end

for k=1:length(inputinfo.Variables)
    if(strcmp(inputinfo.Variables(k).Name,[inputvariable '_mean']))
        avgexist=true;
    elseif(strcmp(inputinfo.Variables(k).Name,[inputvariable '_trend']))
        trendexist=true;
    end
end



allfiles=dir([inputdir fileformatstart '*' fileformatend]);
allfiles=struct2cell(allfiles);
allfiles=char(allfiles(1,:));
indf=find(datenum(allfiles(:,timestart:timestart+length(filetimeformat)-1),filetimeformat)>=starttime&datenum(allfiles(:,timestart:timestart+length(filetimeformat)-1),filetimeformat)<=endtime);

% open file with raw data and get ID #s for variables
ncid=netcdf.open(inputfile,'NC_NOWRITE');
tid=netcdf.inqVarID(ncid,'time');
lonid=netcdf.inqVarID(ncid,'lon');
latid=netcdf.inqVarID(ncid,'lat');
varid=netcdf.inqVarID(ncid,inputvariable);
varfill=netcdf.getAtt(ncid,varid,'_FillValue');

% read in time from raw hourly file
t=netcdf.getVar(ncid,tid);

if(isempty(t))
    fprintf(2,'%s E: %s has no data\n',caller,inputfile);
    netcdf.close(ncid)
    return;
end

% read in lon & lat
lon=netcdf.getVar(ncid,lonid);
lat=netcdf.getVar(ncid,latid);
vardescription=netcdf.getAtt(ncid,varid,'longname');
varunits=netcdf.getAtt(ncid,varid,'units');
netcdf.close(ncid)

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


originaldata=nan(length(lon),length(lat),length(indf));

for i=1:length(indf)
    
    ncid=netcdf.open([inputdir allfiles(indf(i),:)],'NC_NOWRITE');
    tid=netcdf.inqVarID(ncid,'time');
    t=netcdf.getVar(ncid,tid);
    if(~isempty(t))
        inputinfo=ncinfo([inputdir allfiles(indf(i),:)]);
        for k=1:length(inputinfo.Variables)
            if(strcmp(inputinfo.Variables(k).Name,inputvariable))
                varid=netcdf.inqVarID(ncid,inputvariable);
                varfill=netcdf.getAtt(ncid,varid,'_FillValue');
                var=netcdf.getVar(ncid,varid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1]);
                var(var==varfill)=nan;
                originaldata(:,:,i)=var;
            end
        end
    end
    netcdf.close(ncid)
end

meandata=nanmean(originaldata,3);
thresholddata=nan(size(originaldata));
thresholddata(~isnan(originaldata))=0;
thresholddata(originaldata>=threshold(2))=1;
thresholddata(originaldata<=threshold(1))=-1;
trendeddata=nanmean(thresholddata,3);


if(addavg|addtrend)
    outnc=netcdf.open(inputfile,'NC_WRITE');
    outlondimid=netcdf.inqDimID(outnc,'lon');
    outlatdimid=netcdf.inqDimID(outnc,'lat');
    outtdimid=netcdf.inqDimID(outnc,'time');
    
    if(addavg)
        if(~avgexist)
            netcdf.reDef(outnc)
            % define variable
            outmeanid=netcdf.defVar(outnc,[inputvariable '_mean'],'float',[outlondimid outlatdimid outtdimid]);
            netcdf.putAtt(outnc,outmeanid,'longname', ['Mean ' vardescription]);
            netcdf.putAtt(outnc,outmeanid,'units',varunits);
            netcdf.putAtt(outnc,outmeanid,'description',[num2str(avgtime) '-day ' avgformat ' mean ' vardescription]);
            netcdf.putAtt(outnc,outmeanid,'_FillValue',single(-999.0));
            netcdf.defVarFill(outnc,outmeanid,false,single(-999));
            netcdf.endDef(outnc)
        end
        
        if(avgexist&replacevals)
            outmeanid=netcdf.inqVarID(outnc,[inputvariable '_mean']);
            netcdf.reDef(outnc)
            netcdf.delAtt(outnc,outmeanid,'description');
            netcdf.putAtt(outnc,outmeanid,'description',[num2str(avgtime) '-day ' avgformat ' mean ' vardescription]);
            netcdf.endDef(outnc)
        end
        
        if(~avgexist|replacevals)
            meandata(isnan(meandata))=-999;
            netcdf.putVar(outnc,outmeanid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1],single(meandata));
            meandata(meandata==-999)=nan;
        end
    end
    
    if(addtrend)
        if(~trendexist)
            netcdf.reDef(outnc)
            % define vorticity variable
            outtrendid=netcdf.defVar(outnc,[inputvariable '_trend'],'float',[outlondimid outlatdimid outtdimid]);
            netcdf.putAtt(outnc,outtrendid,'longname', ['Trend of ' vardescription]);
            netcdf.putAtt(outnc,outtrendid,'units','unitless');
            netcdf.putAtt(outnc,outtrendid,'description',[num2str(avgtime) '-day ' avgformat ' trended ' vardescription]);
            netcdf.putAtt(outnc,outtrendid,'_FillValue',single(-999.0));
            netcdf.defVarFill(outnc,outtrendid,false,single(-999));
            netcdf.endDef(outnc)
        end
        
        if(trendexist&replacevals)
            outtrendid=netcdf.inqVarID(outnc,[inputvariable '_trend']);
            netcdf.reDef(outnc)
            netcdf.delAtt(outnc,outtrendid,'description');
            netcdf.putAtt(outnc,outtrendid,'description',[num2str(avgtime) '-day ' avgformat ' trended ' vardescription]);
            netcdf.endDef(outnc)
        end
        
        if(~trendexist|replacevals)
            trendeddata(isnan(trendeddata))=-999;
            netcdf.putVar(outnc,outtrendid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1],single(trendeddata));
            trendeddata(trendeddata==-999)=nan;
        end
    end
    
    netcdf.close(outnc)
end

