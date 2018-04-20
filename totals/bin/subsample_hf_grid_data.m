function [BPU_grid_ind]=subsample_hf_grid_data(x,y,spacing)
%% This function subsets the MARCOOS HF Radar grid data
%% Subset the data so only certain vectors are plotted

%-------------------------------------------------------------------------
% Rutgers University Coastal Ocean Observation Lab
%------------------------------------------------------------------------
%
%Usage:
%   [marcoos_grid_ind]=subsample_hf_grid_data(x,y)
%Input:
%   x,y = lon and lat points of the grid you want to subsample
%Output:
%   marcoos_grid_ind = indices of the grid points you want to keep
%Updates:
%   Written by Hugh Roarty 3/16/2010
%Notes:
%   Assumes a 130 by 138 grid



% Find the unique values for the lat and lon grid
unique_x=unique(x);
unique_y=unique(y);

% The indices for the rows and columns to keep
% x_ind=[1;10;19;27;36;44;53;61;70;78;87;96;105;113;122;130;138];
% y_ind=[1;;10;20;29;38;47;57;66;76;85;94;103;113;121;130];

%x_ind=[1:3:138];
%y_ind=[1:3:130];

[NX,NY]=meshgrid(unique_x(1:spacing:end),unique_y(1:spacing:end));

% Vectorize the array
NX=NX(:);
NY=NY(:);

% find the indices of the grid points you want to keep
for i=1:length(NX)
    ind=find(NX(i)==x &NY(i)==y);
    if ~isempty(ind)
        BPU_grid_ind(i)=ind;
    end
end

% remove the zeros from the indices array
BPU_grid_ind(BPU_grid_ind==0)=[];