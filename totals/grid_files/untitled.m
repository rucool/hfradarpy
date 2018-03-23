load('/Users/mikesmith/Documents/data/tuv_oi_MARA_2017_07_11_2300.mat');
grid = load('OI_6km_Grid.txt');

[X, Y] = meshgrid(unique(grid(:,1)), unique(grid(:,2)));

[u_test, notOnGrid] = griddata_nointerp(X, Y, TUV.LonLat(:,1), TUV.LonLat(:,2), TUV.U, .01, .01);