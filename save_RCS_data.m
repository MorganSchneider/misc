DebrisType(1).Name = 'Leaf';
% DebrisType(1).Dimensions.depth_m = 0.002;
DebrisType(1).Dimensions.depth_m = 0.001;
DebrisType(1).Dimensions.height_m = 0.04;
DebrisType(1).Dimensions.width_m = 0.04;
DebrisType(1).Dimensions.density_kg_m3 = 350.0;
DebrisType(1).Dimensions.string = [num2str(100 * DebrisType(1).Dimensions.depth_m), 'cm x ',...
    num2str(100 * DebrisType(1).Dimensions.height_m), 'cm x ',...
    num2str(100 * DebrisType(1).Dimensions.width_m), 'cm'];

DebrisType(2).Name = 'BigLeaf';
DebrisType(2).Dimensions.depth_m = 0.001;
DebrisType(2).Dimensions.height_m = 0.08;
DebrisType(2).Dimensions.width_m = 0.06;
DebrisType(2).Dimensions.density_kg_m3 = 350.0;
DebrisType(2).Dimensions.string = [num2str(100 * DebrisType(2).Dimensions.depth_m), 'cm x ',...
    num2str(100 * DebrisType(2).Dimensions.height_m), 'cm x ',...
    num2str(100 * DebrisType(2).Dimensions.width_m), 'cm'];

DebrisType(3).Name = 'WoodBoard2x4';
DebrisType(3).Dimensions.depth_m = 0.0508; % 2 in
DebrisType(3).Dimensions.height_m = 0.3048; % 12 in
DebrisType(3).Dimensions.width_m = 0.1016; % 4 in
DebrisType(3).Dimensions.density_kg_m3 = 500.0;
DebrisType(3).Dimensions.string = [num2str(100 / 2.54 * DebrisType(3).Dimensions.depth_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(3).Dimensions.height_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(3).Dimensions.width_m), 'in'];

DebrisType(4).Name = 'WoodBoard4x8';
DebrisType(4).Dimensions.depth_m = 0.1016; % 4 in
DebrisType(4).Dimensions.height_m = 0.3048; % 12 in
DebrisType(4).Dimensions.width_m = 0.2032; % 8 in
DebrisType(4).Dimensions.density_kg_m3 = 500.0;
DebrisType(4).Dimensions.string = [num2str(100 / 2.54 * DebrisType(4).Dimensions.depth_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(4).Dimensions.height_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(4).Dimensions.width_m), 'in'];

DebrisType(5).Name = 'MetalSheet';
DebrisType(5).Dimensions.depth_m = 0.001;
DebrisType(5).Dimensions.height_m = 1.0;
DebrisType(5).Dimensions.width_m = 1.0;
% DebrisType(5).Dimensions.density_kg_m3 = 7850.0;
DebrisType(5).Dimensions.density_kg_m3 = 350.0;
DebrisType(5).Dimensions.string = [num2str(100 * DebrisType(5).Dimensions.depth_m), 'cm x ',...
    num2str(100 * DebrisType(5).Dimensions.height_m), 'cm x ',...
    num2str(100 * DebrisType(5).Dimensions.width_m), 'cm'];

DebrisType(6).Name = 'Brick';
DebrisType(6).Dimensions.depth_m = 0.065;
DebrisType(6).Dimensions.height_m = 0.215;
DebrisType(6).Dimensions.width_m = 0.1125;
DebrisType(6).Dimensions.density_kg_m3 = 2200.0;
DebrisType(6).Dimensions.string = [num2str(100 * DebrisType(6).Dimensions.depth_m), 'cm x ',...
    num2str(100 * DebrisType(6).Dimensions.height_m), 'cm x ',...
    num2str(100 * DebrisType(6).Dimensions.width_m), 'cm'];

DebrisType(7).Name = 'WoodBoard1';
% DebrisType(7).Dimensions.depth_m = 0.0381; % 1.5 in
% DebrisType(7).Dimensions.height_m = 2.4384; % 96 in
% DebrisType(7).Dimensions.width_m = 0.0889; % 3.5 in
DebrisType(7).Dimensions.depth_m = 0.0508; % 2 in
DebrisType(7).Dimensions.height_m = 0.0762; % 3 in
DebrisType(7).Dimensions.width_m = 0.1016; % 4 in
DebrisType(7).Dimensions.density_kg_m3 = 500.0;
DebrisType(7).Dimensions.string = [num2str(100 / 2.54 * DebrisType(7).Dimensions.depth_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(7).Dimensions.height_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(7).Dimensions.width_m), 'in'];

DebrisType(8).Name = 'WoodBoard2';
% DebrisType(8).Dimensions.depth_m = 0.0381; % 1.5 in
% DebrisType(8).Dimensions.height_m = 1.2192; % 48 in
% DebrisType(8).Dimensions.width_m = 0.0889; % 3.5 in
DebrisType(8).Dimensions.depth_m = 0.0508; % 2 in
DebrisType(8).Dimensions.height_m = 0.1524; % 6 in
DebrisType(8).Dimensions.width_m = 0.1016; % 4 in
DebrisType(8).Dimensions.density_kg_m3 = 500.0;
DebrisType(8).Dimensions.string = [num2str(100 / 2.54 * DebrisType(8).Dimensions.depth_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(8).Dimensions.height_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(8).Dimensions.width_m), 'in'];

DebrisType(9).Name = 'WoodBoard3';
% DebrisType(9).Dimensions.depth_m = 0.0381; % 1.5 in
% DebrisType(9).Dimensions.height_m = 0.6096; % 24 in
% DebrisType(9).Dimensions.width_m = 0.0889; % 3.5 in
DebrisType(9).Dimensions.depth_m = 0.0508; % 2 in
DebrisType(9).Dimensions.height_m = 0.3048; % 12 in
DebrisType(9).Dimensions.width_m = 0.1016; % 4 in
DebrisType(9).Dimensions.density_kg_m3 = 500.0;
DebrisType(9).Dimensions.string = [num2str(100 / 2.54 * DebrisType(9).Dimensions.depth_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(9).Dimensions.height_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(9).Dimensions.width_m), 'in'];

DebrisType(10).Name = 'WoodBoard4';
% DebrisType(10).Dimensions.depth_m = 0.0381; % 1.5 in
% DebrisType(10).Dimensions.height_m = 0.3048; % 12 in
% DebrisType(10).Dimensions.width_m = 0.0889; % 3.5 in
DebrisType(10).Dimensions.depth_m = 0.0508; % 2 in
DebrisType(10).Dimensions.height_m = 0.6096; % 24 in
DebrisType(10).Dimensions.width_m = 0.1016; % 4 in
DebrisType(10).Dimensions.density_kg_m3 = 500.0;
DebrisType(10).Dimensions.string = [num2str(100 / 2.54 * DebrisType(10).Dimensions.depth_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(10).Dimensions.height_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(10).Dimensions.width_m), 'in'];

DebrisType(11).Name = 'WoodSheet1';
% DebrisType(11).Dimensions.depth_m = 0.00635; % 1/4 in
% DebrisType(11).Dimensions.height_m = 1.2192; % 48 in
% DebrisType(11).Dimensions.width_m = 1.2192;
DebrisType(11).Dimensions.depth_m = 0.00635; % 1/4 in
DebrisType(11).Dimensions.height_m = 0.0381; % 1.5 in
DebrisType(11).Dimensions.width_m = 0.0381;
DebrisType(11).Dimensions.density_kg_m3 = 500.0;
DebrisType(11).Dimensions.string = [num2str(100 / 2.54 * DebrisType(11).Dimensions.depth_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(11).Dimensions.height_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(11).Dimensions.width_m), 'in'];

DebrisType(12).Name = 'WoodSheet2';
% DebrisType(12).Dimensions.depth_m = 0.00635; % 1/4 in
% DebrisType(12).Dimensions.height_m = 0.9144; % 36 in
% DebrisType(12).Dimensions.width_m = 0.9144;
DebrisType(12).Dimensions.depth_m = 0.009525; % 3/8 in
DebrisType(12).Dimensions.height_m = 0.05715; % 2.25 in
DebrisType(12).Dimensions.width_m = 0.05715;
DebrisType(12).Dimensions.density_kg_m3 = 500.0;
DebrisType(12).Dimensions.string = [num2str(100 / 2.54 * DebrisType(12).Dimensions.depth_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(12).Dimensions.height_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(12).Dimensions.width_m), 'in'];

DebrisType(13).Name = 'WoodSheet3';
% DebrisType(13).Dimensions.depth_m = 0.00635; % 1/4 in
% DebrisType(13).Dimensions.height_m = 0.6096; % 24 in
% DebrisType(13).Dimensions.width_m = 0.6096;
DebrisType(13).Dimensions.depth_m = 0.00635; % 1/2 in
DebrisType(13).Dimensions.height_m = 0.0762; % 3 in
DebrisType(13).Dimensions.width_m = 0.0762;
DebrisType(13).Dimensions.density_kg_m3 = 500.0;
DebrisType(13).Dimensions.string = [num2str(100 / 2.54 * DebrisType(13).Dimensions.depth_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(13).Dimensions.height_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(13).Dimensions.width_m), 'in'];

DebrisType(14).Name = 'WoodSheet4';
% DebrisType(14).Dimensions.depth_m = 0.00635; % 1/4 in
% DebrisType(14).Dimensions.height_m = 0.3048; % 12 in
% DebrisType(14).Dimensions.width_m = 0.3048;
DebrisType(14).Dimensions.depth_m = 0.0254; % 1 in
DebrisType(14).Dimensions.height_m = 0.1524; % 12 in
DebrisType(14).Dimensions.width_m = 0.1524;
DebrisType(14).Dimensions.density_kg_m3 = 500.0;
DebrisType(14).Dimensions.string = [num2str(100 / 2.54 * DebrisType(14).Dimensions.depth_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(14).Dimensions.height_m), 'in x ',...
    num2str(100 / 2.54 * DebrisType(14).Dimensions.width_m), 'in'];


save('~/Documents/code/misc/simradar-RCS-data.mat', 'DebrisType');

