%% Engineering optimization 
% Structural model
% Martin Janssens
% Rens Liebrand

function [x_coord,y_coord]=Airfoil_coordinates(r,h,c)
NACAcode='0012';
wd = fileparts(which(mfilename));
fname = mfilename;
file = fopen([wd filesep 'airfoilcommands.inp'],'w'); 

% XFOIL commands
fprintf(file,'NACA %s\n',NACAcode); 
fprintf(file,'GDES\n');
fprintf(file,'TSET\n');
fprintf(file,strcat(num2str(h),'\n\n'));
fprintf(file,'HIGH\n');
fprintf(file,strcat(num2str(r),'\n\n'));
fprintf(file,'\n'); 
fprintf(file,'PCOP\n'); 
fprintf(file,'PPAR\n'); 
fprintf(file,'N\n'); 
fprintf(file,'20\n\n\n'); 
fprintf(file,'save %s\n','xfoilprofile.txt');
fprintf(file,'\nquit\n');
fclose(file);

% run XFOIL
cmd = sprintf('cd %s && xfoil.exe < airfoilcommands.inp > xfoil.out',wd);
system(cmd);

datapoints = importdata('xfoilprofile.txt',' ',1);
datapoints = datapoints.data;

x_coord = datapoints(:,1);
y_coord = datapoints(:,2);

delete('xfoilprofile.txt')
delete('airfoilcommands.inp')

x_coord=x_coord*c;
y_coord=y_coord*c;


end