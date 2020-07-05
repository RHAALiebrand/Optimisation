function [F]=import_database()
% Import parameters
Params;

% Read LHS matrix
% r,h,c 
LHD=dlmread('Cd_database/Data/n=3_p=2000.dat'); % Is not necessary for calculation but import for plotting 
LHD_values=dlmread('Cd_database/Data/n=3_p=2000_values.dat');
LHD_size=size(LHD);  
samples=LHD_size(1);  % Number of samples
%scatter3(LHD(:,1),LHD(:,2),LHD(:,3))

% Import database
Cd_data=zeros(samples,1);
for i=[1:samples]
    % Sometimes Xfoil does not converge, therefore try to do this and
    % otherwise just doe nothing
    try
    xfoil_output=dlmread(strcat('Cd_database/Data/c',num2str(LHD_values(i,3),6),'_h',num2str(LHD_values(i,2),6),'_r',num2str(LHD_values(i,1),6),'.dat'),' ',12,0);
    Cd_index=find(xfoil_output ~= 0);
    Cd_data(i)=xfoil_output(Cd_index(1));
    end 
end

% Now check where XFOIL is not converged and remove corresponding Cd + corresponding datapoint (r,h,c) 
index_non_conv=find(Cd_data<=0.001);
Cd_data(index_non_conv)=[];
LHD_values(index_non_conv,:)=[];

% Making interpolantion surface
F=scatteredInterpolant(LHD_values(:,1),LHD_values(:,2),LHD_values(:,3),Cd_data,'linear','linear');

%% Plot LHD
% figure()
% scatter3(LHD_values(:,1),LHD_values(:,2),LHD_values(:,3))
% xlabel('$r/c$ [-]','fontsize',16,'Interpreter','LaTex')
% ylabel('$h/c$ [-]','fontsize',16,'Interpreter','LaTex')
% zlabel('$c$ [mm]','fontsize',16,'Interpreter','LaTex')
% ylim([0.1 0.6])
% Determine drag force


end 