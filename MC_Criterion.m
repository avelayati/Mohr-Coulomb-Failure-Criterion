%%
%%% Arian Velayati, PhD
%%%% This script is used to develop the Mohr-Cuolomb Failure envelope
   ... and determine rock's shear parameters

clc; clear; close;
%% Input (Principal Stresses)

S3 = [10;20;30;40]; % Confning pressure (MPa)
S1 = [132;183;218;259]; % Axial peak stress (MPa)

%% Plot S1 vs S3 from which parameters such as ucs, cohesion, and friction angle can be determined
% Linear fitting
P = polyfit(S3,S1,1);
% rsqr

% Method 1
Co = P(2); % UCS (MPa)
k = P(1); % Slope of the S1.S3 plot
Phi = 2*(atand(sqrt(k))-45); % Friction Angle
Teta = 45 + Phi/2; % Teta
Si = Co/(2*tand(Teta)); % Cohesion (MPa)

% Method 2
Phi2 = asind((k-1)/(k+1));
Si2 = Co*(1-sind(Phi2))/(2*cosd(Phi2));

%% plot S1 vs. S3; Use basic linear fitting to show the failure envelope
r = linspace(0,max(S3),100);
y = P(1)*r + P(2);
figure(1)
scatter(S3,S1,'r')
hold on
plot(r,y)
xlabel('Conf. Stress (MPa)')
ylabel('S1 (MPa)')
caption = sprintf('y = %.2f * x + %.2f', P(1), P(2));
text(5,240,caption, 'FontSize', 16, 'Color', 'r', 'FontWeight', 'bold');
hold off

%% MC Circle plot; Draw a line to form the failure envelope
D = (S1-S3); rad = D/2; % Radius of the circle
h = S3 + rad; % offset of the centre

n=1000;
h2 = zeros(length(S3),n);

for i =1:length(S3)
X(i,:) = linspace(S3(i),S1(i),n);
end
for i=1:n
    h2(:,i) = h;
    rad2(:,i) = rad;
end

   for i =1:length(S3)
Y(i,:) = sqrt(rad2(i,:).^2 - (X(i,:)-h2(i,:)).^2);
   end
   
   %Plot circles

figure(2)
   plot(X(1,:),Y(1,:),'r', 'LineWidth',2)
   hold on
   plot(X(2,:),Y(2,:),'b','LineWidth',2)
   plot(X(3,:),Y(3,:),'y','LineWidth',2)
   plot(X(4,:),Y(4,:),'g','LineWidth',2)
   xlabel('Sn (MPa)')
   ylabel('SS (MPa)')
   xx = 0:max(S1); yy = tand(Phi).*xx+Si;
   plot(xx,yy)
   %% Output
   T = table(Phi,Si,Co,Teta,'VariableNames',{'Friction_Angle','Cohesion_MPa','UCS_MPa','Teta'})
  
    writetable(T,'MC_Out.csv')
