function [A,Rsq]=f_Gauss2D_contextual(S,Plot,FitOrientation,outputDir,id,row_angle ,col_angle,col_width,row_width,Boundaries,stim_circle)
% Fit a 2D Gaussian Function to Data
% INPUT:
%   S: two-dimensional array of size nxm.

% OUTPUT: 
%	A = [Amplitude, x0, Sigma_x, y0, Sigma_y, Angle(in Radians)]
%                              Residual sum of squares

% 1. 2D Gaussian function ( A requires 5 coefs ).
g = @(A,X) A(1)*exp( -((X(:,:,1)-A(2)).^2/(2*A(3)^2) + (X(:,:,2)-A(4)).^2/(2*A(5)^2)) );
% 2. 2D Rotated Gaussian function ( A requires 6 coefs ).
f = @(A,X) A(1)*exp( -(...
    ( X(:,:,1)*cos(A(6))-X(:,:,2)*sin(A(6)) - A(2)*cos(A(6))+A(4)*sin(A(6)) ).^2/(2*A(3)^2) + ... 
    ( X(:,:,1)*sin(A(6))+X(:,:,2)*cos(A(6)) - A(2)*sin(A(6))-A(4)*cos(A(6)) ).^2/(2*A(5)^2) ) );


%% ---Parameters---
[n,m]=size(S);      % n x m pixels area/data matrix
[~,l]=max(S(:));
[row,col]=ind2sub(size(S),l);

A0 = [1,col,m/5,row,n/5,0];   % Inital (guess) parameters
InterpMethod='nearest'; % 'nearest','linear','spline','cubic'
% Numerical Grid
[x,y]=meshgrid(1:m,1:n); X=zeros(n,m,2); X(:,:,1)=x; X(:,:,2)=y;

% Define lower and upper bounds [Amp,xo,wx,yo,wy,fi]
lb = [0,-1,0,-1,0,0];
ub = [realmax('double'),m+1,m/3,n+1,n/3,pi/4];

% Fit sample data
switch FitOrientation
    case 0, [A,RSS] = lsqcurvefit(g,A0(1:5),X,S,lb(1:5),ub(1:5));A=[A,0];
    case 1, [A,RSS] = lsqcurvefit(f,A0,X,S,lb,ub);
    otherwise, error('invalid entry for Fit Orientation');
end

Rsq=1-(RSS/(rssq(S(:)-mean(S(:)))^2));
if Plot

    %% half width at half maximum
    HWHM_X=sqrt(2*log(2))*A(3);
    HWHM_Y=sqrt(2*log(2))*A(5);
    HWHM_Y_angle = HWHM_Y*row_angle;
    HWHM_X_angle = HWHM_X*col_angle;
    RF_size = (HWHM_Y_angle+HWHM_X_angle)/2;
    %% ellips
    t = linspace(0,2*pi) ;
    ellips_x = (A(2)+ HWHM_X*cos(t))*col_width ;
    ellips_y = (A(4)+ HWHM_Y*sin(t))*row_width ;
    % rotationAngle = -2*A(6);
    % transformMatrix = [cos(rotationAngle), sin(rotationAngle);...
    % 	-sin(rotationAngle), cos(rotationAngle)];
    % xyAligned = [ FWHM_X*cos(t); FWHM_Y*sin(t)]';
    % xyRotated = xyAligned * transformMatrix;
    % xRotated = xyRotated(:, 1)+A(2);
    % yRotated = xyRotated(:, 2)+A(4);
    
    
    %% Plot 3D Data and Fitted curve
%     hf1=figure('visible','off'); set(hf1,'Position',[1000 600 800 500]);
%     switch FitOrientation
%         case 0, C=del2(g(A,X)); mesh(x,y,g(A,X),C); hold on
%         case 1,  C=del2(f(A,X)); mesh(x,y,f(A,X),C); hold on
%     end
%     surface(x,y,S,'EdgeColor','none'); alpha(0.5);
% %     view(-60,20);
%     grid on; hold off
%     title(['RF map, Neuron:',num2str(id),', RF radius:',num2str(RF_size),' deg, RSS:',num2str(RSS)])
%     saveas(hf1,fullfile(outputDir,['RF_3D_cell-',num2str(id),'.png']))
%     savefig(hf1,fullfile(outputDir,['RF_3D_cell-',num2str(id),'.fig']))
    %% Plot Sample Pixels data
    hf2=figure('visible','off'); set(hf2,'Position',[20 20 1000 800]);
    subplot(4,4,[5,6,7,9,10,11,13,14,15]); imagesc(x(1,:)*col_width,y(:,1)*row_width,S);
    hold on;
%     plot(ellips_x,ellips_y,'k', 'LineWidth', 2);
    plot(Boundaries(:,2),Boundaries(:,1),'r','LineWidth',3);
    plot(stim_circle)
    hold off;

    % Output and compare data and fitted function coefs
    text(5,(n+(n/11))*row_width,sprintf('Amplitude \t \t  X_0 \t \t \t sigma_X \t \t \t Y_0 \t \t \t \t \t sigma_Y \t \t Angle'),'Color','blue');
    text(5,(n+(n/8))*row_width,sprintf(' %1.3f \t \t  %1.3f \t \t \t %1.3f \t \t \t %1.3f \t \t \t %1.3f \t \t \t %1.3f',A),'Color','red');
    % Plot vertical and horizontal axis
    vx_h=x(1,:); vy_v=y(:,1);
    switch FitOrientation
        case 1, M=-tan(A(6));
            % generate points along _horizontal & _vertical axis
            vy_h=M*(vx_h-A(2))+A(4); hPoints = interp2(x,y,S,vx_h,vy_h,InterpMethod);
            vx_v=M*(A(4)-vy_v)+A(2); vPoints = interp2(x,y,S,vx_v,vy_v,InterpMethod);
        case 0, A(6)=0;
            % generate points along _horizontal & _vertical axis
            vy_h=A(4)*ones(size(vx_h)); hPoints = interp2(x,y,S,vx_h,vy_h,InterpMethod);
            vx_v=A(2)*ones(size(vy_v)); vPoints = interp2(x,y,S,vx_v,vy_v,InterpMethod);
    end
    % plot lines
%     hold on; plot(A(2),A(4),'+b',vx_h,vy_h,'.r',vx_v,vy_v,'.g'); hold off;
    % Plot cross sections
    dmin=1.1*min(S(:)); xfit=x(1,:); hfit=A(1)*exp(-(xfit-A(2)).^2/(2*A(3)^2));
    dmax=1.1*max(S(:)); yfit=y(:,1); vfit=A(1)*exp(-(yfit-A(4)).^2/(2*A(5)^2));
    subplot(4,4,[1,2,3]); xposh = (vx_h-A(2))/cos(A(6))+A(2);
    plot(xposh,hPoints,'r.',xfit,hfit,'black'); grid on; axis([1,m,dmin,dmax]);xticks([]);
    subplot(4,4,[8,12,16]); xposv = (vy_v-A(4))/cos(A(6))+A(4);
    plot(vPoints,xposv,'g.',vfit,yfit,'black'); grid on; axis([dmin,dmax,1,n]);yticks([]);
    set(gca,'YDir','reverse');
    sgtitle(['RF map, Neuron:',num2str(id),', RF radius:',num2str(RF_size),' deg, Rsq:',num2str(Rsq)])
    saveas(hf2,fullfile(outputDir,['RF_cell-',num2str(id),'.png']))
    
    close all
    
end