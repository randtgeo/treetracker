function [Results] = RT_Convert(PointCloud,RT)

% Read in point cloud and calculate Center Point 
 
% [filename,path]=uigetfile('TT4_2020-05.txt'); 
%     PointCloud=dlmread(filename) ;
    
    CP= [median(PointCloud(:,1));median(PointCloud(:,2));median(PointCloud(:,3));1];
%for i=1:length(file)
% Read in text file of RT matrix:
%     [filename,path]=uigetfile('TT4_2020-05_to_2020-10.txt'); 
%     RT=dlmread(filename,'',2,0); % Use line 9 if header was not deleted
    %RT=dlmread(filename); % Use this one if header was deleted 
 
% User Input of center Point: 
 
    %prompt='What is the center point (format [X;Y;Z;1])'; % Use to prompt for center point 
    %CP=input(prompt);
   % CP=[14.8859;479.7249;39.2344;1]; % Enter CP into script here, or use
   % prompt
 
%Translation vector between initial CP and final CP
 
    T=RT*CP-CP;
    Tlength=(T(1)^2+T(2)^2+T(3)^2)^0.5; % calculates the length of the translation vector from the center point to new center point
   % CPnew=RT*CP % new center point used for scan to scan analysis 
 
%Trend Translation vector (T) 
    if T(1)>=0 && T(2)>=0  
        Trend=atan(T(1)/T(2))*180/pi; % calculates trend of translation if vector is in Quandrant I - +y being North (0-90 degrees) 
 
    elseif T(1)<= 0&& T(2)>=0 
        Trend=atan(T(1)/T(2))*180/pi+360; % calculates trend if translation vector is in Quandrant II (270-360)
    
    else
        Trend=atan(T(1)/T(2))*180/pi+180; % calculates trend if translation vector is in Quandrant III or IV (90-180 & 180-270)
 
    end 
 
% Calculate Plunge of Displacement:
         
         PlungeT=asin(dot(T,[0;0;1;0])/(norm(T)*norm([0;0;1;0])))*180/pi;
          
% Topple Angle Conversion:  
 
    Mrot=[RT(:,1),RT(:,2),RT(:,3),[0;0;0;1]]; % Align initial state scan to sucessive scan
     
    vtop=Mrot*[0;0;1;0]; % Calculates toppling vector from rotation matrix 
     
    Topangle=acos(dot(vtop,[0;0;1;0])/(norm(vtop)*norm([0;0;1;0])))*180/pi; % divide the dot product of the two vectors (vtop and vertical) by the products of magnitudes - gives the cosine angle between the two vectors [AdotB=magA*magBcosAng]
     
 % Azimuth of topple vector:
 
    if vtop(1)>=0 && vtop(2)>=0  
         
        Aztop=atan(vtop(1)/vtop(2))*180/pi;  % calculates azimuth if topple vector is in Quandrant I - +y being north 
     
    elseif vtop(1)<= 0&&vtop(2)>=0 
         
        Aztop=atan(vtop(1)/vtop(2))*180/pi+360; % calculates azimuth if topple vector is in Quandrant II
    
    else
        Aztop =atan(vtop(1)/vtop(2))*180/pi+180;% calculates azimuth if topple vector is in Quandrant III or IV         
    end 
 
% Tilt Angle Conversion 
 
    alpha=-asin(vtop(2)/((vtop(2)^2+vtop(3)^2)^0.5));% Alpha angle required for topple matrix calculation 
 
    Beta=asin(vtop(1)/((vtop(1)^2+vtop(3)^2)^0.5));% Beta angle required for topple matrix calculation 
 
    Mtopple=[cos(Beta),0,-sin(Beta),0;0,1,0,0;sin(Beta),0,cos(Beta),0;0,0,0,1]*[1,0,0,0;0,cos(alpha),sin(alpha),0;0,-sin(alpha),cos(alpha),0;0,0,0,1];% Topple matrix  
 
    Mtilt=inv(Mtopple)*Mrot; % Tilt matrix is the product of the inverse topple matrix with the rotation matrix 
 
    xt=vtop(1); % x component of toppling vector 
 
    yt=vtop(2); % y component of toppling vector 
 
    zt=vtop(3); % z component of toppling vector 
 
    % F1 to F9 correspond to the differences of each element of the two
    % Mtilt matrices 
 
    f1=@(Kappa)  (1+(1-cos(Kappa))*((xt^2)-1)- Mtilt(1,1));
 
    f2=@(Kappa) (-zt*sin(Kappa)+(1-cos(Kappa))*xt*yt-Mtilt(1,2));
 
    f3=@(Kappa)  (yt*sin(Kappa)+(1-cos(Kappa))*xt*zt-Mtilt(1,3));
 
    f4=@(Kappa)  (zt*sin(Kappa)+(1-cos(Kappa))*xt*yt-Mtilt(2,1));
 
    f5=@(Kappa)  (1+(1-cos(Kappa))*((yt^2)-1)-Mtilt(2,2));
 
    f6=@(Kappa)  (-xt*sin(Kappa)+(1-cos(Kappa))*yt*zt-Mtilt(2,3)); 
 
    f7=@(Kappa)  (-yt*sin(Kappa)+(1-cos(Kappa))*xt*zt-Mtilt(3,1));
 
    f8=@(Kappa)  (xt*sin(Kappa)+(1-cos(Kappa))*yt*zt-Mtilt(3,2));
 
    f9=@(Kappa)  (1+(1-cos(Kappa))*((zt^2)-1)-Mtilt(3,3));
 
    f=@(Kappa) f1(Kappa).^2+f2(Kappa).^2+f3(Kappa).^2+f4(Kappa).^2+f5(Kappa).^2+f6(Kappa).^2+f7(Kappa).^2+f8(Kappa).^2+f9(Kappa).^2; % squared difference function for least squares minimization 
 
    options=optimset('TolX',1e-8);
    Tiltangle=fminsearch(f,0,options)*180/pi; % search for min around 0 
 
    %Plot tilt residual function 
        %Kappa=-2*pi:pi/1000:2*pi; % plots function radians on x axis. 
        %plot(Kappa,f(Kappa))
        %grid on
         
      ErrorX=T(1);
      ErrorY=T(2);
      ErrorZ=T(3);
       
      if atan(0.005/ErrorX)>atan(0.005/ErrorY)
      TrendErr= atan(0.005/ErrorX*180/pi());
      else
           TrendErr=atan(0.005/ErrorY)*180/pi();
      end 
       
     PlungeErr= atan(0.005/ErrorZ)*180/pi();
       
     ToppAzErr=asin(0.016/Topangle)*180/pi() ;
      
      %ErrorZ(i)=T(3);
      Results = [T(1:3)',Tlength,PlungeT,Trend,Topangle,Aztop,Tiltangle];
   % Results=[Tlength,PlungeT,PlungeErr,Trend,TrendErr,Topangle,Aztop,ToppAzErr,Tiltangle]
    %Reseltsreal(i,:)=real(Results(i,:))
 
end
