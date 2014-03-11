function EdgeDetection
%This file applies the function imcontour and edge detection 
%using a few different methods.
%************************************************************************
% first read and display original image
I=imread('track.tif');
figure(1);imshow(I)
title('original image')
%Edge detection using 'prewitt' method
I2 = edge(I,'prewitt');
figure(3), imshow(I2)
%************************************************************************
%Creat Look Up Table for edge
[m,n,k]=size(I2);
Lu1=[];
k=0;
for i=1:m,
    for j=1:n,
        if ((I2(i,j,:))> (230) ),
            k=k+1;
            Lu1(k,:)=[i,j];
        end
    end
end


function EmbedeObject(Pic,Lu,PositionVector,rr)
global LargeImage
    % use the displacement values to compute the angle of displacement
    % then rotate Look Up table
    % first calculate angle teta from rr
    if rr(2)==0,
        teta=90;
    else
        teta=atand(rr(1)/rr(2));
    end
    if rr(1)<0,
        teta=teta+180;
    end
    % now find transformation matrix
    T=[cosd(teta) -sind(teta);sind(teta) cosd(teta)];
    % first translate to cartez coordinates
    [m,n,k]=size(Pic);
    R_Lu=size(Lu,1); %equivalent to number of pixels to be copied
    Lu_T(:,1)=Lu(:,1)-m/2; % m/2 is half of horiz  size of but
    Lu_T(:,2)=Lu(:,2)-n/2; % n/2 is half of vert  size of but
    b=T*Lu_T';
    Lu_T=ceil(b');
    % translate back to MATLAB coordinate
    Lu_T(:,1)=Lu_T(:,1)+m/2; % m/2 is half of horiz  size of but
    Lu_T(:,2)=n/2 -Lu_T(:,2); % n/2 is half of vert  size of but
    M1=round(mean(Lu_T(:,1))); % compute the center of the new LU
    M2=round(mean(Lu_T(:,2)));
    
    % Copy each pixel from Pic (using new LU) to LargeImage
    % Subtract the mean of LU 
    for i=1:R_Lu,
        if ( PositionVector(1)>0 && PositionVector(2)>0)
            LargeImage(PositionVector(1)+Lu_T(i,1)- M1 ,PositionVector(2)+Lu_T(i,2)-M2,:)=Pic(Lu(i,1),Lu(i,2),:);
        end
    end
   
b=1;