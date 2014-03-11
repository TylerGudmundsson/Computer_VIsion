function CarMovement
global x y 
OrigLargeImage=imread('BirdsEye.jpg','jpg');
car=imread('Car.jpg','jpg');
car=imresize(car,.1);
car_temp=imread('template_car.jpg');
car_temp=imresize(car_temp,.1);
car_temp=imrotate(car_temp,100);
car=imrotate(car,100);
[m1,n1,k1]=size(car);
rr=[1,2];
Lu=[]; %look up table for butterfly coordinates. Don't know its size
lucount=0;
for i=1:m1,
    for j=1:n1,
        if ((car_temp(i,j,:))> (230) ),
            lucount=lucount+1;
            Lu(lucount,:)=[i,j];
        end
    end
end
x=0;
Car_p=[1:1,1:2];
for t=1:20
  x=x+5;
  y=-2.18*(x+200) + 1184;
  %y1=-2.18*((x-5)+100) + 1184;
% call function to embede butterfly (see function EmbedeObject)
EmbedeObject(car,Lu,Car_p(1,:),rr)
        % update butterfly position
        %rr=((ceil(5)) (ceil(-2.18*5)));
        rr=[5 ceil(2.18*5)];
        Car_p(1,1)=y;
        Car_p(1,2)=x;
     
end
    
end       
        
  function EmbedeObject(Pic,Lu,PositionVector,rr)
global LargeImage
rr=[5 -ceil(2.18*5)];
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
end      