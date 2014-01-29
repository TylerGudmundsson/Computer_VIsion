function DrunkWalk
% generates multiple butterflies using function call that move randomly in
% all directions with 4 fixed lion in 4 corners. If a butterfly reaches
% withing some minimum distance from a Lion's head, a beep sound will be
% heard, the butterfly will disappear and a new one will appear in the
% middle of the window.
% Note: butterflies move in the direction their bodies are aligned
%
%
global LargeImage guy guy_template  %will be used by a function
OrigLargeImage=imread('Tulips.jpg','jpg');
car=imread('Car.jpg','jpg');
car=imresize(car,0.25);
template_car=imread('template_car.jpg');

[m1,n1,k1]=size(car);

% Now create a Look Up table that contain the coordinates of all pixels
% that comprise the butterfly
% embede 4 lion images at 4 corners
%OrigLargeImage(100:100+m1-1,800:800+n1-1,:)=car;
%OrigLargeImage(550:550+m1-1,800:800+n1-1,:)=car;
OrigLargeImage(100:100+m1-1,100:100+n1-1,:)=car(:,n1:-1:1,:);
%OrigLargeImage(550:550+m1-1,100:100+n1-1,:)=car(:,n1:-1:1,:);
% create a small cross in the mouth of the upper left lion
OrigLargeImage(136,215:231,1)=255;
OrigLargeImage(128:144,223,1)=255;

% display the image for debugging purpose
figure(1);imagesc(OrigLargeImage);

% read both the butterfly image and its template, reduce them by a factor
% of 0.5 then extract its dimensions
guy=imread('DRUNKGUY.jpg','jpg');
guy=imresize(car,0.5);
guy_template=imread('DRUNKGUYTEMPLATE.jpg','jpg');
guy_template=imresize(car_template,.5);
[m,n,k]=size(guy_template);

% Now create a Look Up table that contain the coordinates of all pixels
% that comprise the butterfly
Lu=[]; %look up table for butterfly coordinates. Don't know its size
k=0;
for i=1:m,
    for j=1:n,
        if ((guy_template(i,j,:))< (230) ),
            k=k+1;
            Lu(k,:)=[i,j];
        end
    end
end



[M,N,K]=size(OrigLargeImage);
[m,n,k]=size(guy);

NumOfGuy=2; % number of butterflies to be embedded
Guy_Position=400*ones(NumOfGuy,2); %initial positions: coordinates (400,400)

% Main loop of the program
for kk=1:50, % run the loop so many times
    LargeImage=OrigLargeImage; % make a clean copy of the original image
        % create a new position for each butterfly
        for i=1:size(Guy_Position,1),
            rr=ceil(99*rand(1,2))-50;
            x=Guy_Position(i,1)+rr(1);
            y=Guy_Position(i,2)+rr(2);
        % test the new position is not too close to borders
        while ( x<100 | x >(M-100) | y<100 | y>(N-100) ),
             rr=ceil(59*rand(1,2))-30;
             x=Guy_Position(i,1)+rr(1);
             y=Guy_Position(i,2)+rr(2);
        end
        % call function to embede butterfly (see function EmbedeObject)
        EmbedeObject(guy,Lu,Guy_Position(i,:),rr)
        % update butterfly position
        Guy_Position(i,1)=x;
        Guy_Position(i,2)=y;
        % ----- test for collision with lions ----
        L1=[136 224];L2=[590 225]; L3=[136 836]; L4=[586 836]; %lion loc's
        % distance of butterfly from each lion
        D1=L1-Guy_Position(i,:);
        D2=L2-Guy_Position(i,:);
        D3=L3-Guy_Position(i,:);
        D4=L4-Guy_Position(i,:);
        RR=30;  % max allowed distance
        % if butterfly too close to a lion, create a sound, change its
        % position (equivalent to removing it from its location and
        % creating a new one)
        if ( abs(D1)<RR || abs(D2)<RR || abs(D3)<RR || abs(D4)<RR ),
             a=(500:2000);
            soundsc(sin(a));
            Guy_Position(i,:)=[400 500];
            pause(.03); % pause needed for soundsc function
        end
         
        
    end
    
    % now display the composed image. This figure will be display every
    % cycle of the main loop
    figure(1);imagesc(LargeImage);   
    pause(0.1)% needed for display to work properly 
 
end % of main loop

% ------------------------------------------------------------------------
% This function will embede any image with its LU table in the larger
% image, LargeImage
%   Pic image to be embeded (need to be smaller than LargeImage
%   Lu  its look up table (contains coordinates of object in image)
%   PositionVector contains the coord's (row, col)of the current position
%   rr  a row vector of 1x2 containing the displacement of the butterfly
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