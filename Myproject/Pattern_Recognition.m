function Pattern_Recognition
global pattern image Z X Y
pattern=imread('ball.jpg');
image=imread('Background.jpg');
%Scan images into global variables
[m,n,k]=size(pattern);
image=imrotate(image,180);
[M,N,K]=size(image);
figure1;imagsc(image);
%set up variables used through out program.

%failed attempt to use Lu Table.
%Lu=[];
%lucount=1;
%for a=1:m
   % for b=1:n
     %   if ((pattern(m,n)) > (230)),
      %      lucount=lucount+1;
      %      Lu(lucount,:)=[a,b];
     %   end
 %   end
%end

 

r=1;%Used to represent total number of possible locations for object
x=[];%Table of possible X coordinates
y=[];%Table of possible Y coordinates
tic
for t=1:(M-m)
     for c=1:(N-n)
        for i=1:(m)
            for j=1:(n)
 %Loop to compare pattern to background pixel by pixel. 
 %Move over a pixel in the background and repeat
        
 Z(i,j,1)=sum(sum(abs(image((i+t),(j+c))-(pattern(i,j)))));
  S=sum(sum(Z));

 
            end
        end
     %Records location of center of matching comparisons
        if S<4500,
     x(r,1)=t+(m/2);
     y(r,1)=c+(m/2);
     r=r+1;     
        end
    end
end
toc
X=sum(sum(x)); %Sum all possible X locations of ball
Y=sum(sum(y)); %Sum all possible Y locations of ball
X=X/(r-1); %Divide by number of possible x locations to find average
Y=Y/(r-1); %Divide by number of possible y locations to find average

image(ceil(X),:,1)=255; %Display crosshair at center of the suspected ball location
image(:,ceil(Y),1)=255;
figure(1);imagesc(image);
k=1;
end

    
    
    
    
    