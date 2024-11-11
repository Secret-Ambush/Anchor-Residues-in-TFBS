% PROGRAM TO GENERATE CIRCULAR LANDSCAPE

function [X,Y,Z,f]= draw_ssl_v6(Uint3,Us3,lm,sm,l_l,op,n,l,ln,ls2,ring_start,ring_stop,n_pt,mx)

fid = fopen(op,'w');
fprintf(fid,'X\t\t\t\t Y\t\t Intensity(Z)\tSmoothInt(Z)\tSequence\t mismatch\n');
f=figure;
hold on;
%colorbar;
 set(f,'Renderer','zbuffer');
%set(f,'OuterPosition',[500 200 500 700])
% grid on;

% KEEPING TRACK OF X,Y COORDINATES AND UN-SMOOTHENED & SMOOTHENED INTENSITY
X=zeros(n,1);
Y=zeros(n,1);
Z=zeros(n,1);
Zs=zeros(n,1);

% SELECT WHICH RINGS TO DISPLAY
if ring_start==0
    im=1;
else
    im=ring_start;
end
if ring_stop==0
    ie=lm+1;
else
    ie=ring_stop;
end

if (ie>3)
    ie=3;
end

% FOR EACH RING PLOT IT
for i=im:ie
if i==1 || l(i)>=10
    theta=linspace(pi,-pi,l(i));
    r=[10*i-l_l; 10*i; 10*i+l_l];
    X2=r*cos(theta);
    Y2=r*sin(theta);
    Z2=[0 ;1 ;0]*Uint3(ln(i,1):ln(i,2))';
    %Z2=[0 ;1 ;0]*(Uint3(ln(i,1):ln(i,2))-1)'; FOR CSI-SEQ
    Zs2=Z2;
    
    % SMOOTHEN THE INTENSITY
%     Zs2(2,:)=smooth(Z2(2,:),round((l(i)*sm/100)+1.5),'loess');
    Zs2(2,:)=smooth(Z2(2,:),round((l(i)*sm/100)+1.5));
    
    surf(X2,Y2,Zs2);
    
    X(ln(i,1):ln(i,2))=X2(2,:);
    Y(ln(i,1):ln(i,2))=Y2(2,:);
    Z(ln(i,1):ln(i,2))=Z2(2,:);
    
    Zs(ln(i,1):ln(i,2))=Zs2(2,:);
    
    if n_pt~=0
     tt=ceil(l(i)/n_pt);
    end
    
    
    % WRITE ALL THE POINTS IN A TEXT FILE
    for j=ln(i,1):ln(i,2)
       fprintf(fid,'%8.4f\t %8.4f\t %8.4f\t %8.4f\t %s\t %d\n',X(j),Y(j),Z(j),Zs(j),char(Us3(j,lm:ls2-lm+1)),(i-1));
%        if  (~isempty(strfind(char(Us3(j,lm:ls2-lm+1)),'AGATTA')) || ~isempty(strfind(char(Us3(j,lm:ls2-lm+1)),'CGATTA')) || ~isempty(strfind(char(Us3(j,lm:ls2-lm+1)),'GGATTA')))
%            text(X(j),Y(j),Z(j),char(Us3(j,lm:ls2-lm+1)),'FontSize',15);
%        end
       if n_pt~=0 % TO WRITE n_pt SEQUENCES ON THE RINGS
           if (rem(j-ln(i,1),tt)==0)
               text(X(j),Y(j),Z(j),char(Us3(j,lm:ls2-lm+1)),'FontSize',15);
           end    
       end
    end
end
end
fclose(fid);
shading interp;
set(gca, 'XTickLabelMode', 'Manual')
set(gca, 'XTick', [])
set(gca, 'YTickLabelMode', 'Manual')
set(gca, 'YTick', [])
set(gca, 'ZTickLabelMode', 'Manual')
set(gca, 'ZTick', [])
% zlabel('Z(Binding Intensity/Z-Score)');
% view(90,90);
view(90,70);
% view(3);
colormp2(mx);
hold off;
datacursormode on;
%dcm_obj = datacursormode(f);
%set(dcm_obj,'DisplayStyle','window');