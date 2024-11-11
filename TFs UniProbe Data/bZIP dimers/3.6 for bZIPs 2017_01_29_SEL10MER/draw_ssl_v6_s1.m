% PROGRAM TO GENERATE CIRCULAR LANDSCAPE

function [X,Y,Z,f]= draw_ssl_v6_s1(Uint3,Us3,lm,sm,l_l,op,n,l,ln,ls2,ring_start,ring_stop,mx)

f=figure;
hold on;
set(f,'Renderer','zbuffer');

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

% WRITING TO OUTPUT FILE
fid = fopen(op,'w');
fprintf(fid,'X\tY\tIntensity(Z)\tSmoothInt(Z)\tSequence\tSequenceReverse\t mismatch\n');

% PLOT EACH MISMATCH RING
for i=im:ie
    if i==1 || l(i)>=10
        
        % DEFINING COORDINATES TO BE PLOTTED
        theta=linspace(pi,-pi,l(i));
        r=[10*i-l_l; 10*i; 10*i+l_l];
        X2=r*cos(theta);
        Y2=r*sin(theta);
        Z2=[0 ;1 ;0]*Uint3(ln(i,1):ln(i,2))';
        Zs2=Z2;
        
        % SMOOTHEN THE Z-SCORE
        Zs2(2,:)=smooth(Z2(2,:),round((l(i)*sm/100)+1.5));
        
        % PLOT i-th RING
        surf(X2,Y2,Zs2);
        
        X(ln(i,1):ln(i,2))=X2(2,:);
        Y(ln(i,1):ln(i,2))=Y2(2,:);
        Z(ln(i,1):ln(i,2))=Z2(2,:);
        Zs(ln(i,1):ln(i,2))=Zs2(2,:);
        
        for j=ln(i,1):ln(i,2)
            % WRITE ALL THE POINTS IN A TEXT FILE
            fprintf(fid,'%8.4f\t %8.4f\t %8.4f\t %8.4f\t %s\t %s\t %d\n',...
                X(j),Y(j),Z(j),Zs(j),char(Us3(j,lm:ls2-lm+1)),seqrcomplement(char(Us3(j,lm:ls2-lm+1))),(i-1));
        end
    end
end
fclose(fid);

% CHANGE SHADING
shading interp;

% SETTING AXIS Ticks AND LABELs
set(gca, 'XTickLabelMode', 'Manual')
set(gca, 'XTick', [])
set(gca, 'YTickLabelMode', 'Manual')
set(gca, 'YTick', [])
set(gca, 'ZTickLabelMode', 'Manual')
set(gca, 'ZTick', [])

% CHANGING THE VIEW ANGLE
view(90,70);
% view(3);

% CHANGE IN COLORMAP TO ASSIGN MAX TO mx
colormp2(mx);

hold off;
datacursormode on;
