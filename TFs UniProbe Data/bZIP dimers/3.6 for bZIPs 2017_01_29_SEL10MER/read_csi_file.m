% FUNCTION TO READ INPUT FILE

function  [F,R,Int,ls,n]=read_csi_file(inp,colum)

% IF FRMT.MAT EXISTS IT READS COLUMN 1 OR 2 ACCORDING TO VALUE OF frmt
% IT TAKES CARE OF 1 HEADER LINE BY ITSELF
if exist('frmt.mat')
    load frmt;
    if frmt==1
        [F,R,Int]=textread(inp,'%s	%s	%s	%*s	%*s	%*s	%*s	%*s',1);
        if isempty(str2num(char(Int(1))))
            [F,R,Int]=textread(inp,'%s	%s	%f	%*s	%*s	%*s	%*s	%*s','headerlines',1);
        else
            [F,R,Int]=textread(inp,'%s	%s	%f	%*s	%*s	%*s	%*s	%*s');
        end
    else 
        [F,R,Int]=textread(inp,'%s	%s	%*s	%s	%*s	%*s	%*s	%*s	%*s',1);
        if isempty(str2num(char(Int(1))))
            [F,R,Int]=textread(inp,'%s	%s	%*f	%f	%*s	%*s	%*s	%*s	%*s','headerlines',1);
        else
            [F,R,Int]=textread(inp,'%s	%s	%*f	%f	%*s	%*s	%*s	%*s	%*s');
        end 
    end
else
    % IF FRMT.MAT IS NOT THERE IT READS COLUMN NUMBER colum
    xx='[F,R,Int]=textread(inp,''%s %s';
    for i=2:colum
        xx=[xx ' %*s'];
    end
    xx=[xx ' %s'];
	for i=1:8
        xx=[xx ' %*s'];
    end
    xx=[xx ' '',1); '];
    eval(xx);
    
    if isempty(str2num(char(Int(1)))) % IF THERE IS A HEADER LINE
        xx='[F,R,Int]=textread(inp,''%s %s';
        for i=2:colum
            xx=[xx ' %*s'];
        end
            xx=[xx ' %f'];
        for i=1:8
            xx=[xx ' %*s'];
        end
        xx=[xx ' '',''headerlines'',1); '];
        eval(xx);
    else   % IF THERE IS NO HEADER LINE
        xx='[F,R,Int]=textread(inp,''%s %s';
        for i=2:colum
            xx=[xx ' %*s'];
        end
            xx=[xx ' %f'];
        for i=1:8
            xx=[xx ' %*s'];
        end
        xx=[xx ' ''); '];
        eval(xx);
    end
end

F=upper(char(F));
R=upper(char(R));

n=size(F,1);
ls=size(F,2);
ls2=size(R,2);

if ls~=ls2 || ~isempty(find(R==' ',1)) || ~isempty(find(F==' ',1)) % IF LENGTH OF FORWARD AND REVERSE SEQUENCE NOT EQUAL
    errordlg('Error in the format of input file','Bad Input File','modal');
end

end