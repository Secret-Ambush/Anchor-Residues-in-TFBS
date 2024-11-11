% FUNCTION TO GET MISMATCHES AND ARRANGEMENT FOR THE MOTIF SEQUENCE

function  [Int,Int_indx,Ucode1,Us1,l,ln,n,ls2]=get_mismatches_v6_perl_v2_s1(inp,motif,lm,ls,sel,colum)

perl('ssl_5.pl',inp,motif,num2str(ls),num2str(colum));
A=importdata(['perl_' motif '_' inp]);

n=length(A.rowheaders);
seq=['X'*ones(n,lm-1) char(A.rowheaders) 'X'*ones(n,lm-1)];
ls2=ls+2*(lm-1); % +4 for XX XX

dat=A.data; % size = n,lm+3

l=zeros(3,1); % NUMBER OF SEQUENCES IN EACH RING
ln=zeros(3,2);% START AND STOP OF SEQUENCES IN THE SORTED LIST FOR EACH RING
Us1=zeros(n,ls2);% KEEPS TRACH OF THE SEQUENCE
Ucode1=zeros(n,lm+2+ls2+lm+1); % ALL THE SEQUENCE INFORMATION IS STORED IN IT

for i=1:n
    k=dat(i,lm+3);
    p2=dat(i,2:lm+1);
    pos=dat(i,lm+2)+(lm-1)-2;
    
    % GENERATING Ucode
    Us1(i,:)=seq(i,:);
    XX=char('X'*ones(1,lm));
    YY=Us1(i,pos:lm+pos-1);
    XX(p2==1)=YY(p2==1);
    l(k+1)=l(k+1)+1;
    
    % DIFFERENT ARRANGEMENTS FOR SEQUENCES IN A RING
    switch sel
        case 1
            Ucode1(i,:)=[k p2(lm:-1:1) XX Us1(i,lm+pos:ls2) 'X' Us1(i,pos-1:-1:1) motif pos];
        case 2
            Ucode1(i,:)=[k p2(lm:-1:1) XX Us1(i,pos-1:-1:1) 'X' Us1(i,lm+pos:ls2) motif pos];
        case 3
            Ucode1(i,:)=[k p2(lm:-1:1) Us1(i,lm+pos:ls2) 'X' Us1(i,pos-1:-1:1) XX motif pos];
        case 4
            Ucode1(i,:)=[k p2(lm:-1:1) Us1(i,pos-1:-1:1) 'X' Us1(i,lm+pos:ls2) XX motif pos];
        case 5
            Ucode1(i,:)=[k p2(lm:-1:1) XX pos Us1(i,lm+pos:ls2) 'X' Us1(i,pos-1:-1:1) motif];
        case 6
            Ucode1(i,:)=[k p2(lm:-1:1) XX pos Us1(i,pos-1:-1:1) 'X' Us1(i,lm+pos:ls2) motif];
        case 7
            Ucode1(i,:)=[k p2(lm:-1:1) pos Us1(i,lm+pos:ls2) 'X' Us1(i,pos-1:-1:1) XX motif];
        case 8
            Ucode1(i,:)=[k p2(lm:-1:1) pos Us1(i,pos-1:-1:1) 'X' Us1(i,lm+pos:ls2) XX motif];
        otherwise
            error('Incorrect sel input should be between 1-8');
    end
    
end

%FINDING NUMBER OF SEQUENCES IN EACH RING
ln(1,1)=1;
ln(1,2)=l(1);
for i=2:3
    ln(i,1)=ln(i-1,2)+1;
    ln(i,2)=ln(i,1)+l(i)-1;
end

[Ucode1,Int_indx]=sortrows(Ucode1);% Ucode1 is PASSED TO REMEMBER MISMATCH
Us1=Us1(Int_indx,:);
Int=dat(:,1);


