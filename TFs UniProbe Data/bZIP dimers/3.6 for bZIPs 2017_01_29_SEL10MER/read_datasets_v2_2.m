% FUNCTION TO READ DATASET AND INPUT FILE

function  [n,mism,pos_match,match_to_mot,mot_fr,F,R,Int,ls2,lm]=read_datasets_v2_2(inp,motif,colum)

%ns=NUMBER OF MOTIFS lm=LENGTH OF MOTIF
[ns,lm]=size(motif);

%LOAD DATASET FOR LENGTH OF MOTIF
xx=['datasetxx_' int2str(lm)];
if exist(xx,'file')
    error('Motif length should be >=3 & <=10');
end
load(xx);

nn=n; % TOTAL NMERS (HERE NMER=lm)
clear n;

pos_mis_mot=zeros(2*ns,lm);% position where mismatch occurred for each motif
pos_match=zeros(nn,lm);% final positions where match occurred for each sequence of length lm (length of motif)

mism=zeros(nn,1);% KEEPS TRACK OF MINIMUM MISMATCH
match_to_mot=zeros(nn,1); %KEEPS TRACK OF MAXIMUM MATCH WITH WHICH MOTIF
mot_fr=zeros(nn,1); %MATCHED MOTIF IS FORWARD OR REVERSE

mm=zeros(2*ns,1);% TEMPORARY VARIABLE STORES THE NUMBER OF MISMATCHES FOR EACH MOTIF FOR EACH SEQUENCE

% FOR EACH NMER (HERE NMER=lm)
for i=1:nn
    
    % FINDING WHICH MOTIF CAUSED MINIMUM MISMATCH & STORING IT & MATCH POS.
    for nm=1:ns
        pos_mis_mot(2*nm-1,:)=(seq_frwd(i,:)-motif(nm,:))~=0;
        mm(2*nm-1)=length(find(pos_mis_mot(2*nm-1,:)));
        pos_mis_mot(2*nm,:)=(seq_rev(i,:)-motif(nm,:))~=0;
        mm(2*nm)=length(find(pos_mis_mot(2*nm,:)));
    end
    
    [mism(i), ind]=min(mm);
    match_to_mot(i)=ceil(ind/2);
    mot_fr(i)=rem(ind,2); % 1== frwd ; 0== rev
    
    ind_t=2*match_to_mot(i)-1; % temporary variable DOUBLE OF FROM WHICH MOTIF IT CAME
    
    % IF MISMATCH IS SAME FOR FRWD & REV GET THE ONE IN WHICH MISMATCH
    % OCCURRED FIRST (SPECIALLY FOR PALINDROME MOTIFS)
    if mm(ind_t)== mm(ind_t+1)
        if (find(pos_mis_mot(ind_t,:)>pos_mis_mot(ind_t+1,:),1)>find(pos_mis_mot(ind_t,:)<pos_mis_mot(ind_t+1,:),1))% IF REV MISMATCH CAME FIRST
            mot_fr(i)=0;
            pos_match(i,:)=pos_mis_mot(ind_t+1,:);
        else% IF FRWD MISMATCH CAME FIRST OR EQUAL
            mot_fr(i)=1;
            pos_match(i,:)=pos_mis_mot(ind_t,:);
        end
    else
        pos_match(i,:)=pos_mis_mot(ind,:);
    end
    
end

%READ CSI FILE
[Fx,Rx,Int,ls,n]=read_csi_file(inp,colum);

% ADD lm-1 Xs TO START AND END OF SEQUENCE
ls2=ls+2*(lm-1);% LENGTH OF SEQUENCE AFTER ADDING Xs
X_ascii='X'+0;
F=X_ascii*ones(n,ls2);
R=X_ascii*ones(n,ls2);
F(:,lm:lm+ls-1)=Fx;
R(:,lm:lm+ls-1)=Rx;

end