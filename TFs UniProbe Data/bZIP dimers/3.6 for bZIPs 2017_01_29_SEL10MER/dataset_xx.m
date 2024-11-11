% program to generate dataset mat files

nmers=[3 4 5 6 7 8 9 10]; % space gapped nmers... to create dataset.mat file
len_nmers=length(nmers);

seq_map='ACGTX';
rev_map='TGCAX';

for j=1:len_nmers
    
    nmer=nmers(j); %nmer FOR CURRENT LOOP
    disp(nmer);
    
    n=5^nmer; % TOTAL NUMBER OF NMER BASE=5
        
    rev_seq_map=zeros(n,1); 
    % == 0      IF PALINDROME 
    % ==-1      IF IT IS REVERSE COMPLEMENT OF SOMETHING ALREADY OCCURRED, AND
    % ==pos_no  TELLS THE INDEX OF IT'S REVERSE COMPLEMENT
    
    seq_frwd=zeros(n,nmer); %FORWARD SEQUENCE
    seq_rev=zeros(n,nmer);  %FORWARD SEQUENCE
    
    for i=1:n
        x=i-1;
        indn_rev=0;
        
        % FINDING FORWARD SEQUENCE seq_frwd & FORWARD SEQUENCE seq_rev
        for ii=1:nmer
            reme=rem(x,5);
            seq_frwd(i,ii)=seq_map(reme+1);
            x=(x-reme)/5;
            seq_rev(i,nmer-ii+1)=rev_map(reme+1);
            indn_rev=4*indn_rev+((4-1)-reme);
        end
        
        if indn_rev<i-1 % ALREADY OCCURRED
            rev_seq_map(i)=-1;
        elseif indn_rev==i-1 % PALINDROME
            rev_seq_map(i)=0;
        else
            rev_seq_map(i)=indn_rev+1; %INDEX OF REVERSE COMPLEMENT
        end
    end
    xx=['datasetxx_' int2str(nmer)];
    save(xx,'n','rev_seq_map','seq_frwd','seq_rev');
    
end
