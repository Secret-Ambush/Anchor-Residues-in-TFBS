% FUNCTION TO CONVERT MOTIFS INTO ACGTs
function [motif_n]=motif_degenerate2_2(motif)
[a,lm]=size(motif);
motif_n=[];
for i=1:a
    m=motif(i,:);
    r=1;
    for j=1:lm
        if motif(i,j)=='A' || motif(i,j)=='C' || motif(i,j)=='G' || motif(i,j)=='T';
        elseif motif(i,j)=='W'
            m=[m;m];
            m(1:r,j)='A';
            m(r+1:2*r,j)='T';
            r=r*2;
        elseif motif(i,j)=='S'
            m=[m;m];
            m(1:r,j)='C';
            m(r+1:2*r,j)='G';
            r=r*2;
        elseif motif(i,j)=='Y'
            m=[m;m];
            m(1:r,j)='C';
            m(r+1:2*r,j)='T';
            r=r*2;
        elseif motif(i,j)=='R'
            m=[m;m];
            m(1:r,j)='A';
            m(r+1:2*r,j)='G';
            r=r*2;
        elseif motif(i,j)=='K'
            m=[m;m];
            m(1:r,j)='G';
            m(r+1:2*r,j)='T';
            r=r*2;
        elseif motif(i,j)=='M'
            m=[m;m];
            m(1:r,j)='A';
            m(r+1:2*r,j)='C';
            r=r*2;
        elseif motif(i,j)=='B'
            m=[m;m;m];
            m(1:r,j)='C';
            m(r+1:2*r,j)='G';
            m(2*r+1:3*r,j)='T';
            r=r*3;
        elseif motif(i,j)=='D'
            m=[m;m;m];
            m(1:r,j)='A';
            m(r+1:2*r,j)='G';
            m(2*r+1:3*r,j)='T';
            r=r*3;
        elseif motif(i,j)=='H'
            m=[m;m;m];
            m(1:r,j)='A';
            m(r+1:2*r,j)='C';
            m(2*r+1:3*r,j)='T';
            r=r*3;
        elseif motif(i,j)=='V'
            m=[m;m;m];
            m(1:r,j)='A';
            m(r+1:2*r,j)='C';
            m(2*r+1:3*r,j)='G';
            r=r*3;
        elseif motif(i,j)=='N'
            m=[m;m;m;m];
            m(1:r,j)='A';
            m(r+1:2*r,j)='C';
            m(2*r+1:3*r,j)='G';
            m(3*r+1:4*r,j)='T';
            r=r*4;
        else
            error('Incorrect input for motif: Motif should only contain ACGTWSYRKMBDHVN');
        end
    end
    motif_n=[motif_n; m];
end
motif_n=sortrows(motif_n);
