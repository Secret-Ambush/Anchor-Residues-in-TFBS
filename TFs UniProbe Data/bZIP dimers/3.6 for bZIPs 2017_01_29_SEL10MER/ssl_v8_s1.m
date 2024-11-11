% AUTHOR OF THE PROGRAM - DEVESH BHIMSARIA


% PROGRAM TO GENERATE SEQUENCE SPECIFICITY AND ENERGY LANDSCAPES USING
% 10mer FILES AS INPUT

% inp         = INPUT FILE
% op          = OUTPUT FILE
% motif       = MOTIF
% l_l         = THICKNESS OF EACH SSL RING
% sm          = SMOOTHENING 
% ring_start  =1 (0 mismatch ring) , 2 (1 mismatch ring), 0 (default i.e. 1)%
% ring_stop   =1 (0 mismatch ring) , 2 (1 mismatch ring), 0 (default i.e. till last ring)
% colum       =1 (1st column after forward and reverse sequence), 2(2nd)

% sel         = SELECTING ORDERING OF SEQUENCES IN DIFFERENT RIGNS RING
%               DIFF INPUTS AND CORRESPONDING PREFERENCES IN ORDERING 
% 	1- 1st- flanking on 3? end of motif. 2nd- flanking on 5? end of motif. 3rd- mismatched base. 4th- position of motif in sequence.(default)
% 	2- 1st- flanking on 5? end of motif. 2nd- flanking on 3? end of motif. 3rd- mismatched base. 4th- position of motif in sequence.
% 	3- 1st- mismatched base. 2nd- flanking on 3? end of motif. 3rd- flanking on 5? end of motif. 4th- position of motif in sequence.
% 	4- 1st- mismatched base. 1st- flanking on 5? end of motif. 2nd- flanking on 3? end of motif. 4th- position of motif in sequence.
% 	5- 1st- position of motif in sequence. 2nd- flanking on 3? end of motif. 3rd- flanking on 5? end of motif. 4th- mismatched base. 
% 	6- 1st- position of motif in sequence. 2nd- flanking on 5? end of motif. 3rd- flanking on 3? end of motif. 4th- mismatched base. 
% 	7- 1st- position of motif in sequence. 2nd- mismatched base. 3rd- flanking on 3? end of motif. 4th- flanking on 5? end of motif.
% 	8- 1st- position of motif in sequence. 2nd- mismatched base. 3rd- flanking on 5? end of motif. 4th- flanking on 3? end of motif.

% EXAMPLE - ssl_v8_s1('BATF2_JUN.txt', 'abc.txt', 'TGACGTCA',1,0,1,1,4,1)
function ssl_v8_s1(inp,op,motif,l_l,sm,sel,ring_start,ring_stop,colum)

motif=upper(motif);
motif=motif_degenerate2_2(motif); % TO CONVERT DEGENERATE NUCLEOTIDES TO ACGT

[n,mism,pos_mis_mot,match_to_mot,mot_fr,F,R,Int,ls2,lm]=read_datasets_v2_2(inp,motif,colum);%READ FILE AND MOTIF INFORMATION
[Int_indx,~,Us3,l,ln]=get_mismatches_v8(n,mism,pos_mis_mot,match_to_mot,mot_fr,F,R,ls2,lm,motif,sel);%ARRANGE THE SEQUENCES AND GET MISMATCH

Intx=Int(Int_indx);% SORTING INTENSITY IN ORDER OF SEQUENCES
Int_sort=sort(Int,'descend'); mx=Int_sort(10);

[~,~,~,f]=draw_ssl_v6_s1(Intx,Us3,lm,sm,l_l,op,n,l,ln,ls2,ring_start,ring_stop,mx);
set(f,'OuterPosition',[500 200 400 500])
    
end
