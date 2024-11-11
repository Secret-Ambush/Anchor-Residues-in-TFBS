% AUTHOR OF THE PROGRAM - DEVESH BHIMSARIA

% PROGRAM TO GENERATE SEQUENCE SPECIFICITY AND ENERGY LANDSCAPES FOR 12MERS


% inp         = INPUT FILE - THIS IS THE OUTPUT FILE FROM ssl_4.pl OR ssl_5.pl 
% op          = OUTPUT FILE
% motif       = MOTIF SHOULD CONTAIN A,C,G OR T ONLY
% l_l         = THICKNESS OF EACH SSL RING
% sm          = SMOOTHENING 
% ring_start  =1 (0 mismatch ring) , 2 (1 mismatch ring), 0 (default i.e. 1)%
% ring_stop   =1 (0 mismatch ring) , 2 (1 mismatch ring), 0 (default i.e. till last ring)
% colum       =1 (1st column after forward and reverse sequence), 2(2nd)
% ls          = 10 or 12 (length of k-mer input file)

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

% EXAMPLE - ssl_v8_perl_v2_s1('ATF3.txt', 'abc.txt', 'TGACGTCA',1,0,1,1,12)
function ssl_v8_perl_v2_s1(inp,op,motif,l_l,sm,sel,colum,ls)

motif=upper(motif);
lm=length(motif);
ring_start=1;
ring_stop=3;

[Int,Int_indx,~,Us3,l,ln,n,ls2]=...
    get_mismatches_v6_perl_v2_s1(inp,motif,lm,ls,sel,colum);%ARRANGE THE SEQUENCES AND GET MISMATCH

Intx=Int(Int_indx);% SORTING INTENSITY IN ORDER OF SEQUENCES
Int_sort=sort(Int,'descend'); mx=Int_sort(10);

[~,~,~,f]=draw_ssl_v6_s1(Intx,Us3,lm,sm,l_l,op,n,l,ln,ls2,ring_start,ring_stop,mx);
set(f,'OuterPosition',[500 200 400 500])
    
end
