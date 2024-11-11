% FUNCTION TO CHANGE COLORMAP FOR THE LANDSCAPES

function colormp2(mx)
colormap(jet(1024));
cx=caxis;
lmt=round((0-cx(1))/(cx(2)-0)*1024*(cx(2)/mx))+1; % lmt FOR THE NUMBER OF FEATURES AS GREY
lmt2=min(round((cx(2))/mx*1024),1024); % lmt2 FOR EXTENDING THE MAX COLOR LIMIT TO MAX INTENSITY
a=colormap;
ncmap=[kron([0.9529 0.9529 0.9529 ],ones(lmt,1));a(1:lmt2,:)];
colormap(ncmap);
