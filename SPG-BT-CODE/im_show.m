function im_show(I)

figure;

iptsetpref('ImshowBorder','tight');
imshow(I,[0 1]);

[m n]=size(I);

set(gcf, 'WindowStyle', 'normal');
set(gca, 'Unit', 'inches'); 
set(gca, 'Position', [0 0 4.5*(n/m) 4.5]); % image position and size
set(gcf, 'Unit', 'inches'); 
set(gcf, 'Position', [2 2 4.5*(n/m) 4.5]); % figure position and size
