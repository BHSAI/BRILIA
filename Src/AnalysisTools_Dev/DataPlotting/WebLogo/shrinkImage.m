%function IMr = imshrink(IM,H,W) will shrink image IM into another image of
%size HxW pixels.

function IMr = shrinkImage(IM,H,W)
H1 = size(IM,1);
W1 = size(IM,2);
Hidx = round(linspace(1,H1,H));
Widx = round(linspace(1,W1,W));
IMtemp = IM(Hidx,:,:);
IMr = IMtemp(:,Widx,:);
