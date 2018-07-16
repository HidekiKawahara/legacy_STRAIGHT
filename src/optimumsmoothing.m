function ovc=optimumsmoothing(eta,pc)
%	ovc=optimumsmoothing(eta,pc)
%	Calculate the optimum smoothing function
%	ovc	: coefficients for 2nd order cardinal B-spline
%	eta	: temporal stretch factor
%	pc	: power exponent for nonlinearity

%   05/July/2006 dirty patch. This routine has to be re-programmed.

fx=-8:0.05:8;
cb=max(0,1-abs(fx));
gw=exp(-pi*(fx*eta*1.4).^2).^pc;
cmw=conv(cb,gw);
bb=(1:length(cb));
bbc=bb+(length(cb)-1)/2;
cmw=cmw(bbc)/max(cmw);
ss=(abs(fx-round(fx))<0.025).*(1:length(cb));
ss=ss(ss>0);
cmws=cmw(ss);

nn=length(cmws);
idv=1:nn;

hh=zeros(2*nn,nn);
for ii=1:nn
  hh((ii-1)+idv,ii)=cmws';
end;
bv=zeros(2*nn,1);
bv(nn+1)=1;             % This is the original unit impulse.
h=hh'*hh;
ov = h \ (hh'*bv);

idc=(nn-1)/2+2;

ovc=ov(idc+(0:3));
