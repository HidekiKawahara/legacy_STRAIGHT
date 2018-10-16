function y=getvalufromedit(co,defv)

ss=get(gco,'String');
y=str2num(ss);
if (length(y) <1) | (length(y)>1)
  y=defv;
end;
