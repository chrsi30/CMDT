function incr(x)
% x++
evalin("caller",sprintf('%s = %s + 1;',x,x)  );
end