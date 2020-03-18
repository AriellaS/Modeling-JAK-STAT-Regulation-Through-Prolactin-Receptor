figure(1)
first = 20;
last = 30;
numvars = last - first + 1;
for i = first : last
	for j = first : i
		pos = (i-first) * (numvars) + (j-first) + 1;
		subplot(numvars,numvars,pos);
		a = results(:,length(predTime),i);
		b = results(:,length(predTime),j);
		scatter(b,a);
		if j==first
			ylabel(i);
		end
		if i==last
			xlabel(j);
		end
	end
end
