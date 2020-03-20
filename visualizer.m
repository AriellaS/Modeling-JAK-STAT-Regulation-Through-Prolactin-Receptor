figure(1)
startA = 1;
startB = 50;
numvars = 10;
for i = 0 : numvars - 1
	for j = 0 : i
		pos = i * numvars + j + 1;
		subplot(numvars, numvars, pos);
		a = results(:, length(predTime), i + startA);
		b = results(:, length(predTime), j + startB);
		scatter(b, a);
		if j == 0
			ylabel(i + startA);
		end
		if i == numvars - 1
			xlabel(j + startB);
		end
	end
end
