% returns an array of n random numbers between x/10 and x*10 on a log-uniform distribution
function random = random(x, n)
	random = x/10 * 10.^(2 * rand(n, 1));
end

