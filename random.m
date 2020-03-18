% returns a random number between n/10 and n*10 on a uniform distribution
function random = random(n)
	random = n/10 + ((n*10 - n/10) * rand);
end

