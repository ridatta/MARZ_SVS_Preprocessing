function out = gaussianFn(x,mu,sigma,A)
 out = A .* exp( -1 * (x - mu).^2 / (2 * sigma.^2)); 
end