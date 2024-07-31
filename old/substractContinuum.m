function out = substractContinuum(Lq,u)
    [fitresult, ~] = createFit(Lq,u);
    for ii = 1:6
     C(ii) = getfield(fitresult,['c', num2str(ii)]);
     A(ii) = getfield(fitresult,['a', num2str(ii)]);
     B(ii) = getfield(fitresult,['b', num2str(ii)]);
    end
[C,idx] = sort(C,'descend'); A = A(idx); B = B(idx);
% Keep only first 3
f = @(x) A(1) * exp(-((x-B(1))/C(1)).^2) + A(2) * exp(-((x-B(2))/C(2)).^2) + A(3) * exp(-((x-B(3))/C(3)).^2);

out = u - f(Lq);
end