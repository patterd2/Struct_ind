function [x,deltas] = opi_dyn_fernando_uncoupled(maxt,Na,x0,gamma,eta,eta2,sigma,sigma_ND,W1,u,b,alpha,delta0)

deltas=cell([1,maxt+1]);
x=zeros(Na,maxt+1);
x(:,1) = x0;
dBt = randn(Na,maxt);
deltas(1,1) = {delta0};

for t=1:maxt
    deltas{1,t+1} = delta0;

    for j = 1:Na
        DaDo = 0;
        SaDo = alpha(j,1) * x(j,t);
        for k = 1:Na
            if k ~= j
                del = deltas{1,t};
                DaDo = DaDo + del(j,k)*x(k,t);
            end
        end
        y = SaDo + DaDo;
        S = tanh(W1 * y);
        x(j,t+1) = (gamma * x(j,t) + b(j,1) + u * S) ...
            + sigma*dBt(j,t);

    end

    %% ddelta/dt
    %          fx = tanh(W2*abs(x(:,t) - x(:,t)'));
    %keyboard;
    %          deltas(1,t+1) = {deltas{1,t} + dt * (-deltas{1,t}/tau + fx/tau)};
    %deltas(1,t+1) = {delta0};
end

end
