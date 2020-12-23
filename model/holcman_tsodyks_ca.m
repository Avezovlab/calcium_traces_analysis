tau = 0.5; %s
tr = 5; %s
U = 0.5;
sigma = 2.2; %mV
wt = 12.6; %mV/Hz
T = 2; %mV
alpha = 5; %Hz/mV


dt = 0.0001;
DT = 0.33;
trat = DT / dt;

N = 1750000;

Is = 1;
mus = 1;
cpt = 1;

cur_I = Is(end);
cur_mu = mus(end);
for n=2:N
    R = 0;
    if cur_I > T
        R = alpha * (cur_I - T);
    end

    cur_I = cur_I + (-cur_I(end) + cur_mu * U * wt * R) * dt/tau + sqrt(dt/tau) * sigma * randn;
    cur_mu = cur_mu + ((1 - cur_mu) / tr - U * cur_mu * R) * dt;
    
    if mod(cpt, trat) == 0
        Is = [Is; cur_I];
        mus = [mus; cur_mu];
        cpt = 1;
    else
        cpt = cpt + 1;
    end
end

figure
subplot(1,2,1)
plot((0:(length(Is)-1)) * DT, Is)
axis square
subplot(1,2,2)
plot((0:(length(Is)-1)) * DT, mus)
axis square

figure
plot(mus, Is)
axis square