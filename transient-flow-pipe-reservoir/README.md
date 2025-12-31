This solves transient flow through pipe connected to a reservoir and is open to
atmosphere. 

For this, we write force equation in a pipe control volume that reveals
V_(i+1) = V_(i) + Delta t * ( H - (1 + f * L/D) * V_i^2/(2g) ) * g/L.

For friction factor, we use Swamee-Jain formula,
f = (0.25) / (log(k_s/(3.7 D) + 5.74/(Re^(0.9))))^2

We use a ratio that calculates % left to reach steady value. This value can be
used to figure out when to stop above calculations.
