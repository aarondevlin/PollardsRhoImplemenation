POLLARD'S RHO ALGORITHM IMPLEMENTATION
	Pollard's Rho considers (x, y) points on an elliptic curve ( ie. x^2 + y^2 = 1 + d*(x^2)(y^2) )
where d is given. Pollard's Rho algorithm recovers the discrete logarithm ("m" in b = a^m) when
given two cartesian points "a" and "b".  From "a" and "b", alpha, beta and z = (x,y) are created.
These values are updated based on the x value of z in the previous iteration (as shown in the
newabz() function. Values for k and for 2k are maintained at each iteration, so that a check can
be performed to see if zk = z2k. The 2k values of alpha, beta and z values are calculated by
using newabz() twice on the previous 2k values, where the first set of z2k values are initialized
by making them equal to the same initialization points as the k set: alpha = 0, beta = 0, z = (0, 1).
Once zk = z2k, Pollard's Rho Algorithm has found "m", and it is calculated with the k and 2k alpha
and beta values from that iteration.
	This program lets the user enter d, n, p, & a = (x, y).  The average number of steps, "k", that are
needed to find the discrete logarithm m is calculated Based off of the average amount of steps for
"N" runs.
