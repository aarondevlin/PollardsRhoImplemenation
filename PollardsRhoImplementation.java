import java.math.BigInteger;
import java.util.Random;
import java.util.Scanner;

/*
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
*/

public class PollardsRho {

	public static void main(String[] args) {


		// Testing exp()

		BigInteger au[] = new BigInteger[2];
		au[0] = BigInteger.valueOf(12);
		au[1] = BigInteger.valueOf(61833);

		BigInteger mu = BigInteger.valueOf(2);
		BigInteger du = BigInteger.valueOf(154);
		BigInteger pu = BigInteger.valueOf(65519);

		BigInteger testexp[] = new BigInteger[2];
		testexp = exp(au, mu, du, pu);



		// USER INPUT PARAMETERS FOR a, d, n, & p

		Scanner reader = new Scanner(System.in);  // Reading from System.in

		System.out.println("Enter d (for example, 154):");
		int dd = reader.nextInt(); // Scans the next token of the input as an int.
		BigInteger d = BigInteger.valueOf(dd);

		System.out.println("Enter n (for example, 16339):");
		int nn = reader.nextInt(); // Scans the next token of the input as an int.
		BigInteger n = BigInteger.valueOf(nn);

		System.out.println("Enter p (for example, 65519):");
		int pp = reader.nextInt(); // Scans the next token of the input as an int.
		BigInteger p = BigInteger.valueOf(pp);

		System.out.println("Enter ax (for example, 12):");
		int aex = reader.nextInt(); // Scans the next token of the input as an int.
		BigInteger a[] = new BigInteger[2];
		a[0] = BigInteger.valueOf(aex);

		System.out.println("Enter ay (for example, 61833):");
		int awhy = reader.nextInt();
		a[1] = BigInteger.valueOf(awhy);

		reader.close();


		// Find the average of k steps needed to find m' = m in check() for N random discrete logarithms
		BigInteger testcheck;
		BigInteger sumk = BigInteger.ZERO;
		BigInteger avgk;

		int N = 100;
		BigInteger bign = BigInteger.valueOf(N);

		for (int i = 1;  i<=N; i++) {
			testcheck = check(a, d, p, n);
			sumk = sumk.add(testcheck);
		}

		avgk = sumk.divide(bign);

		System.out.println("Average # of K steps to find m' = m for N random discrete logarithms: " + avgk);
	}



	// mul() finds a3, which is the product of a1 and a2
	private static BigInteger[] mul(BigInteger[] a1var, BigInteger[] a2var, BigInteger dvar, BigInteger pvar) {

		BigInteger[] a3 = new BigInteger[2];
		BigInteger x1 = a1var[0];
		BigInteger y1 = a1var[1];
		BigInteger x2 = a2var[0];
		BigInteger y2 = a2var[1];
		BigInteger one = BigInteger.ONE;

		//calculate a3x
		BigInteger x1y2 = x1.multiply(y2).mod(pvar);
		BigInteger y1x2 = y1.multiply(x2).mod(pvar);
		BigInteger x3numr8r = x1y2.add(y1x2).mod(pvar);
		BigInteger dx1x2y1y2 = dvar.multiply(x1).mod(pvar).multiply(x2).mod(pvar).multiply(y1).mod(pvar).multiply(y2).mod(pvar);
		BigInteger x3dnmn8r = dx1x2y1y2.add(one).mod(pvar);
		BigInteger x3 = x3numr8r.multiply(x3dnmn8r.modInverse(pvar)).mod(pvar);
		a3[0] = x3;

		//calculate a3y
		BigInteger y1y2 = y1.multiply(y2).mod(pvar);
		BigInteger x1x2 = x1.multiply(x2).mod(pvar);
		BigInteger y3numr8r = y1y2.subtract(x1x2).mod(pvar);
		BigInteger y3dnmn8r = one.subtract(dx1x2y1y2).mod(pvar);
		BigInteger y3 = y3numr8r.multiply(y3dnmn8r.modInverse(pvar)).mod(pvar);  //BigInteger y3 = y3numr8r.divide(y3dnmn8r).mod(pvar); is wrong, cannot divide.
		a3[1] = y3;


		return a3;
	}



	// Exponentiation Algorithm
	// Requires O(k) group operations, where k is the number of bits in m
	// This is much faster than a function that multiplies a*a, looping until the desired exponent "m" is reached
	private static BigInteger[] exp(BigInteger[] avar, BigInteger mvar, BigInteger dvar, BigInteger pvar) {

		BigInteger bee[] = new BigInteger[2];
		int numofbits = mvar.bitLength();   // Find k
		bee[0] = BigInteger.ZERO;    // Initialize b = ( 0, 1 )
		bee[1] = BigInteger.ONE;

		for (int i = numofbits - 1; i >= 0; i--) {
			bee = mul(bee, bee, dvar, pvar);      // For every bit i in m, b = b * b

			if (mvar.testBit(i)) {      // If the i-th bit of m = 1, then b = b * a
				bee = mul(bee, avar, dvar, pvar);
			}

		}
		return bee;

	}



	// Function to create alphak+1, betak+1, zk+1 from alphak, betak, and zk based on the xcoord(zk)mod3 switch statements
	// Returns a BigInt array of size 4 containing alphak+1, betak+1, zk+1's x and y values
	private static BigInteger[] newabz(BigInteger alphak, BigInteger betak, BigInteger zkx, BigInteger zky, BigInteger[] avar, BigInteger[] bvar, BigInteger dvar, BigInteger pvar) {

		BigInteger abz[] = new BigInteger[4];  //stores the output alphak+1, betak+1, zk+1's x and y values
		BigInteger three = BigInteger.valueOf(3); //used for mod 3
		BigInteger two = BigInteger.valueOf(2); //used for * 2
		BigInteger zk[] = new BigInteger[2];  //stores the x and y values of zk
		zk[0] = zkx;
		zk[1] = zky;

		//Take the x coord of zk mod(3), look at the 3 cases in if statement to update alphak, betak and zk
		BigInteger zkxcoordm3 = zkx.mod(three);


		// Computing alphak, betak, and zk

		// xcoord(z)mod(3) if statements
		// If x coord of zk mod(3) = 0.  compareTo returns -1 if zkxcoordm3 < 0, 0 if zkxcoordm3 = 0, 1 if zkxcoordm3 > 0. So we want compareTo to return 0
		if (zkxcoordm3.compareTo(BigInteger.ZERO) == 0) {
			zk = mul(bvar, zk, dvar, pvar);
			alphak = alphak.add(BigInteger.ONE);
		}
		else if (zkxcoordm3.compareTo(BigInteger.ONE) == 0) {
			zk = mul(zk, zk, dvar, pvar);
			alphak = alphak.multiply(two);
			betak = betak.multiply(two);
		}
		else {
			zk = mul(avar, zk, dvar, pvar);
			betak = betak.add(BigInteger.ONE);
		}

		abz[0] = alphak;
		abz[1] = betak;
		abz[2] = zk[0];
		abz[3] = zk[1];

		return abz;
	}



	// RHO(): given a and b, find (m), where b = a ^ m, and (k) the number of steps of pollard's rho that it took to find m.
	private static BigInteger[] rho(BigInteger[] avar, BigInteger[] bvar, BigInteger dvar, BigInteger pvar, BigInteger nvar) {
		BigInteger[] mk = new BigInteger[2];  //BigInt array for the output m and k
		BigInteger m;
		BigInteger mnmr8r;
		BigInteger mdnm8r;

		// Initialize alpha0, beta0, and z0, where z0 is a BigInt array of size 2 which contains x and y values
		BigInteger alphakay = BigInteger.ZERO;  //initialize alphak as alpha0 = 0
		BigInteger betakay = BigInteger.ZERO;  //initialize betak as beta0 = 0
		BigInteger zkay[] = new BigInteger[2];
		zkay[0]= BigInteger.ZERO;  //initialize zk as z0 = (0, 1)
		zkay[1]= BigInteger.ONE;


		// Initialize variables for k and 2k values
		BigInteger kvals[] = new BigInteger[4];
		BigInteger ktimes2vals[] = new BigInteger[4];

		// Initialize alphak, betak, zk at k = 0
		kvals[0] = alphakay;
		kvals[1] = betakay;
		kvals[2] = zkay[0];
		kvals[3] = zkay[1];


		// Initialize alpha2k, beta2k, z2k at k = 0
		ktimes2vals = kvals;


		//Find the Discrete Logarithm

		int n = nvar.intValue(); //convert nvar to an integer for comparison with k in the coming for loop
		int steps = 0;

		for (int k = 1; k<=n; k++) {
			steps++;
			kvals = newabz(kvals[0], kvals[1], kvals[2], kvals[3], avar, bvar, dvar, pvar);  //find alphak, betak, zk at k
			ktimes2vals = newabz(ktimes2vals[0], ktimes2vals[1], ktimes2vals[2], ktimes2vals[3], avar, bvar, dvar, pvar); //find //find alpha2k, beta2k, z2k at k
			ktimes2vals = newabz(ktimes2vals[0], ktimes2vals[1], ktimes2vals[2], ktimes2vals[3], avar, bvar, dvar, pvar);

			// Check to see if zk = z2k
			if (kvals[2].compareTo(ktimes2vals[2]) == 0) {
				if (kvals[3].compareTo(ktimes2vals[3]) == 0) {
					mnmr8r = ktimes2vals[1].subtract(kvals[1]).mod(nvar);
					mdnm8r = kvals[0].subtract(ktimes2vals[0]).mod(nvar);
					m = mnmr8r.multiply(mdnm8r.modInverse(nvar)).mod(nvar); // used to be mnmr8r.multiply(mdnm8r.modInverse(pvar)).mod(pvar).mod(nvar);

					mk[0] = m;
					// K must be of BigInt type...
					BigInteger bigk = BigInteger.valueOf(steps);
					mk[1] = bigk;

					break;
				}
			}
		}

		return mk;
	}



	// Check() creates a random m, then from this it creates b. Then m' is computed with rho(), then validates m' is correct.
	// The number of steps required, "k", is returned.
	private static BigInteger check(BigInteger[] avar, BigInteger dvar, BigInteger pvar, BigInteger nvar) {

			Random r = new Random();
			int intr = r.nextInt(nvar.intValue())+1; //make a new random number 1 <= intr <= n, since nextInt makes a random # from 0 <= # < arg, so we need the +1 so that we dont get 0
			BigInteger randm = BigInteger.valueOf(intr);

			System.out.println("------------------------random m:");
			System.out.println(randm);

			// Compute b with the new m and exp(a, m, d, p)
			BigInteger[] b = new BigInteger[2];
			b = exp(avar, randm, dvar, pvar);

			// Compute the discrete logarithm "m'" with rho(a, b, d, p, n)
			BigInteger[] mprime = new BigInteger[2]; //store the output of rho
			mprime = rho(avar, b, dvar, pvar, nvar);

			System.out.println("------------------------from rho(), m':");
			System.out.println(mprime[0]);
			System.out.println("k steps:");
			System.out.println(mprime[1]);

			// Check to see if the discrete logarithm m' = m, where m is the actual m
			if (mprime[0].compareTo(randm) != 0) {
				throw new RuntimeException();
			}

			return mprime[1];
	}


}
