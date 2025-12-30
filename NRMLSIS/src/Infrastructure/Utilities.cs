// #######################################################################
// MSIS® (NRL-SOF-014-1) SOFTWARE
// NRLMSIS® empirical atmospheric model software. Use is governed by the
// Open Source Academic Research License Agreement contained in the file
// nrlmsis2.1_license.txt, which is part of this software package. BY
// USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
// CONDITIONS OF THE LICENSE.
// #######################################################################

// ===========================================================================
// NRLMSIS 2.1:
// Neutral atmosphere empirical model from the surface to lower exosphere
// ===========================================================================

// **************************************************************************************************
// MSIS_Utils: Contains the following auxiliary subroutines:
//  Alt2Gph:  Converts geodetic altitude to geopotential height
//  Gph2Alt:  Converts geopotential height to geodetic altitude
//  BSpline:  Computes B-splines using input nodes and up to specified order
//  Dilog:    Computes dilogarithm function (expansion truncated at order 3, error < 1E-5)
// **************************************************************************************************

using System;

namespace NRLMSIS.Infrastructure
{
    /// <summary>
    /// Utility functions for NRLMSIS 2.1
    /// </summary>
    public static class Utilities
    {
        // ==================================================================================================
        // ALT2GPH: Altitude to Geopotential Height
        // References:
        //   DMA Technical Report TR8350.2 (1987),
        //     http://earth-info.nga.mil/GandG/publications/historic/historic.html
        //   Featherstone, W. E., and S. J. Claessens (2008), Closed-form transformation between
        //     geodetic and ellipsoidal coordinates, Studia Geophysica et Geodaetica, 52, 1-18
        //   Jekeli, C. (2009), Potential theory and static gravity field of the Earth, in
        //     Treatise on Geophysics, ed. T. Herring, vol 3, 11-42
        //   NIMA Technical Report TR8350.2 (2000, 3rd edition, Amendment1),
        //     http://earth-info.nga.mil/GandG/publications/tr8350.2/tr8350_2.html
        // ==================================================================================================
        /// <summary>
        /// Converts geodetic altitude to geopotential height
        /// </summary>
        /// <param name="lat">Geodetic latitude (degrees)</param>
        /// <param name="alt">Geodetic altitude (km)</param>
        /// <returns>Geopotential height (km)</returns>
        public static double Alt2Gph(double lat, double alt)
        {
            const double Deg2Rad = 0.017453292519943295;

            // WGS84 Defining parameters
            const double A = 6378.1370e3;           // Semi-major axis of reference ellipsoid (m)
            const double FInv = 298.257223563;      // 1/f = Reciprocal of flattening
            const double W = 7292115e-11;           // Angular velocity of Earth rotation (rad/s)
            const double GM = 398600.4418e9;        // Gravitational constant x Earth mass (m^3/s^2)

            // WGS84 Derived parameters
            double aSq = A * A;
            double wSq = W * W;
            double f = 1.0 / FInv;
            double eSq = 2 * f - f * f;
            double e = Math.Sqrt(eSq);              // Ellipsoid eccentricity
            double elin = A * e;                    // Linear eccentricity of ellipsoid
            double elinSq = elin * elin;
            double epr = e / (1 - f);               // Second eccentricity
            double q0 = ((1.0 + 3.0 / (epr * epr)) * Math.Atan(epr) - 3.0 / epr) / 2.0; // DMA Technical Report tr8350.2, Eq. 3-25
            double u0 = -GM * Math.Atan(epr) / elin - wSq * aSq / 3.0; // Theoretical potential of reference ellipsoid (m^2/s^2)
            const double G0 = 9.80665;              // Standard gravity (m/s^2), CGPM 1901; WMO
            double gmDivElin = GM / elin;

            // Parameters for centrifugal potential taper
            double x0Sq = Math.Pow(2e7, 2);         // Axial distance squared at which tapering begins (m^2)
            double hSq = Math.Pow(1.2e7, 2);        // Relaxation scale length of taper (m^2)

            // Working variables
            double altM, sinSqLat, v, xSq, zSq;
            double rSqMinElinSq, uSq, cosSqDelta, eprU, atanEprU, q, U, Vc;

            // Compute Cartesian and ellipsoidal coordinates
            altM = alt * 1000.0;
            sinSqLat = Math.Pow(Math.Sin(lat * Deg2Rad), 2);
            v = A / Math.Sqrt(1 - eSq * sinSqLat);           // Radius of curvature of the reference ellipsoid, Featherstone eq. 4
            xSq = Math.Pow(v + altM, 2) * (1 - sinSqLat);    // Squared x-coordinate of geocentric system, Featherstone eq. 1
            zSq = Math.Pow(v * (1 - eSq) + altM, 2) * sinSqLat; // Squared z-coordinate of geocentric system, Featherstone eq. 3
            rSqMinElinSq = xSq + zSq - elinSq;
            uSq = rSqMinElinSq / 2.0 + Math.Sqrt(Math.Pow(rSqMinElinSq, 2) / 4.0 + elinSq * zSq); // Ellipsoidal distance coordinate, Featherstone eq. 19
            cosSqDelta = zSq / uSq;                          // Ellipsoidal polar angle, Featherstone eq. 21

            // Compute gravitational potential
            eprU = elin / Math.Sqrt(uSq);                    // Second eccentricity at ellipsoidal coordinate u
            atanEprU = Math.Atan(eprU);
            q = ((1 + 3.0 / (eprU * eprU)) * atanEprU - 3.0 / eprU) / 2.0; // Jekeli, eq. 114
            U = -gmDivElin * atanEprU - wSq * (aSq * q * (cosSqDelta - 1.0 / 3.0) / q0) / 2.0; // Jekeli, eq. 113

            // Compute centrifugal potential and adjust total potential
            if (xSq <= x0Sq)
            {
                Vc = (wSq / 2.0) * xSq;
            }
            else
            {
                Vc = (wSq / 2.0) * (hSq * Math.Tanh((xSq - x0Sq) / hSq) + x0Sq); // Centrifugal potential taper
            }
            U = U - Vc;

            // Compute geopotential height
            return (U - u0) / G0 / 1000.0;
        }

        // ==================================================================================================
        // GPH2ALT: Geopotential Height to Altitude
        // ==================================================================================================
        /// <summary>
        /// Converts geopotential height to geodetic altitude using Newton-Raphson iteration
        /// </summary>
        /// <param name="theta">Geodetic latitude (degrees)</param>
        /// <param name="gph">Geopotential height (km)</param>
        /// <returns>Geodetic altitude (km)</returns>
        public static double Gph2Alt(double theta, double gph)
        {
            const int MaxN = 10;
            const double Epsilon = 0.0005;

            double x = gph;
            int n = 0;
            double dx = Epsilon + Epsilon;

            while (Math.Abs(dx) > Epsilon && n < MaxN)
            {
                double y = Alt2Gph(theta, x);
                double dydz = (Alt2Gph(theta, x + dx) - y) / dx;
                dx = (gph - y) / dydz;
                x = x + dx;
                n = n + 1;
            }

            return x;
        }

        // ==================================================================================================
        // BSPLINE: Returns array of nonzero b-spline values, for all orders up to specified order (max 6)
        // ==================================================================================================
        /// <summary>
        /// Computes B-splines using input nodes and up to specified order
        /// </summary>
        /// <param name="x">Location at which splines are to be evaluated</param>
        /// <param name="nodes">Spline node locations (0-based array)</param>
        /// <param name="nd">Number of spline nodes minus one (nodes is 0:nd)</param>
        /// <param name="kmax">Maximum order (up to 6 allowed) of evaluated splines</param>
        /// <param name="eta">Precomputed weights for recursion (reciprocals of node differences) [0:30, 2:6]</param>
        /// <param name="S">Output: b-spline values [spline index relative to i (-5:0), spline order (2:6)]</param>
        /// <param name="i">Output: Index of last nonzero b-spline</param>
        public static void BSpline(double x, double[] nodes, int nd, int kmax,
                                   double[,] eta, out double[,] S, out int i)
        {
            // Initialize output array S[-5:0, 2:6]
            // In C#, we'll use [0:5, 0:4] and map indices: S[j+5, k-2]
            S = new double[6, 5]; // [-5:0] maps to [0:5], [2:6] maps to [0:4]

            // Find index of last (rightmost) nonzero spline
            if (x >= nodes[nd])
            {
                i = nd;
                return;
            }
            if (x <= nodes[0])
            {
                i = -1;
                return;
            }

            int low = 0;
            int high = nd;
            i = (low + high) / 2;
            while (x < nodes[i] || x >= nodes[i + 1])
            {
                if (x < nodes[i])
                {
                    high = i;
                }
                else
                {
                    low = i;
                }
                i = (low + high) / 2;
            }

            // Initialize with linear splines (k=2)
            S[5, 0] = (x - nodes[i]) * eta[i, 0]; // S[0,2] -> S[5,0]
            if (i > 0) S[4, 0] = 1 - S[5, 0];      // S[-1,2] -> S[4,0]
            if (i >= nd - 1) S[5, 0] = 0.0;        // Reset out-of-bounds spline to zero

            // k = 3 (quadratic splines)
            double[] w = new double[5]; // w[-4:0] maps to w[0:4]
            w[4] = (x - nodes[i]) * eta[i, 1];     // w[0] -> w[4]
            if (i != 0) w[3] = (x - nodes[i - 1]) * eta[i - 1, 1]; // w[-1] -> w[3]

            if (i < (nd - 2)) S[5, 1] = w[4] * S[5, 0]; // S[0,3]
            if (((i - 1) >= 0) && ((i - 1) < (nd - 2)))
                S[4, 1] = w[3] * S[4, 0] + (1.0 - w[4]) * S[5, 0]; // S[-1,3]
            if ((i - 2) >= 0) S[3, 1] = (1.0 - w[3]) * S[4, 0]; // S[-2,3]

            // k = 4 (cubic splines)
            for (int l = 0; l >= -2; l--)
            {
                int j = i + l;
                if (j < 0) break; // Skip out-of-bounds splines
                w[l + 4] = (x - nodes[j]) * eta[j, 2]; // w[l] -> w[l+4]
            }

            if (i < (nd - 3)) S[5, 2] = w[4] * S[5, 1]; // S[0,4]
            for (int l = -1; l >= -2; l--)
            {
                if (((i + l) >= 0) && ((i + l) < (nd - 3)))
                    S[l + 5, 2] = w[l + 4] * S[l + 5, 1] + (1.0 - w[l + 1 + 4]) * S[l + 1 + 5, 1];
            }
            if ((i - 3) >= 0) S[2, 2] = (1.0 - w[2]) * S[3, 1]; // S[-3,4]

            // k = 5
            for (int l = 0; l >= -3; l--)
            {
                int j = i + l;
                if (j < 0) break; // Skip out-of-bounds splines
                w[l + 4] = (x - nodes[j]) * eta[j, 3]; // w[l] -> w[l+4]
            }

            if (i < (nd - 4)) S[5, 3] = w[4] * S[5, 2]; // S[0,5]
            for (int l = -1; l >= -3; l--)
            {
                if (((i + l) >= 0) && ((i + l) < (nd - 4)))
                    S[l + 5, 3] = w[l + 4] * S[l + 5, 2] + (1.0 - w[l + 1 + 4]) * S[l + 1 + 5, 2];
            }
            if ((i - 4) >= 0) S[1, 3] = (1.0 - w[1]) * S[2, 2]; // S[-4,5]
            if (kmax == 5) return; // Exit if only 5th order spline is needed

            // k = 6
            for (int l = 0; l >= -4; l--)
            {
                int j = i + l;
                if (j < 0) break; // Skip out-of-bounds splines
                w[l + 4] = (x - nodes[j]) * eta[j, 4]; // w[l] -> w[l+4]
            }

            if (i < (nd - 5)) S[5, 4] = w[4] * S[5, 3]; // S[0,6]
            for (int l = -1; l >= -4; l--)
            {
                if (((i + l) >= 0) && ((i + l) < (nd - 5)))
                    S[l + 5, 4] = w[l + 4] * S[l + 5, 3] + (1.0 - w[l + 1 + 4]) * S[l + 1 + 5, 3];
            }
            if ((i - 5) >= 0) S[0, 4] = (1.0 - w[0]) * S[1, 3]; // S[-5,6]
        }

        // ==================================================================================================
        // DILOG: Calculate dilogarithm in the domain [0,1)
        // Retains terms up to order 3 in the expansion, which results in relative errors less than 1E-5.
        // Reference:
        //   Ginsberg, E. S., and D. Zaborowski (1975), The Dilogarithm function of a real argument,
        //   Commun. ACM, 18, 200–202.
        // ==================================================================================================
        /// <summary>
        /// Computes the dilogarithm function for arguments in [0,1)
        /// </summary>
        /// <param name="x0">Input argument in domain [0,1)</param>
        /// <returns>Dilogarithm value</returns>
        public static double Dilog(double x0)
        {
            double pi2_6 = Constants.Pi * Constants.Pi / 6.0;

            double x = x0;
            double dilog;

            if (x > 0.5)
            {
                double lnx = Math.Log(x);
                x = 1.0 - x; // Reflect argument into [0,0.5] range
                double xx = x * x;
                double x4 = 4.0 * x;
                dilog = pi2_6 - lnx * Math.Log(x)
                        - (4.0 * xx * (23.0 / 16.0 + x / 36.0 + xx / 576.0 + xx * x / 3600.0)
                            + x4 + 3.0 * (1.0 - xx) * lnx) / (1.0 + x4 + xx);
            }
            else
            {
                double xx = x * x;
                double x4 = 4.0 * x;
                dilog = (4.0 * xx * (23.0 / 16.0 + x / 36.0 + xx / 576.0 + xx * x / 3600.0)
                            + x4 + 3.0 * (1.0 - xx) * Math.Log(1.0 - x)) / (1.0 + x4 + xx);
            }

            return dilog;
        }

        // ==================================================================================================
        // Helper method to access S array with Fortran-style indexing
        // ==================================================================================================
        /// <summary>
        /// Helper to access S array with Fortran-style negative indexing
        /// </summary>
        /// <param name="S">The spline array</param>
        /// <param name="splineIdx">Spline index (-5 to 0)</param>
        /// <param name="order">Spline order (2 to 6)</param>
        /// <returns>The spline value</returns>
        public static double GetSValue(double[,] S, int splineIdx, int order)
        {
            return S[splineIdx + 5, order - 2];
        }

        /// <summary>
        /// Helper to set S array with Fortran-style negative indexing
        /// </summary>
        public static void SetSValue(double[,] S, int splineIdx, int order, double value)
        {
            S[splineIdx + 5, order - 2] = value;
        }
    }
}