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
// MSIS_Gfn: Contains subroutines to calculate global (horizontal and time-dependent) model
//           basis functions
// **************************************************************************************************

using System;

namespace NRLMSIS
{
    /// <summary>
    /// Global basis functions for NRLMSIS 2.1
    /// </summary>
    public static class BasisFunctions
    {
        // Module-level variables for caching
        private static double[,] plg = new double[Constants.MaxN + 1, Constants.MaxN + 1]; // [0:maxn, 0:maxn]
        private static double[] cdoy = new double[2]; // [0:1] for indices 1-2 in Fortran
        private static double[] sdoy = new double[2]; // [0:1] for indices 1-2 in Fortran
        private static double[] clst = new double[3]; // [0:2] for indices 1-3 in Fortran
        private static double[] slst = new double[3]; // [0:2] for indices 1-3 in Fortran
        private static double[] clon = new double[2]; // [0:1] for indices 1-2 in Fortran
        private static double[] slon = new double[2]; // [0:1] for indices 1-2 in Fortran
        private static double sfluxAvgRef = 150.0; // Reference F10.7 value (=150 in NRLMSISE-00)
        private static double sfluxAvgQuadCutoff = 150.0; // Cutoff F10.7 for truncated quadratic F10.7a function
        private static double lastLat = -999.9;
        private static double lastDoy = -999.9;
        private static double lastLst = -999.9;
        private static double lastLon = -999.9;

        // ==================================================================================================
        // GLOBE: Calculate horizontal and time-dependent basis functions
        //        (Same purpose as NRLMSISE-00 "GLOBE7" subroutine)
        // ==================================================================================================
        /// <summary>
        /// Calculate horizontal and time-dependent basis functions
        /// </summary>
        /// <param name="doy">Day of year</param>
        /// <param name="utsec">Universal time in seconds</param>
        /// <param name="lat">Latitude (degrees)</param>
        /// <param name="lon">Longitude (degrees)</param>
        /// <param name="sfluxavg">81-day average F10.7</param>
        /// <param name="sflux">Daily F10.7</param>
        /// <param name="ap">Ap geomagnetic activity index history array [0:6] (7 elements)</param>
        /// <param name="bf">Output: basis function terms [0:maxnbf-1]</param>
        public static void Globe(double doy, double utsec, double lat, double lon,
                                double sfluxavg, double sflux, double[] ap, out double[] bf)
        {
            bf = new double[Constants.MaxNbf];

            double lst;
            double slat, clat, clat2, clat4, slat2;
            double cosdoy, sindoy;
            double coslon, sinlon;
            double pl;
            double coslst, sinlst;
            double dfa, df;
            double sza;
            int c;

            // Associated Legendre polynomials
            if (lat != lastLat)
            {
                clat = Math.Sin(lat * Constants.Deg2Rad);  // clat <=> sin, Legendre polynomial defined in colat
                slat = Math.Cos(lat * Constants.Deg2Rad);  // slat <=> cos, Legendre polynomial defined in colat
                clat2 = clat * clat;
                clat4 = clat2 * clat2;
                slat2 = slat * slat;

                plg[0, 0] = 1.0;
                plg[1, 0] = clat;
                plg[2, 0] = 0.5 * (3.0 * clat2 - 1.0);
                plg[3, 0] = 0.5 * (5.0 * clat * clat2 - 3.0 * clat);
                plg[4, 0] = (35.0 * clat4 - 30.0 * clat2 + 3.0) / 8.0;
                plg[5, 0] = (63.0 * clat2 * clat2 * clat - 70.0 * clat2 * clat + 15.0 * clat) / 8.0;
                plg[6, 0] = (11.0 * clat * plg[5, 0] - 5.0 * plg[4, 0]) / 6.0;

                plg[1, 1] = slat;
                plg[2, 1] = 3.0 * clat * slat;
                plg[3, 1] = 1.5 * (5.0 * clat2 - 1.0) * slat;
                plg[4, 1] = 2.5 * (7.0 * clat2 * clat - 3.0 * clat) * slat;
                plg[5, 1] = 1.875 * (21.0 * clat4 - 14.0 * clat2 + 1.0) * slat;
                plg[6, 1] = (11.0 * clat * plg[5, 1] - 6.0 * plg[4, 1]) / 5.0;

                plg[2, 2] = 3.0 * slat2;
                plg[3, 2] = 15.0 * slat2 * clat;
                plg[4, 2] = 7.5 * (7.0 * clat2 - 1.0) * slat2;
                plg[5, 2] = 3.0 * clat * plg[4, 2] - 2.0 * plg[3, 2];
                plg[6, 2] = (11.0 * clat * plg[5, 2] - 7.0 * plg[4, 2]) / 4.0;

                plg[3, 3] = 15.0 * slat2 * slat;
                plg[4, 3] = 105.0 * slat2 * slat * clat;
                plg[5, 3] = (9.0 * clat * plg[4, 3] - 7.0 * plg[3, 3]) / 2.0;
                plg[6, 3] = (11.0 * clat * plg[5, 3] - 8.0 * plg[4, 3]) / 3.0;

                lastLat = lat;
            }

            // Fourier harmonics of day of year
            if (doy != lastDoy)
            {
                cdoy[0] = Math.Cos(Constants.Doy2Rad * doy);
                sdoy[0] = Math.Sin(Constants.Doy2Rad * doy);
                cdoy[1] = Math.Cos(Constants.Doy2Rad * doy * 2.0);
                sdoy[1] = Math.Sin(Constants.Doy2Rad * doy * 2.0);
                lastDoy = doy;
            }

            // Fourier harmonics of local time
            lst = (utsec / 3600.0 + lon / 15.0) % 24.0;
            if (lst < 0) lst += 24.0;
            if (lst != lastLst)
            {
                clst[0] = Math.Cos(Constants.Lst2Rad * lst);
                slst[0] = Math.Sin(Constants.Lst2Rad * lst);
                clst[1] = Math.Cos(Constants.Lst2Rad * lst * 2.0);
                slst[1] = Math.Sin(Constants.Lst2Rad * lst * 2.0);
                clst[2] = Math.Cos(Constants.Lst2Rad * lst * 3.0);
                slst[2] = Math.Sin(Constants.Lst2Rad * lst * 3.0);
                lastLst = lst;
            }

            // Fourier harmonics of longitude
            if (lon != lastLon)
            {
                clon[0] = Math.Cos(Constants.Deg2Rad * lon);
                slon[0] = Math.Sin(Constants.Deg2Rad * lon);
                clon[1] = Math.Cos(Constants.Deg2Rad * lon * 2.0);
                slon[1] = Math.Sin(Constants.Deg2Rad * lon * 2.0);
                lastLon = lon;
            }

            //---------------------------------------------
            // Coupled Linear Terms
            //---------------------------------------------

            // Reset basis functions
            Array.Clear(bf, 0, bf.Length);

            // Time-independent (pure latitude dependence)
            c = Constants.CTimeInd;
            for (int n = 0; n <= Constants.AMaxN; n++)
            {
                bf[c] = plg[n, 0];
                c++;
            }

            // Intra-annual (annual and semiannual)
            if (c != Constants.CIntAnn) throw new InvalidOperationException("Problem with basis definitions");
            for (int s = 1; s <= Constants.AMaxS; s++)
            {
                cosdoy = cdoy[s - 1];
                sindoy = sdoy[s - 1];
                for (int n = 0; n <= Constants.AMaxN; n++)
                {
                    pl = plg[n, 0];
                    bf[c] = pl * cosdoy;
                    bf[c + 1] = pl * sindoy;
                    c += 2;
                }
            }

            // Migrating Tides (local time dependence)
            if (c != Constants.CTide) throw new InvalidOperationException("Problem with basis definitions");
            for (int l = 1; l <= Constants.TMaxL; l++)
            {
                coslst = clst[l - 1];
                sinlst = slst[l - 1];
                for (int n = l; n <= Constants.TMaxN; n++)
                {
                    pl = plg[n, l];
                    bf[c] = pl * coslst;
                    bf[c + 1] = pl * sinlst;
                    c += 2;
                }
                // Intra-annual modulation of tides
                for (int s = 1; s <= Constants.TMaxS; s++)
                {
                    cosdoy = cdoy[s - 1];
                    sindoy = sdoy[s - 1];
                    for (int n = l; n <= Constants.TMaxN; n++)
                    {
                        pl = plg[n, l];
                        bf[c] = pl * coslst * cosdoy;
                        bf[c + 1] = pl * sinlst * cosdoy;
                        bf[c + 2] = pl * coslst * sindoy;
                        bf[c + 3] = pl * sinlst * sindoy;
                        c += 4;
                    }
                }
            }

            // Stationary Planetary Waves (longitude dependence)
            if (c != Constants.CSpw) throw new InvalidOperationException("Problem with basis definitions");
            for (int m = 1; m <= Constants.PMaxM; m++)
            {
                coslon = clon[m - 1];
                sinlon = slon[m - 1];
                for (int n = m; n <= Constants.PMaxN; n++)
                {
                    pl = plg[n, m];
                    bf[c] = pl * coslon;
                    bf[c + 1] = pl * sinlon;
                    c += 2;
                }
                // Intra-annual modulation of SPWs
                for (int s = 1; s <= Constants.PMaxS; s++)
                {
                    cosdoy = cdoy[s - 1];
                    sindoy = sdoy[s - 1];
                    for (int n = m; n <= Constants.PMaxN; n++)
                    {
                        pl = plg[n, m];
                        bf[c] = pl * coslon * cosdoy;
                        bf[c + 1] = pl * sinlon * cosdoy;
                        bf[c + 2] = pl * coslon * sindoy;
                        bf[c + 3] = pl * sinlon * sindoy;
                        c += 4;
                    }
                }
            }

            // Linear solar flux terms
            if (c != Constants.CSfx) throw new InvalidOperationException("Problem with basis definitions");
            dfa = sfluxavg - sfluxAvgRef;
            df = sflux - sfluxavg;
            bf[c] = dfa;
            bf[c + 1] = dfa * dfa;
            bf[c + 2] = df;
            bf[c + 3] = df * df;
            bf[c + 4] = df * dfa;
            c += Constants.NSfx;

            // Additional linear terms
            if (c != Constants.CExtra) throw new InvalidOperationException("Problem with basis definitions");
            sza = SolZen(doy, lst, lat, lon);
            bf[c] = -0.5 * Math.Tanh((sza - 98.0) / 6.0);      // Solar zenith angle logistic function for O, H
            bf[c + 1] = -0.5 * Math.Tanh((sza - 101.5) / 20.0); // Solar zenith angle logistic function for NO
            bf[c + 2] = dfa * bf[c];                            // Solar flux modulation of logistic sza term
            bf[c + 3] = dfa * bf[c + 1];                        // Solar flux modulation of logistic sza term
            bf[c + 4] = dfa * plg[2, 0];                        // Solar flux modulation of P(2,0) term
            bf[c + 5] = dfa * plg[4, 0];                        // Solar flux modulation of P(4,0) term
            bf[c + 6] = dfa * plg[0, 0] * cdoy[0];              // Solar flux modulation of global AO
            bf[c + 7] = dfa * plg[0, 0] * sdoy[0];              // Solar flux modulation of global AO
            bf[c + 8] = dfa * plg[0, 0] * cdoy[1];              // Solar flux modulation of global SAO
            bf[c + 9] = dfa * plg[0, 0] * sdoy[1];              // Solar flux modulation of global SAO

            if (sfluxavg <= sfluxAvgQuadCutoff)
            {
                bf[c + 10] = dfa * dfa;
            }
            else
            {
                bf[c + 10] = (sfluxAvgQuadCutoff - sfluxAvgRef) * (2.0 * dfa - (sfluxAvgQuadCutoff - sfluxAvgRef));
            }
            bf[c + 11] = bf[c + 10] * plg[2, 0];                // P(2,0) modulation of truncated quadratic F10.7a term
            bf[c + 12] = bf[c + 10] * plg[4, 0];                // P(4,0) modulation of truncated quadratic F10.7a term
            bf[c + 13] = df * plg[2, 0];                        // P(2,0) modulation of df
            bf[c + 14] = df * plg[4, 0];                        // P(4,0) modulation of df

            //---------------------------------------------
            // Nonlinear Terms
            //---------------------------------------------

            c = Constants.CNonLin;

            // Solar flux modulation terms
            if (c != Constants.CSfxMod) throw new InvalidOperationException("Problem with basis definitions");
            bf[c] = dfa;
            bf[c + 1] = dfa * dfa;
            bf[c + 2] = df;
            bf[c + 3] = df * df;
            bf[c + 4] = df * dfa;
            c += Constants.NSfxMod;

            // Terms needed for legacy geomagnetic activity dependence
            if (c != Constants.CMag) throw new InvalidOperationException("Problem with basis set");
            for (int i = 0; i < 7; i++)
            {
                bf[c + i] = ap[i] - 4.0; // ap array is 1-based in Fortran, 0-based here
            }
            bf[c + 8] = Constants.Doy2Rad * doy;
            bf[c + 9] = Constants.Lst2Rad * lst;
            bf[c + 10] = Constants.Deg2Rad * lon;
            bf[c + 11] = Constants.Lst2Rad * utsec / 3600.0;
            bf[c + 12] = Math.Abs(lat);
            c += 13;
            for (int m = 0; m <= 1; m++)
            {
                for (int n = 0; n <= Constants.AMaxN; n++)
                {
                    bf[c] = plg[n, m];
                    c++;
                }
            }

            // Terms needed for legacy UT dependence
            c = Constants.CUt;
            bf[c] = Constants.Lst2Rad * utsec / 3600.0;
            bf[c + 1] = Constants.Doy2Rad * doy;
            bf[c + 2] = dfa;
            bf[c + 3] = Constants.Deg2Rad * lon;
            bf[c + 4] = plg[1, 0];
            bf[c + 5] = plg[3, 0];
            bf[c + 6] = plg[5, 0];
            bf[c + 7] = plg[3, 2];
            bf[c + 8] = plg[5, 2];

            //---------------------------------------------
            // Apply Switches
            //---------------------------------------------
            for (int i = 0; i <= Constants.Mbf; i++)
            {
                if (!Initialization.Swg[i]) bf[i] = 0.0;
            }
        }

        // ==================================================================================================
        // SOLZEN: Calculate solar zenith angle (adapted from IRI subroutine)
        // ==================================================================================================
        /// <summary>
        /// Calculate solar zenith angle
        /// </summary>
        private static double SolZen(double ddd, double lst, double lat, double lon)
        {
            const double Humr = Constants.Pi / 12.0;
            // const double Dumr = Constants.Pi / 182.5;
            double[] p = { 0.017203534, 0.034407068, 0.051610602, 0.068814136, 0.103221204 };

            double wlon = 360.0 - lon;
            double teqnx = ddd + (lst + wlon / 15.0) / 24.0 + 0.9369;
            teqnx = ddd + 0.9369;

            // Solar declination
            double dec = 23.256 * Math.Sin(p[0] * (teqnx - 82.242)) + 0.381 * Math.Sin(p[1] * (teqnx - 44.855))
                       + 0.167 * Math.Sin(p[2] * (teqnx - 23.355)) - 0.013 * Math.Sin(p[3] * (teqnx + 11.97))
                       + 0.011 * Math.Sin(p[4] * (teqnx - 10.410)) + 0.339137;
            dec = dec * Constants.Deg2Rad;

            // Equation of time
            double tf = teqnx - 0.5;
            double teqt = -7.38 * Math.Sin(p[0] * (tf - 4.0)) - 9.87 * Math.Sin(p[1] * (tf + 9.0))
                        + 0.27 * Math.Sin(p[2] * (tf - 53.0)) - 0.2 * Math.Cos(p[3] * (tf - 17.0));

            double phi = Humr * (lst - 12.0) + teqt * Constants.Deg2Rad / 4.0;
            double rlat = lat * Constants.Deg2Rad;

            // Cosine of solar zenith angle
            double cosx = Math.Sin(rlat) * Math.Sin(dec) + Math.Cos(rlat) * Math.Cos(dec) * Math.Cos(phi);
            if (Math.Abs(cosx) > 1.0) cosx = Math.Sign(cosx) * 1.0;

            return Math.Acos(cosx) / Constants.Deg2Rad;
        }

        // ==================================================================================================
        // SFLUXMOD: Legacy nonlinear modulation of intra-annual, tide, and SPW terms
        // ==================================================================================================
        /// <summary>
        /// Legacy nonlinear modulation of intra-annual, tide, and SPW terms
        /// </summary>
        public static double SFluxMod(int iz, double[] gf, BasisSubset parmset, double dffact)
        {
            double f1, f2, f3, sum;

            // Intra-annual modulation factor
            if (Initialization.Swg[Constants.CSfxMod])
            {
                f1 = parmset.Beta[Constants.CSfxMod, iz] * gf[Constants.CSfxMod]
                   + (parmset.Beta[Constants.CSfx + 2, iz] * gf[Constants.CSfxMod + 2]
                    + parmset.Beta[Constants.CSfx + 3, iz] * gf[Constants.CSfxMod + 3]) * dffact;
            }
            else
            {
                f1 = 0.0;
            }

            // Migrating tide (local time) modulation factor
            if (Initialization.Swg[Constants.CSfxMod + 1])
            {
                f2 = parmset.Beta[Constants.CSfxMod + 1, iz] * gf[Constants.CSfxMod]
                   + (parmset.Beta[Constants.CSfx + 2, iz] * gf[Constants.CSfxMod + 2]
                    + parmset.Beta[Constants.CSfx + 3, iz] * gf[Constants.CSfxMod + 3]) * dffact;
            }
            else
            {
                f2 = 0.0;
            }

            // SPW (longitude) modulation factor
            if (Initialization.Swg[Constants.CSfxMod + 2])
            {
                f3 = parmset.Beta[Constants.CSfxMod + 2, iz] * gf[Constants.CSfxMod];
            }
            else
            {
                f3 = 0.0;
            }

            sum = 0.0;
            for (int j = 0; j <= Constants.Mbf; j++)
            {
                // Apply intra-annual modulation
                if (Initialization.Zsfx[j])
                {
                    sum += parmset.Beta[j, iz] * gf[j] * f1;
                    continue;
                }
                // Apply migrating tide modulation
                if (Initialization.Tsfx[j])
                {
                    sum += parmset.Beta[j, iz] * gf[j] * f2;
                    continue;
                }
                // Apply SPW modulation
                if (Initialization.Psfx[j])
                {
                    sum += parmset.Beta[j, iz] * gf[j] * f3;
                    continue;
                }
            }

            return sum;
        }

        // ==================================================================================================
        // GEOMAG: Legacy nonlinear ap dependence (daily ap mode and ap history mode), including mixed
        //         ap/UT/Longitude terms.
        // Master switch control is as follows:
        //   !swg(cmag) && !swg(cmag+1)     Do nothing: Return zero
        //   swg(cmag) && swg(cmag+1)       Daily Ap mode
        //   swg(cmag) != swg(cmag+1)       3-hour ap history mode
        // ==================================================================================================
        /// <summary>
        /// Legacy nonlinear ap dependence
        /// </summary>
        public static double GeoMag(double[] p0, double[] bf, double[,] plg)
        {
            // Return zero if both master switches are off
            if (!(Initialization.Swg[Constants.CMag] || Initialization.Swg[Constants.CMag + 1]))
            {
                return 0.0;
            }

            // Copy parameters
            double[] p = (double[])p0.Clone();
            bool[] swg1 = new bool[Constants.NMag];
            Array.Copy(Initialization.Swg, Constants.CMag, swg1, 0, Constants.NMag);

            double geomag;

            // Calculate function
            if (swg1[0] == swg1[1])
            {
                // Daily Ap mode
                if (p[1] == 0) // If k00s is zero, then cannot compute function
                {
                    return 0.0;
                }

                // Apply switches
                for (int i = 2; i <= 25; i++)
                {
                    if (!swg1[i]) p[i] = 0.0;
                }
                p[8] = p0[8]; // Need doy phase term

                double delA = G0Fn(bf[0], p[0], p[1]);
                geomag = (p[2] * plg[0, 0] + p[3] * plg[2, 0] + p[4] * plg[4, 0]                               // time independent
                    + (p[5] * plg[1, 0] + p[6] * plg[3, 0] + p[7] * plg[5, 0]) * Math.Cos(bf[8] - p[8])        // doy modulation
                    + (p[9] * plg[1, 1] + p[10] * plg[3, 1] + p[11] * plg[5, 1]) * Math.Cos(bf[9] - p[12])     // local time modulation
                    + (1.0 + p[13] * plg[1, 0]) *
                      (p[14] * plg[2, 1] + p[15] * plg[4, 1] + p[16] * plg[6, 1]) * Math.Cos(bf[10] - p[17])   // longitude effect
                    + (p[18] * plg[1, 1] + p[19] * plg[3, 1] + p[20] * plg[5, 1]) * Math.Cos(bf[10] - p[21]) *
                      Math.Cos(bf[8] - p[8])                                                                    // longitude with doy modulation
                    + (p[22] * plg[1, 0] + p[23] * plg[3, 0] + p[24] * plg[5, 0]) * Math.Cos(bf[11] - p[25]))  // universal time
                    * delA;
            }
            else
            {
                // 3-hour ap history mode
                if (p[28] == 0) // If beta00 is zero, then cannot compute function
                {
                    return 0.0;
                }

                // Apply switches
                for (int i = 30; i < Constants.NMag; i++)
                {
                    if (!swg1[i]) p[i] = 0.0;
                }
                p[36] = p0[36]; // Need doy phase term

                double gbeta = p[28] / (1 + p[29] * (45.0 - bf[12]));
                double ex = Math.Exp(-10800.0 * gbeta);
                double sumex = 1 + (1 - Math.Pow(ex, 19.0)) * Math.Pow(ex, 0.5) / (1 - ex);
                double[] G = new double[7]; // G[1:6], index 0 unused
                for (int i = 1; i <= 6; i++)
                {
                    G[i] = G0Fn(bf[i], p[26], p[27]);
                }
                double delA = (G[1]
                            + (G[2] * ex + G[3] * ex * ex + G[4] * Math.Pow(ex, 3.0)
                             + (G[5] * Math.Pow(ex, 4.0) + G[6] * Math.Pow(ex, 12.0)) * (1 - Math.Pow(ex, 8.0)) / (1 - ex))) / sumex;

                geomag = (p[30] * plg[0, 0] + p[31] * plg[2, 0] + p[32] * plg[4, 0]                              // time independent
                    + (p[33] * plg[1, 0] + p[34] * plg[3, 0] + p[35] * plg[5, 0]) * Math.Cos(bf[8] - p[36])      // doy modulation
                    + (p[37] * plg[1, 1] + p[38] * plg[3, 1] + p[39] * plg[5, 1]) * Math.Cos(bf[9] - p[40])      // local time modulation
                    + (1.0 + p[41] * plg[1, 0]) *
                      (p[42] * plg[2, 1] + p[43] * plg[4, 1] + p[44] * plg[6, 1]) * Math.Cos(bf[10] - p[45])     // longitude effect
                    + (p[46] * plg[1, 1] + p[47] * plg[3, 1] + p[48] * plg[5, 1]) * Math.Cos(bf[10] - p[49]) *
                      Math.Cos(bf[8] - p[36])                                                                     // longitude with doy modulation
                    + (p[50] * plg[1, 0] + p[51] * plg[3, 0] + p[52] * plg[5, 0]) * Math.Cos(bf[11] - p[53]))    // universal time
                    * delA;
            }

            return geomag;
        }

        // Helper function for GeoMag
        private static double G0Fn(double a, double k00r, double k00s)
        {
            return a + (k00r - 1.0) * (a + (Math.Exp(-a * k00s) - 1.0) / k00s);
        }

        // ==================================================================================================
        // UTDEP: Legacy nonlinear UT dependence
        // ==================================================================================================
        /// <summary>
        /// Legacy nonlinear UT dependence
        /// </summary>
        public static double UtDep(double[] p0, double[] bf)
        {
            // Copy parameters
            double[] p = (double[])p0.Clone();
            bool[] swg1 = new bool[Constants.NUt];
            Array.Copy(Initialization.Swg, Constants.CUt, swg1, 0, Constants.NUt);

            // Apply switches
            for (int i = 3; i < Constants.NUt; i++)
            {
                if (!swg1[i]) p[i] = 0.0;
            }

            // Calculate function
            double utdep = Math.Cos(bf[0] - p[0]) *
                          (1 + p[3] * bf[4] * Math.Cos(bf[1] - p[1])) *
                          (1 + p[4] * bf[2]) * (1 + p[5] * bf[4]) *
                          (p[6] * bf[4] + p[7] * bf[5] + p[8] * bf[6]) +
                          Math.Cos(bf[0] - p[2] + 2 * bf[3]) * (p[9] * bf[7] + p[10] * bf[8]) * (1 + p[11] * bf[2]);

            return utdep;
        }
    }
}