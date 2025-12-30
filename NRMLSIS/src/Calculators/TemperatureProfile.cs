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
// MSIS_Tfn: Contains vertical temperature profile parameters and subroutines, including
//           temperature integration terms.
//
// Key Translation Notes:
// - TemperatureProfile structure contains all temperature profile parameters
// - Array indexing: Fortran's Beta[0:maxnbf-1, bl:nl] becomes C# Beta[maxnbf, nl-bl+1]
// - Fortran array slicing in dot_product calls handled with helper methods
// - wght parameter in TfnX: Fortran wght(-3:0) maps to C# wght[0:3]
// **************************************************************************************************

using System;
using System.Linq;

namespace NRLMSIS.Calculators
{
    /// <summary>
    /// Temperature profile parameters structure
    /// </summary>
    public class TemperatureProfile
    {
        public double[] Cf { get; set; } = new double[Constants.Nl + 1];      // Spline coefficients [0:nl]
        public double TZetaF { get; set; }                                          // Tn at zetaF
        public double TZetaA { get; set; }                                          // Tn at zetaA (reference altitude for O1, H1)
        public double DlnTdzA { get; set; }                                         // log-temperature gradient at zetaA (km^-1)
        public double LnDTotF { get; set; }                                         // ln total number density at zetaF (m^-3)
        public double Tex { get; set; }                                             // Exospheric temperature
        public double Tgb0 { get; set; }                                            // Temperature gradient at zetaB
        public double Tb0 { get; set; }                                             // Temperature at zetaB
        public double Sigma { get; set; }                                           // Shape factor
        public double SigmaSq { get; set; }                                         // Sigma squared
        public double B { get; set; }                                               // b = 1-tb0/tex
        public double[] Beta { get; set; } = new double[Constants.Nl + 1];     // 1st integration coefficients on k=5 splines [0:nl]
        public double[] Gamma { get; set; } = new double[Constants.Nl + 1];    // 2nd integration coefficients on k=6 splines [0:nl]
        public double CVs { get; set; }                                             // 1st integration constant (spline portion)
        public double CVb { get; set; }                                             // 1st integration constant (Bates portion)
        public double CWs { get; set; }                                             // 2nd integration constant (spline portion)
        public double CWb { get; set; }                                             // 2nd integration constant (Bates portion)
        public double VZetaF { get; set; }                                          // 1st indefinite integral at zetaF
        public double VZetaA { get; set; }                                          // 1st indefinite integral at zetaA
        public double WZetaA { get; set; }                                          // 2nd indefinite integral at zetaA
        public double VZeta0 { get; set; }                                          // 1st indefinite integral at zeta=0
    }

    /// <summary>
    /// Temperature and species-independent profile functions
    /// </summary>
    public static class TemperatureFunction
    {
        // ==================================================================================================
        // TemperatureParameters: Compute the vertical temperature and species-independent profile parameters
        // ==================================================================================================
        /// <summary>
        /// Compute the vertical temperature and species-independent profile parameters
        /// </summary>
        /// <param name="gf">Array of horizontal and temporal basis function terms [0:maxnbf-1]</param>
        /// <param name="tpro">Output structure containing temperature vertical profile parameters</param>
        public static void TemperatureParameters(double[] gf, out TemperatureProfile tpro)
        {
            tpro = new TemperatureProfile();

            // Unconstrained spline coefficients
            for (int ix = 0; ix <= Constants.ItB0 - 1; ix++)
            {
                tpro.Cf[ix] = DotProduct(Initialization.TN.Beta, 0, Constants.Mbf, ix, gf, 0, Constants.Mbf);
            }

            for (int ix = 0; ix <= Constants.ItB0 - 1; ix++)
            {
                if (Initialization.Smod[ix])
                {
                    // sfluxmod adds F10.7 modulation of tides
                    tpro.Cf[ix] = tpro.Cf[ix] + BasisFunctions.SFluxMod(ix, gf, Initialization.TN, 1.0 / Initialization.TN.Beta[0, ix - Initialization.TN.Bl]);
                }
            }

            // Exospheric temperature
            tpro.Tex = DotProduct(Initialization.TN.Beta, 0, Constants.Mbf, Constants.ItEx, gf, 0, Constants.Mbf);
            tpro.Tex = tpro.Tex + BasisFunctions.SFluxMod(Constants.ItEx, gf, Initialization.TN, 1.0 / Initialization.TN.Beta[0, Constants.ItEx - Initialization.TN.Bl]);

            // Extract geomag parameters and basis functions
            double[] geomagParms = new double[Constants.NMag];
            for (int i = 0; i < Constants.NMag; i++)
            {
                geomagParms[i] = Initialization.TN.Beta[Constants.CMag + i, Constants.ItEx - Initialization.TN.Bl];
            }
            double[] geomagBf1 = new double[13];
            Array.Copy(gf, Constants.CMag, geomagBf1, 0, 13);
            // Extract geomag basis functions - need to construct plg array from gf
            double[,] geomagBf2 = new double[7, 2]; // plg(0:6, 0:1)
            int gfOffset = Constants.CMag + 13;
            for (int n = 0; n <= 6; n++)
            {
                geomagBf2[n, 0] = gf[gfOffset + n];
            }
            for (int n = 0; n <= 6; n++)
            {
                geomagBf2[n, 1] = gf[gfOffset + 7 + n];
            }
            tpro.Tex = tpro.Tex + BasisFunctions.GeoMag(geomagParms, geomagBf1, geomagBf2);

            // Extract utdep parameters and basis functions
            double[] utdepParms = new double[Constants.NUt];
            for (int i = 0; i < Constants.NUt; i++)
            {
                utdepParms[i] = Initialization.TN.Beta[Constants.CUt + i, Constants.ItEx - Initialization.TN.Bl];
            }
            double[] utdepBf = new double[9];
            Array.Copy(gf, Constants.CUt, utdepBf, 0, 9);
            tpro.Tex = tpro.Tex + BasisFunctions.UtDep(utdepParms, utdepBf);

            // Temperature gradient at zetaB (122.5 km)
            tpro.Tgb0 = DotProduct(Initialization.TN.Beta, 0, Constants.Mbf, Constants.ItGb0, gf, 0, Constants.Mbf);
            if (Initialization.Smod[Constants.ItGb0])
            {
                tpro.Tgb0 = tpro.Tgb0 + BasisFunctions.SFluxMod(Constants.ItGb0, gf, Initialization.TN, 1.0 / Initialization.TN.Beta[0, Constants.ItGb0 - Initialization.TN.Bl]);
            }

            // Extract geomag parameters for Tgb0
            for (int i = 0; i < Constants.NMag; i++)
            {
                geomagParms[i] = Initialization.TN.Beta[Constants.CMag + i, Constants.ItGb0 - Initialization.TN.Bl];
            }
            tpro.Tgb0 = tpro.Tgb0 + BasisFunctions.GeoMag(geomagParms, geomagBf1, geomagBf2);

            // Temperature at zetaB (122.5 km)
            tpro.Tb0 = DotProduct(Initialization.TN.Beta, 0, Constants.Mbf, Constants.ItB0, gf, 0, Constants.Mbf);
            if (Initialization.Smod[Constants.ItB0])
            {
                tpro.Tb0 = tpro.Tb0 + BasisFunctions.SFluxMod(Constants.ItB0, gf, Initialization.TN, 1.0 / Initialization.TN.Beta[0, Constants.ItB0 - Initialization.TN.Bl]);
            }

            // Extract geomag parameters for Tb0
            for (int i = 0; i < Constants.NMag; i++)
            {
                geomagParms[i] = Initialization.TN.Beta[Constants.CMag + i, Constants.ItB0 - Initialization.TN.Bl];
            }
            tpro.Tb0 = tpro.Tb0 + BasisFunctions.GeoMag(geomagParms, geomagBf1, geomagBf2);

            // Shape factor
            tpro.Sigma = tpro.Tgb0 / (tpro.Tex - tpro.Tb0);

            // Constrain top three spline coefficients for C2 continuity
            double[] bc = new double[3];
            bc[0] = 1.0 / tpro.Tb0;
            bc[1] = -tpro.Tgb0 / (tpro.Tb0 * tpro.Tb0);
            bc[2] = -bc[1] * (tpro.Sigma + 2.0 * tpro.Tgb0 / tpro.Tb0);

            // Matrix multiplication: bc * c2tn
            for (int i = Constants.ItB0; i <= Constants.ItEx; i++)
            {
                int col = i - Constants.ItB0;
                tpro.Cf[i] = 0.0;
                for (int row = 0; row < 3; row++)
                {
                    tpro.Cf[i] += bc[row] * Constants.C2Tn[row, col];
                }
            }

            // Reference temperature at zetaF (70 km)
            tpro.TZetaF = 1.0 / DotProduct(tpro.Cf, Constants.IzFx, Constants.IzFx + 2,
                                           Constants.S4ZetaF, 0, 2);

            // Reference temperature and gradient at zetaA (85 km)
            tpro.TZetaA = 1.0 / DotProduct(tpro.Cf, Constants.IzAx, Constants.IzAx + 2,
                                           Constants.S4ZetaA, 0, 2);
            tpro.DlnTdzA = -DotProduct(tpro.Cf, Constants.IzAx, Constants.IzAx + 2,
                                       Constants.WghtAxdz, 0, 2) * tpro.TZetaA;

            // Calculate spline coefficients for first and second 1/T integrals
            tpro.Beta[0] = tpro.Cf[0] * Constants.WBeta[0];
            for (int ix = 1; ix <= Constants.Nl; ix++)
            {
                tpro.Beta[ix] = tpro.Beta[ix - 1] + tpro.Cf[ix] * Constants.WBeta[ix];
            }

            tpro.Gamma[0] = tpro.Beta[0] * Constants.WGamma[0];
            for (int ix = 1; ix <= Constants.Nl; ix++)
            {
                tpro.Gamma[ix] = tpro.Gamma[ix - 1] + tpro.Beta[ix] * Constants.WGamma[ix];
            }

            // Integration terms and constants
            tpro.B = 1 - tpro.Tb0 / tpro.Tex;
            tpro.SigmaSq = tpro.Sigma * tpro.Sigma;
            tpro.CVs = -DotProduct(tpro.Beta, Constants.ItB0 - 1, Constants.ItB0 + 2,
                                   Constants.S5ZetaB, 0, 3);
            tpro.CWs = -DotProduct(tpro.Gamma, Constants.ItB0 - 2, Constants.ItB0 + 2,
                                   Constants.S6ZetaB, 0, 4);
            tpro.CVb = -Math.Log(1 - tpro.B) / (tpro.Sigma * tpro.Tex);
            tpro.CWb = -Utilities.Dilog(tpro.B) / (tpro.SigmaSq * tpro.Tex);
            tpro.VZetaF = DotProduct(tpro.Beta, Constants.IzFx - 1, Constants.IzFx + 2,
                                     Constants.S5ZetaF, 0, 3) + tpro.CVs;
            tpro.VZetaA = DotProduct(tpro.Beta, Constants.IzAx - 1, Constants.IzAx + 2,
                                     Constants.S5ZetaA, 0, 3) + tpro.CVs;
            tpro.WZetaA = DotProduct(tpro.Gamma, Constants.IzAx - 2, Constants.IzAx + 2,
                                     Constants.S6ZetaA, 0, 4) + tpro.CVs * (Constants.ZetaA - Constants.ZetaB) + tpro.CWs;
            tpro.VZeta0 = DotProduct(tpro.Beta, 0, 2,
                                     Constants.S5Zeta0, 0, 2) + tpro.CVs;

            // Compute total number density at zetaF
            tpro.LnDTotF = Constants.LnP0 - Constants.MbarG0DivKB * (tpro.VZetaF - tpro.VZeta0)
                         - Math.Log(Constants.KB * tpro.TZetaF);
        }

        // ==================================================================================================
        // TFNX: Compute the temperature at specified geopotential height
        // ==================================================================================================
        /// <summary>
        /// Compute the temperature at specified geopotential height
        /// </summary>
        /// <param name="z">Geopotential height (km)</param>
        /// <param name="iz">Bspline reference index</param>
        /// <param name="wght">Bspline weights - array[4] representing Fortran wght(-3:0)</param>
        /// <param name="tpro">Structure containing temperature vertical profile parameters</param>
        /// <returns>Temperature at height z</returns>
        public static double TfnX(double z, int iz, double[] wght, TemperatureProfile tpro)
        {
            double tfnx;

            if (z < Constants.ZetaB)
            {
                // Spline region
                int i = Math.Max(iz - 3, 0);
                int jStart;

                if (iz < 3)
                {
                    jStart = -iz;
                }
                else
                {
                    jStart = -3;
                }

                // Fortran: dot_product(tpro%cf(i:iz), wght(jStart:0))
                // wght in Fortran is wght(-3:0), in C# wght[0:3]
                // wght(-3) maps to wght[0], wght(0) maps to wght[3]
                // jStart ranges from -3 to 0
                tfnx = 0.0;
                int wghtIdx = jStart + 3; // Map Fortran index to C# index
                for (int cfIdx = i; cfIdx <= iz; cfIdx++)
                {
                    tfnx += tpro.Cf[cfIdx] * wght[wghtIdx];
                    wghtIdx++;
                }
                tfnx = 1.0 / tfnx;
            }
            else
            {
                // Bates profile region
                tfnx = tpro.Tex - (tpro.Tex - tpro.Tb0) * Math.Exp(-tpro.Sigma * (z - Constants.ZetaB));
            }

            return tfnx;
        }

        // Helper method for dot product between subset Beta array and gf array
        private static double DotProduct(double[,] beta, int betaRowStart, int betaRowEnd, int betaCol,
                                        double[] gf, int gfStart, int gfEnd)
        {
            double sum = 0.0;
            int gfIdx = gfStart;
            for (int i = betaRowStart; i <= betaRowEnd; i++)
            {
                sum += beta[i, betaCol - Initialization.TN.Bl] * gf[gfIdx];
                gfIdx++;
            }
            return sum;
        }

        // Helper method for dot product between two arrays
        private static double DotProduct(double[] a, int aStart, int aEnd, double[] b, int bStart, int bEnd)
        {
            double sum = 0.0;
            int bIdx = bStart;
            for (int i = aStart; i <= aEnd; i++)
            {
                sum += a[i] * b[bIdx];
                bIdx++;
            }
            return sum;
        }
    }
}