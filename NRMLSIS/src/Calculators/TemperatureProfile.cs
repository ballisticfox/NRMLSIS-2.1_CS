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
using NRLMSIS.Infrastructure;
using static NRLMSIS.Infrastructure.BasisFunctionExtractor;

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
        /// Compute the vertical temperature and species-independent profile parameters.
        /// This computes spline coefficients, reference temperatures, and integration terms.
        /// </summary>
        /// <param name="basisFunctions">Array of horizontal and temporal basis function terms [0:maxnbf-1]</param>
        /// <param name="temperatureProfile">Output structure containing temperature vertical profile parameters</param>
        public static void TemperatureParameters(double[] basisFunctions, out TemperatureProfile temperatureProfile)
        {
            temperatureProfile = new TemperatureProfile();

            // Extract basis function subsets once (avoid duplication)
            double[] geomagPrimary = ExtractGeomagneticPrimary(basisFunctions);
            double[,] geomagSecondary = ExtractGeomagneticSecondary(basisFunctions);
            double[] utDependent = ExtractUniversalTimeDependent(basisFunctions);

            // Compute unconstrained spline coefficients
            ComputeUnconstrainedSplines(basisFunctions, temperatureProfile);

            // Compute exospheric temperature
            ComputeExosphericTemperature(basisFunctions, geomagPrimary, geomagSecondary,
                                        utDependent, temperatureProfile);

            // Compute temperature gradient at zetaB
            ComputeTemperatureGradientAtZetaB(basisFunctions, geomagPrimary, geomagSecondary,
                                             temperatureProfile);

            // Compute temperature at zetaB
            ComputeTemperatureAtZetaB(basisFunctions, geomagPrimary, geomagSecondary,
                                     temperatureProfile);

            // Compute shape factor
            temperatureProfile.Sigma = temperatureProfile.Tgb0 / (temperatureProfile.Tex - temperatureProfile.Tb0);

            // Apply C2 continuity constraint for top spline coefficients
            ApplyC2ContinuityConstraint(temperatureProfile);

            // Compute reference temperatures
            ComputeReferenceTemperatures(temperatureProfile);

            // Compute integration coefficients (Beta and Gamma)
            ComputeIntegrationCoefficients(temperatureProfile);

            // Compute integration constants and indefinite integrals
            ComputeIntegrationConstants(temperatureProfile);

            // Compute total number density at zetaF
            ComputeTotalDensityAtZetaF(temperatureProfile);
        }

        // ==================================================================================================
        // UNCONSTRAINED SPLINE COMPUTATION
        // ==================================================================================================

        /// <summary>
        /// Computes unconstrained spline coefficients for temperature profile below zetaB.
        /// Includes solar flux modulation for tidal components.
        /// </summary>
        private static void ComputeUnconstrainedSplines(double[] gf, TemperatureProfile tpro)
        {
            var subset = Initialization.TN;

            // Compute base coefficients
            for (int ix = 0; ix <= Constants.ItB0 - 1; ix++)
            {
                tpro.Cf[ix] = DotProduct(subset.Beta, 0, Constants.Mbf, ix, gf, 0, Constants.Mbf);
            }

            // Add solar flux modulation for tidal components
            for (int ix = 0; ix <= Constants.ItB0 - 1; ix++)
            {
                if (Initialization.Smod[ix])
                {
                    // SFluxMod adds F10.7 modulation of tides
                    double normalizationFactor = 1.0 / subset.Beta[0, ix - subset.Bl];
                    tpro.Cf[ix] += BasisFunctions.SFluxMod(ix, gf, subset, normalizationFactor);
                }
            }
        }

        // ==================================================================================================
        // BOUNDARY CONDITION COMPUTATION
        // ==================================================================================================

        /// <summary>
        /// Computes exospheric temperature (temperature at infinite altitude).
        /// Includes contributions from basis functions, solar flux, geomagnetic activity, and UT.
        /// </summary>
        private static void ComputeExosphericTemperature(double[] gf,
                                                         double[] geomagPrimary,
                                                         double[,] geomagSecondary,
                                                         double[] utDependent,
                                                         TemperatureProfile tpro)
        {
            var subset = Initialization.TN;
            int columnIndex = Constants.ItEx;

            // Base contribution
            tpro.Tex = DotProduct(subset.Beta, 0, Constants.Mbf, columnIndex, gf, 0, Constants.Mbf);

            // Solar flux modulation
            double normalizationFactor = 1.0 / subset.Beta[0, columnIndex - subset.Bl];
            tpro.Tex += BasisFunctions.SFluxMod(columnIndex, gf, subset, normalizationFactor);

            // Geomagnetic activity contribution
            double[] geomagParameters = ExtractGeomagneticParameters(subset, columnIndex);
            tpro.Tex += BasisFunctions.GeoMag(geomagParameters, geomagPrimary, geomagSecondary);

            // Universal time dependence
            double[] utParameters = ExtractUTParameters(subset, columnIndex);
            tpro.Tex += BasisFunctions.UtDep(utParameters, utDependent);
        }

        /// <summary>
        /// Computes temperature gradient at zetaB (122.5 km).
        /// Includes solar flux modulation and geomagnetic effects.
        /// </summary>
        private static void ComputeTemperatureGradientAtZetaB(double[] gf,
                                                              double[] geomagPrimary,
                                                              double[,] geomagSecondary,
                                                              TemperatureProfile tpro)
        {
            var subset = Initialization.TN;
            int columnIndex = Constants.ItGb0;

            // Base contribution
            tpro.Tgb0 = DotProduct(subset.Beta, 0, Constants.Mbf, columnIndex, gf, 0, Constants.Mbf);

            // Solar flux modulation (if enabled)
            if (Initialization.Smod[columnIndex])
            {
                double normalizationFactor = 1.0 / subset.Beta[0, columnIndex - subset.Bl];
                tpro.Tgb0 += BasisFunctions.SFluxMod(columnIndex, gf, subset, normalizationFactor);
            }

            // Geomagnetic activity contribution
            double[] geomagParameters = ExtractGeomagneticParameters(subset, columnIndex);
            tpro.Tgb0 += BasisFunctions.GeoMag(geomagParameters, geomagPrimary, geomagSecondary);
        }

        /// <summary>
        /// Computes temperature at zetaB (122.5 km).
        /// Includes solar flux modulation and geomagnetic effects.
        /// </summary>
        private static void ComputeTemperatureAtZetaB(double[] gf,
                                                      double[] geomagPrimary,
                                                      double[,] geomagSecondary,
                                                      TemperatureProfile tpro)
        {
            var subset = Initialization.TN;
            int columnIndex = Constants.ItB0;

            // Base contribution
            tpro.Tb0 = DotProduct(subset.Beta, 0, Constants.Mbf, columnIndex, gf, 0, Constants.Mbf);

            // Solar flux modulation (if enabled)
            if (Initialization.Smod[columnIndex])
            {
                double normalizationFactor = 1.0 / subset.Beta[0, columnIndex - subset.Bl];
                tpro.Tb0 += BasisFunctions.SFluxMod(columnIndex, gf, subset, normalizationFactor);
            }

            // Geomagnetic activity contribution
            double[] geomagParameters = ExtractGeomagneticParameters(subset, columnIndex);
            tpro.Tb0 += BasisFunctions.GeoMag(geomagParameters, geomagPrimary, geomagSecondary);
        }

        // ==================================================================================================
        // SPLINE CONSTRAINT APPLICATION
        // ==================================================================================================

        /// <summary>
        /// Applies C2 continuity constraint at zetaB for top three spline coefficients.
        /// Ensures smooth transition between spline and Bates profile regions.
        /// </summary>
        private static void ApplyC2ContinuityConstraint(TemperatureProfile tpro)
        {
            // Boundary condition vector
            double[] bc = new double[3];
            bc[0] = 1.0 / tpro.Tb0;
            bc[1] = -tpro.Tgb0 / (tpro.Tb0 * tpro.Tb0);
            bc[2] = -bc[1] * (tpro.Sigma + 2.0 * tpro.Tgb0 / tpro.Tb0);

            // Matrix multiplication: Cf[ItB0:ItEx] = bc * C2Tn
            for (int i = Constants.ItB0; i <= Constants.ItEx; i++)
            {
                int columnIndex = i - Constants.ItB0;
                tpro.Cf[i] = 0.0;
                for (int row = 0; row < 3; row++)
                {
                    tpro.Cf[i] += bc[row] * Constants.C2Tn[row, columnIndex];
                }
            }
        }

        // ==================================================================================================
        // REFERENCE TEMPERATURE COMPUTATION
        // ==================================================================================================

        /// <summary>
        /// Computes reference temperatures and temperature gradient at key altitudes.
        /// TZetaF at 70 km, TZetaA at 85 km, and gradient DlnTdzA at 85 km.
        /// </summary>
        private static void ComputeReferenceTemperatures(TemperatureProfile tpro)
        {
            // Temperature at zetaF (70 km)
            tpro.TZetaF = 1.0 / DotProduct(tpro.Cf, Constants.IzFx, Constants.IzFx + 2,
                                           Constants.S4ZetaF, 0, 2);

            // Temperature at zetaA (85 km)
            tpro.TZetaA = 1.0 / DotProduct(tpro.Cf, Constants.IzAx, Constants.IzAx + 2,
                                           Constants.S4ZetaA, 0, 2);

            // Logarithmic temperature gradient at zetaA
            tpro.DlnTdzA = -DotProduct(tpro.Cf, Constants.IzAx, Constants.IzAx + 2,
                                       Constants.WghtAxdz, 0, 2) * tpro.TZetaA;
        }

        // ==================================================================================================
        // INTEGRATION COEFFICIENT COMPUTATION
        // ==================================================================================================

        /// <summary>
        /// Computes Beta and Gamma integration coefficients for hydrostatic integration.
        /// Beta: first integral of spline basis functions (for V integral).
        /// Gamma: second integral of spline basis functions (for W integral).
        /// </summary>
        private static void ComputeIntegrationCoefficients(TemperatureProfile tpro)
        {
            // First integration coefficients (Beta) - cumulative sum
            tpro.Beta[0] = tpro.Cf[0] * Constants.WBeta[0];
            for (int ix = 1; ix <= Constants.Nl; ix++)
            {
                tpro.Beta[ix] = tpro.Beta[ix - 1] + tpro.Cf[ix] * Constants.WBeta[ix];
            }

            // Second integration coefficients (Gamma) - cumulative sum
            tpro.Gamma[0] = tpro.Beta[0] * Constants.WGamma[0];
            for (int ix = 1; ix <= Constants.Nl; ix++)
            {
                tpro.Gamma[ix] = tpro.Gamma[ix - 1] + tpro.Beta[ix] * Constants.WGamma[ix];
            }
        }

        // ==================================================================================================
        // INTEGRATION CONSTANT COMPUTATION
        // ==================================================================================================

        /// <summary>
        /// Computes integration constants and indefinite integrals at reference altitudes.
        /// These are used in hydrostatic equilibrium calculations for density profiles.
        /// </summary>
        private static void ComputeIntegrationConstants(TemperatureProfile tpro)
        {
            // Shape parameters
            tpro.B = 1.0 - tpro.Tb0 / tpro.Tex;
            tpro.SigmaSq = tpro.Sigma * tpro.Sigma;

            // First integration constants
            tpro.CVs = -DotProduct(tpro.Beta, Constants.ItB0 - 1, Constants.ItB0 + 2,
                                   Constants.S5ZetaB, 0, 3);
            tpro.CVb = -Math.Log(1.0 - tpro.B) / (tpro.Sigma * tpro.Tex);

            // Second integration constants
            tpro.CWs = -DotProduct(tpro.Gamma, Constants.ItB0 - 2, Constants.ItB0 + 2,
                                   Constants.S6ZetaB, 0, 4);
            tpro.CWb = -Utilities.Dilog(tpro.B) / (tpro.SigmaSq * tpro.Tex);

            // First indefinite integral (V) at reference altitudes
            tpro.VZetaF = DotProduct(tpro.Beta, Constants.IzFx - 1, Constants.IzFx + 2,
                                     Constants.S5ZetaF, 0, 3) + tpro.CVs;
            tpro.VZetaA = DotProduct(tpro.Beta, Constants.IzAx - 1, Constants.IzAx + 2,
                                     Constants.S5ZetaA, 0, 3) + tpro.CVs;
            tpro.VZeta0 = DotProduct(tpro.Beta, 0, 2,
                                     Constants.S5Zeta0, 0, 2) + tpro.CVs;

            // Second indefinite integral (W) at zetaA
            tpro.WZetaA = DotProduct(tpro.Gamma, Constants.IzAx - 2, Constants.IzAx + 2,
                                     Constants.S6ZetaA, 0, 4)
                        + tpro.CVs * (Constants.ZetaA - Constants.ZetaB)
                        + tpro.CWs;
        }

        // ==================================================================================================
        // TOTAL DENSITY COMPUTATION
        // ==================================================================================================

        /// <summary>
        /// Computes total atmospheric number density at zetaF (70 km).
        /// Uses hydrostatic equilibrium with mean molecular mass and reference pressure.
        /// </summary>
        private static void ComputeTotalDensityAtZetaF(TemperatureProfile tpro)
        {
            tpro.LnDTotF = Constants.LnP0
                         - Constants.MbarG0DivKB * (tpro.VZetaF - tpro.VZeta0)
                         - Math.Log(Constants.KB * tpro.TZetaF);
        }

        // ==================================================================================================
        // TFNX: Compute the temperature at specified geopotential height
        // ==================================================================================================

        /// <summary>
        /// Compute the temperature at specified geopotential height.
        /// Uses spline interpolation below zetaB (122.5 km) and Bates profile above.
        /// </summary>
        /// <param name="altitude">Geopotential height (km)</param>
        /// <param name="intervalIndex">B-spline reference index</param>
        /// <param name="splineWeights">B-spline weights - array[4] representing Fortran wght(-3:0)</param>
        /// <param name="temperatureProfile">Structure containing temperature vertical profile parameters</param>
        /// <returns>Temperature at height z (K)</returns>
        public static double TfnX(double altitude, int intervalIndex, double[] splineWeights,
                                 TemperatureProfile temperatureProfile)
        {
            if (altitude < Constants.ZetaB)
            {
                // Spline region: evaluate spline using precomputed weights
                return EvaluateSplineTemperature(intervalIndex, splineWeights, temperatureProfile);
            }
            else
            {
                // Bates profile region: exponential approach to exospheric temperature
                return EvaluateBatesTemperature(altitude, temperatureProfile);
            }
        }

        /// <summary>
        /// Evaluates temperature in spline region using B-spline interpolation.
        /// </summary>
        private static double EvaluateSplineTemperature(int iz, double[] wght, TemperatureProfile tpro)
        {
            // Determine valid coefficient range
            int coefficientStart = Math.Max(iz - 3, 0);

            // Determine weight starting index
            int weightStartOffset;
            if (iz < 3)
            {
                weightStartOffset = -iz;
            }
            else
            {
                weightStartOffset = -3;
            }

            // Compute weighted sum: dot_product(Cf[i:iz], wght[jStart:0])
            // Note: wght in Fortran is wght(-3:0), in C# wght[0:3]
            // wght(-3) maps to wght[0], wght(0) maps to wght[3]
            double inverseTemperature = 0.0;
            int weightIndex = weightStartOffset + 3; // Map Fortran index to C# index

            for (int cfIndex = coefficientStart; cfIndex <= iz; cfIndex++)
            {
                inverseTemperature += tpro.Cf[cfIndex] * wght[weightIndex];
                weightIndex++;
            }

            // Return temperature (reciprocal of sum)
            return 1.0 / inverseTemperature;
        }

        /// <summary>
        /// Evaluates temperature in Bates profile region.
        /// Exponential approach to exospheric temperature.
        /// </summary>
        private static double EvaluateBatesTemperature(double altitude, TemperatureProfile tpro)
        {
            double temperatureDeficit = tpro.Tex - tpro.Tb0;
            double exponentialDecay = Math.Exp(-tpro.Sigma * (altitude - Constants.ZetaB));

            return tpro.Tex - temperatureDeficit * exponentialDecay;
        }

        // ==================================================================================================
        // HELPER METHODS
        // ==================================================================================================

        /// <summary>
        /// Extracts UT-dependent parameters for a given column from the temperature subset.
        /// </summary>
        private static double[] ExtractUTParameters(BasisSubset subset, int column)
        {
            double[] parameters = new double[Constants.NUt];
            int columnIndex = column - subset.Bl;

            for (int i = 0; i < Constants.NUt; i++)
            {
                parameters[i] = subset.Beta[Constants.CUt + i, columnIndex];
            }

            return parameters;
        }

        /// <summary>
        /// Computes dot product of Beta subset array with basis function array.
        /// </summary>
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

        /// <summary>
        /// Computes dot product between two 1D arrays.
        /// </summary>
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