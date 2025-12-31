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
// MSIS_Dfn: Contains vertical species density profile parameters and subroutines
//
// Key Translation Notes:
// - DensityParameters structure contains all density profile parameters for each species
// - Array indexing: Fortran arrays with negative indices mapped to C# 0-based arrays
// - Species indices: 2=N2, 3=O2, 4=O, 5=He, 6=H, 7=Ar, 8=N, 9=Anomalous O, 10=NO
// **************************************************************************************************

using System;
using NRLMSIS.Infrastructure;
using static NRLMSIS.Infrastructure.DensityProfileColumns;
using static NRLMSIS.Infrastructure.BasisFunctionExtractor;

namespace NRLMSIS.Calculators
{
    /// <summary>
    /// Density profile parameters structure containing all vertical profile parameters for a species
    /// </summary>
    public class DensityParameters
    {
        public double LnPhiF { get; set; }                                        // Natural log of mixing ratio at zetaF
        public double LnDRef { get; set; }                                        // Natural log of number density at reference height
        public double ZetaM { get; set; }                                         // Turbopause height
        public double HML { get; set; }                                           // Scale height of lower portion of effective mass profile
        public double HMU { get; set; }                                           // Scale height of upper portion of effective mass profile
        public double C { get; set; }                                             // Chapman term coefficient
        public double ZetaC { get; set; }                                         // Chapman term reference height
        public double HC { get; set; }                                            // Chapman term scale height
        public double R { get; set; }                                             // Chemical/dynamical term coefficient
        public double ZetaR { get; set; }                                         // Chemical/dynamical term reference height
        public double HR { get; set; }                                            // Chemical/dynamical term scale height
        public double[] Cf { get; set; } = new double[Constants.NsplO1 + 2]; // Merged spline coefficients [0:nsplO1+1]
        public double ZRef { get; set; }                                          // Reference height for hydrostatic integral
        public double[] Mi { get; set; } = new double[5];                         // Effective mass at nodes [0:4]
        public double[] ZetaMi { get; set; } = new double[5];                     // Height of nodes [0:4]
        public double[] AMi { get; set; } = new double[5];                        // Slopes of piecewise mass profile [0:4]
        public double[] WMi { get; set; } = new double[5];                        // 2nd indefinite integral at nodes [0:4]
        public double[] XMi { get; set; } = new double[5];                        // Cumulative adjustment [0:4]
        public double IzRef { get; set; }                                         // Indefinite hydrostatic integral at reference height
        public double TRef { get; set; }                                          // Temperature at reference height
        public double ZMin { get; set; }                                          // Minimum height of profile
        public double ZHyd { get; set; }                                          // Hydrostatic terms needed above this height
        public int ISpec { get; set; }                                            // Species index
    }

    /// <summary>
    /// Density profile functions for computing vertical species distributions
    /// </summary>
    public static class DensityProfile
    {
        // ==================================================================================================
        // DensityParameters: Compute the species density profile parameters
        // ==================================================================================================
        /// <summary>
        /// Compute the species density profile parameters for a given species.
        /// This extracts parameters from the basis function coefficients and temperature profile.
        /// </summary>
        /// <param name="speciesIndex">Species index (2=N2, 3=O2, 4=O, 5=He, 6=H, 7=Ar, 8=N, 9=AnomalousO, 10=NO)</param>
        /// <param name="basisFunctions">Array of horizontal and temporal basis function terms [0:maxnbf-1]</param>
        /// <param name="temperatureProfile">Structure containing temperature vertical profile parameters</param>
        /// <param name="densityProfile">Output: density vertical profile parameters</param>
        public static void DensityParameters(int speciesIndex, double[] basisFunctions,
                                            TemperatureProfile temperatureProfile,
                                            out DensityParameters densityProfile)
        {
            densityProfile = new DensityParameters();
            densityProfile.ISpec = speciesIndex;

            // Extract basis function subsets used across multiple species
            double[] geomagPrimary = ExtractGeomagneticPrimary(basisFunctions);
            double[,] geomagSecondary = ExtractGeomagneticSecondary(basisFunctions);
            double[] utDependent = ExtractUniversalTimeDependent(basisFunctions);

            // Compute parameters based on species type
            switch (speciesIndex)
            {
                case 2: // Molecular Nitrogen (N2)
                    ComputeN2Parameters(densityProfile, basisFunctions, temperatureProfile,
                                       geomagPrimary, geomagSecondary, utDependent);
                    break;

                case 3: // Molecular Oxygen (O2)
                    ComputeO2Parameters(densityProfile, basisFunctions, temperatureProfile,
                                       geomagPrimary, geomagSecondary, utDependent);
                    break;

                case 4: // Atomic Oxygen (O)
                    ComputeO1Parameters(densityProfile, basisFunctions, temperatureProfile,
                                       geomagPrimary, geomagSecondary, utDependent);
                    break;

                case 5: // Helium (He)
                    ComputeHeParameters(densityProfile, basisFunctions, temperatureProfile,
                                       geomagPrimary, geomagSecondary, utDependent);
                    break;

                case 6: // Atomic Hydrogen (H)
                    ComputeH1Parameters(densityProfile, basisFunctions, temperatureProfile,
                                       geomagPrimary, geomagSecondary, utDependent);
                    break;

                case 7: // Argon (Ar)
                    ComputeArParameters(densityProfile, basisFunctions, temperatureProfile,
                                       geomagPrimary, geomagSecondary, utDependent);
                    break;

                case 8: // Atomic Nitrogen (N)
                    ComputeN1Parameters(densityProfile, basisFunctions, temperatureProfile,
                                       geomagPrimary, geomagSecondary, utDependent);
                    break;

                case 9: // Anomalous Oxygen (O*)
                    ComputeOAParameters(densityProfile, basisFunctions, temperatureProfile,
                                       geomagPrimary, geomagSecondary, utDependent);
                    break;

                case 10: // Nitric Oxide (NO)
                    ComputeNOParameters(densityProfile, basisFunctions, temperatureProfile,
                                       geomagPrimary, geomagSecondary, utDependent);
                    break;

                default:
                    throw new ArgumentException($"Invalid species index: {speciesIndex}");
            }

            // Compute piecewise mass profile and integration terms
            ComputeMassProfileAndIntegrals(speciesIndex, basisFunctions, temperatureProfile, densityProfile);
        }

        // ==================================================================================================
        // POST-PROCESSING: Mass Profile and Integration Terms
        // ==================================================================================================

        /// <summary>
        /// Computes piecewise mass profile nodes, slopes, and integration terms.
        /// This must be called after species-specific parameters are set.
        /// </summary>
        private static void ComputeMassProfileAndIntegrals(int speciesIndex, double[] basisFunctions,
                                                          TemperatureProfile tpro, DensityParameters dpro)
        {
            // Compute piecewise mass profile node heights
            dpro.ZetaMi[0] = dpro.ZetaM - 2.0 * dpro.HML;
            dpro.ZetaMi[1] = dpro.ZetaM - dpro.HML;
            dpro.ZetaMi[2] = dpro.ZetaM;
            dpro.ZetaMi[3] = dpro.ZetaM + dpro.HMU;
            dpro.ZetaMi[4] = dpro.ZetaM + 2.0 * dpro.HMU;

            // Compute effective masses at nodes
            dpro.Mi[0] = Constants.Mbar;
            dpro.Mi[4] = Constants.SpecMass[speciesIndex];
            dpro.Mi[2] = (dpro.Mi[0] + dpro.Mi[4]) / 2.0;
            double delM = Constants.Tanh1 * (dpro.Mi[4] - dpro.Mi[0]) / 2.0;
            dpro.Mi[1] = dpro.Mi[2] - delM;
            dpro.Mi[3] = dpro.Mi[2] + delM;

            // Compute slopes between nodes
            for (int i = 0; i <= 3; i++)
            {
                dpro.AMi[i] = (dpro.Mi[i + 1] - dpro.Mi[i]) / (dpro.ZetaMi[i + 1] - dpro.ZetaMi[i]);
            }

            // Compute second indefinite integral W at nodes
            for (int i = 0; i <= 4; i++)
            {
                double delz = dpro.ZetaMi[i] - Constants.ZetaB;
                if (dpro.ZetaMi[i] < Constants.ZetaB)
                {
                    // Below turbopause: use B-spline interpolation
                    Utilities.BSpline(dpro.ZetaMi[i], Constants.NodesTN, Constants.Nd + 2, 6,
                                    Initialization.EtaTN, out double[,] Si, out int iz);

                    // Extract Si[-5:0, 6] -> Si[0:5, 4]
                    double Wi = 0.0;
                    for (int j = 0; j <= 5; j++)
                    {
                        int gammaIndex = iz - 5 + j;
                        if (gammaIndex >= 0 && gammaIndex <= Constants.Nl)
                        {
                            Wi += tpro.Gamma[gammaIndex] * Si[j, 4];
                        }
                    }
                    dpro.WMi[i] = Wi + tpro.CVs * delz + tpro.CWs;
                }
                else
                {
                    // Above turbopause: analytic formula
                    dpro.WMi[i] = (0.5 * delz * delz + Utilities.Dilog(tpro.B * Math.Exp(-tpro.Sigma * delz)) / tpro.SigmaSq) / tpro.Tex
                                + tpro.CVb * delz + tpro.CWb;
                }
            }

            // Compute cumulative adjustment terms X at nodes
            dpro.XMi[0] = -dpro.AMi[0] * dpro.WMi[0];
            for (int i = 1; i <= 3; i++)
            {
                dpro.XMi[i] = dpro.XMi[i - 1] - dpro.WMi[i] * (dpro.AMi[i] - dpro.AMi[i - 1]);
            }
            dpro.XMi[4] = dpro.XMi[3] + dpro.WMi[4] * dpro.AMi[3];

            // Calculate hydrostatic integral at reference height and copy temperature
            double Mzref;
            if (dpro.ZRef == Constants.ZetaF)
            {
                Mzref = Constants.Mbar;
                dpro.TRef = tpro.TZetaF;
                dpro.IzRef = Constants.Mbar * tpro.VZetaF;
            }
            else if (dpro.ZRef == Constants.ZetaB)
            {
                Mzref = ComputePiecewiseMass(dpro.ZRef, dpro.ZetaMi, dpro.Mi, dpro.AMi);
                dpro.TRef = tpro.Tb0;
                dpro.IzRef = 0.0;

                if ((Constants.ZetaB > dpro.ZetaMi[0]) && (Constants.ZetaB < dpro.ZetaMi[4]))
                {
                    int nodeIndex = FindNodeIndex(Constants.ZetaB, dpro.ZetaMi);
                    dpro.IzRef = dpro.IzRef - dpro.XMi[nodeIndex];
                }
                else
                {
                    dpro.IzRef = dpro.IzRef - dpro.XMi[4];
                }
            }
            else if (dpro.ZRef == Constants.ZetaA)
            {
                Mzref = ComputePiecewiseMass(dpro.ZRef, dpro.ZetaMi, dpro.Mi, dpro.AMi);
                dpro.TRef = tpro.TZetaA;
                dpro.IzRef = Mzref * tpro.VZetaA;

                if ((Constants.ZetaA > dpro.ZetaMi[0]) && (Constants.ZetaA < dpro.ZetaMi[4]))
                {
                    int nodeIndex = FindNodeIndex(Constants.ZetaA, dpro.ZetaMi);
                    dpro.IzRef = dpro.IzRef - (dpro.AMi[nodeIndex] * tpro.WZetaA + dpro.XMi[nodeIndex]);
                }
                else
                {
                    dpro.IzRef = dpro.IzRef - dpro.XMi[4];
                }
            }
            else
            {
                throw new InvalidOperationException("Integrals at reference height not available");
            }

            // C1 constraint for O1 at 85 km
            if (speciesIndex == 4)
            {
                ApplyO1Constraint(basisFunctions, tpro, dpro, Mzref);
            }

            // C1 constraint for NO at 122.5 km
            if (speciesIndex == 10)
            {
                ApplyNOConstraint(basisFunctions, tpro, dpro, Mzref);
            }
        }

        /// <summary>
        /// Applies C1 continuity constraint for atomic oxygen at 85 km
        /// </summary>
        private static void ApplyO1Constraint(double[] gf, TemperatureProfile tpro, DensityParameters dpro, double Mzref)
        {
            double Cterm = dpro.C * Math.Exp(-(dpro.ZRef - dpro.ZetaC) / dpro.HC);
            double Rterm0 = Math.Tanh((dpro.ZRef - dpro.ZetaR) / (Initialization.HRFactO1Ref * dpro.HR));
            double Rterm = dpro.R * (1.0 + Rterm0);

            double[] bc = new double[2];
            bc[0] = dpro.LnDRef - Cterm + Rterm - dpro.Cf[7] * Constants.C1O1Adj[0];
            bc[1] = -Mzref * Constants.G0DivKB / tpro.TZetaA
                   - tpro.DlnTdzA
                   + Cterm / dpro.HC
                   + Rterm * (1.0 - Rterm0) / dpro.HR * Initialization.DHRFactO1Ref
                   - dpro.Cf[7] * Constants.C1O1Adj[1];

            // Compute coefficients for constrained splines: bc * C1O1
            for (int idx = 8; idx <= 9; idx++)
            {
                dpro.Cf[idx] = 0.0;
                for (int row = 0; row < 2; row++)
                {
                    dpro.Cf[idx] += bc[row] * Constants.C1O1[row, idx - 8];
                }
            }
        }

        /// <summary>
        /// Applies C1 continuity constraint for nitric oxide at 122.5 km
        /// </summary>
        private static void ApplyNOConstraint(double[] gf, TemperatureProfile tpro, DensityParameters dpro, double Mzref)
        {
            double Cterm = dpro.C * Math.Exp(-(dpro.ZRef - dpro.ZetaC) / dpro.HC);
            double Rterm0 = Math.Tanh((dpro.ZRef - dpro.ZetaR) / (Initialization.HRFactNORef * dpro.HR));
            double Rterm = dpro.R * (1.0 + Rterm0);

            double[] bc = new double[2];
            bc[0] = dpro.LnDRef - Cterm + Rterm - dpro.Cf[7] * Constants.C1NOAdj[0];
            bc[1] = -Mzref * Constants.G0DivKB / tpro.Tb0
                   - tpro.Tgb0 / tpro.Tb0
                   + Cterm / dpro.HC
                   + Rterm * (1.0 - Rterm0) / dpro.HR * Initialization.DHRFactNORef
                   - dpro.Cf[7] * Constants.C1NOAdj[1];

            // Compute coefficients for constrained splines: bc * C1NO
            for (int idx = 8; idx <= 9; idx++)
            {
                dpro.Cf[idx] = 0.0;
                for (int row = 0; row < 2; row++)
                {
                    dpro.Cf[idx] += bc[row] * Constants.C1NO[row, idx - 8];
                }
            }
        }

        // ==================================================================================================
        // SPECIES-SPECIFIC PARAMETER COMPUTATION
        // ==================================================================================================

        /// <summary>
        /// Computes density profile parameters for Molecular Nitrogen (N2)
        /// </summary>
        private static void ComputeN2Parameters(DensityParameters dpro, double[] gf,
                                                TemperatureProfile tpro,
                                                double[] geomagPrimary, double[,] geomagSecondary,
                                                double[] utDependent)
        {
            var subset = Initialization.N2;
            int bl = subset.Bl;

            // Reference parameters
            dpro.LnPhiF = Constants.LnVmr[2];
            dpro.LnDRef = tpro.LnDTotF + dpro.LnPhiF;
            dpro.ZRef = Constants.ZetaF;
            dpro.ZMin = -1.0;
            dpro.ZHyd = Constants.ZetaF;

            // Turbopause and scale heights
            dpro.ZetaM = DotProduct(subset.Beta, 0, Constants.Mbf, TurbopauseHeight, gf, 0, Constants.Mbf);
            dpro.HML = subset.Beta[0, LowerScaleHeight - bl];
            dpro.HMU = subset.Beta[0, UpperScaleHeight - bl];

            // Chemical/dynamical term
            dpro.R = 0.0;
            if (Initialization.N2RFlag)
            {
                dpro.R = DotProduct(subset.Beta, 0, Constants.Mbf, ChemicalDynamicCoefficient, gf, 0, Constants.Mbf);
            }
            dpro.ZetaR = subset.Beta[0, ChemicalDynamicReferenceHeight - bl];
            dpro.HR = subset.Beta[0, ChemicalDynamicScaleHeight - bl];
        }

        /// <summary>
        /// Computes density profile parameters for Molecular Oxygen (O2)
        /// </summary>
        private static void ComputeO2Parameters(DensityParameters dpro, double[] gf,
                                                TemperatureProfile tpro,
                                                double[] geomagPrimary, double[,] geomagSecondary,
                                                double[] utDependent)
        {
            var subset = Initialization.O2;
            int bl = subset.Bl;

            // Reference parameters
            dpro.LnPhiF = Constants.LnVmr[3];
            dpro.LnDRef = tpro.LnDTotF + dpro.LnPhiF;
            dpro.ZRef = Constants.ZetaF;
            dpro.ZMin = -1.0;
            dpro.ZHyd = Constants.ZetaF;

            // Turbopause and scale heights
            dpro.ZetaM = subset.Beta[0, TurbopauseHeight - bl];
            dpro.HML = subset.Beta[0, LowerScaleHeight - bl];
            dpro.HMU = subset.Beta[0, UpperScaleHeight - bl];

            // Chemical/dynamical term (includes geomagnetic effect)
            dpro.R = DotProduct(subset.Beta, 0, Constants.Mbf, ChemicalDynamicCoefficient, gf, 0, Constants.Mbf);
            dpro.R += BasisFunctions.GeoMag(
                ExtractGeomagneticParameters(subset, ChemicalDynamicCoefficient),
                geomagPrimary, geomagSecondary);

            dpro.ZetaR = subset.Beta[0, ChemicalDynamicReferenceHeight - bl];
            dpro.HR = subset.Beta[0, ChemicalDynamicScaleHeight - bl];
        }

        /// <summary>
        /// Computes density profile parameters for Atomic Oxygen (O)
        /// </summary>
        private static void ComputeO1Parameters(DensityParameters dpro, double[] gf,
                                                TemperatureProfile tpro,
                                                double[] geomagPrimary, double[,] geomagSecondary,
                                                double[] utDependent)
        {
            var subset = Initialization.O1;
            int bl = subset.Bl;

            // Reference parameters
            dpro.LnPhiF = 0.0;
            dpro.LnDRef = DotProduct(subset.Beta, 0, Constants.Mbf, LogReferenceDensity, gf, 0, Constants.Mbf);
            dpro.ZRef = Constants.ZetaRefO1;
            dpro.ZMin = Constants.NodesO1[3];
            dpro.ZHyd = Constants.ZetaRefO1;

            // Turbopause and scale heights
            dpro.ZetaM = subset.Beta[0, TurbopauseHeight - bl];
            dpro.HML = subset.Beta[0, LowerScaleHeight - bl];
            dpro.HMU = subset.Beta[0, UpperScaleHeight - bl];

            // Chapman term
            dpro.C = DotProduct(subset.Beta, 0, Constants.Mbf, ChapmanCoefficient, gf, 0, Constants.Mbf);
            dpro.ZetaC = subset.Beta[0, ChapmanReferenceHeight - bl];
            dpro.HC = subset.Beta[0, ChapmanScaleHeight - bl];

            // Chemical/dynamical term (includes solar flux, geomagnetic, and UT effects)
            dpro.R = DotProduct(subset.Beta, 0, Constants.Mbf, ChemicalDynamicCoefficient, gf, 0, Constants.Mbf);
            dpro.R += BasisFunctions.SFluxMod(ChemicalDynamicCoefficient, gf, subset, 0.0);
            dpro.R += BasisFunctions.GeoMag(
                ExtractGeomagneticParameters(subset, ChemicalDynamicCoefficient),
                geomagPrimary, geomagSecondary);
            dpro.R += BasisFunctions.UtDep(
                ExtractUTParameters(subset, ChemicalDynamicCoefficient),
                utDependent);

            dpro.ZetaR = subset.Beta[0, ChemicalDynamicReferenceHeight - bl];
            dpro.HR = subset.Beta[0, ChemicalDynamicScaleHeight - bl];

            // Unconstrained splines for O
            for (int splineIndex = 0; splineIndex < Constants.NsplO1; splineIndex++)
            {
                int columnIndex = splineIndex + 10;  // Spline coefficients start at column 10
                dpro.Cf[splineIndex] = DotProduct(subset.Beta, 0, Constants.Mbf, columnIndex, gf, 0, Constants.Mbf);
            }
        }

        /// <summary>
        /// Computes density profile parameters for Helium (He)
        /// </summary>
        private static void ComputeHeParameters(DensityParameters dpro, double[] gf,
                                                TemperatureProfile tpro,
                                                double[] geomagPrimary, double[,] geomagSecondary,
                                                double[] utDependent)
        {
            var subset = Initialization.HE;
            int bl = subset.Bl;

            // Reference parameters
            dpro.LnPhiF = Constants.LnVmr[5];
            dpro.LnDRef = tpro.LnDTotF + dpro.LnPhiF;
            dpro.ZRef = Constants.ZetaF;
            dpro.ZMin = -1.0;
            dpro.ZHyd = Constants.ZetaF;

            // Turbopause and scale heights
            dpro.ZetaM = subset.Beta[0, TurbopauseHeight - bl];
            dpro.HML = subset.Beta[0, LowerScaleHeight - bl];
            dpro.HMU = subset.Beta[0, UpperScaleHeight - bl];

            // Chemical/dynamical term (includes solar flux, geomagnetic, and UT effects)
            dpro.R = DotProduct(subset.Beta, 0, Constants.Mbf, ChemicalDynamicCoefficient, gf, 0, Constants.Mbf);
            dpro.R += BasisFunctions.SFluxMod(ChemicalDynamicCoefficient, gf, subset, 1.0);
            dpro.R += BasisFunctions.GeoMag(
                ExtractGeomagneticParameters(subset, ChemicalDynamicCoefficient),
                geomagPrimary, geomagSecondary);
            dpro.R += BasisFunctions.UtDep(
                ExtractUTParameters(subset, ChemicalDynamicCoefficient),
                utDependent);

            dpro.ZetaR = subset.Beta[0, ChemicalDynamicReferenceHeight - bl];
            dpro.HR = subset.Beta[0, ChemicalDynamicScaleHeight - bl];
        }

        /// <summary>
        /// Computes density profile parameters for Atomic Hydrogen (H)
        /// </summary>
        private static void ComputeH1Parameters(DensityParameters dpro, double[] gf,
                                                TemperatureProfile tpro,
                                                double[] geomagPrimary, double[,] geomagSecondary,
                                                double[] utDependent)
        {
            var subset = Initialization.H1;
            int bl = subset.Bl;

            // Reference parameters
            dpro.LnPhiF = 0.0;
            dpro.LnDRef = DotProduct(subset.Beta, 0, Constants.Mbf, LogReferenceDensity, gf, 0, Constants.Mbf);
            dpro.ZRef = Constants.ZetaA;
            dpro.ZMin = 75.0;
            dpro.ZHyd = Constants.ZetaF;

            // Turbopause and scale heights
            dpro.ZetaM = subset.Beta[0, TurbopauseHeight - bl];
            dpro.HML = subset.Beta[0, LowerScaleHeight - bl];
            dpro.HMU = subset.Beta[0, UpperScaleHeight - bl];

            // Chapman term
            dpro.C = DotProduct(subset.Beta, 0, Constants.Mbf, ChapmanCoefficient, gf, 0, Constants.Mbf);
            dpro.ZetaC = DotProduct(subset.Beta, 0, Constants.Mbf, ChapmanReferenceHeight, gf, 0, Constants.Mbf);
            dpro.HC = subset.Beta[0, ChapmanScaleHeight - bl];

            // Chemical/dynamical term (includes solar flux, geomagnetic, and UT effects)
            dpro.R = DotProduct(subset.Beta, 0, Constants.Mbf, ChemicalDynamicCoefficient, gf, 0, Constants.Mbf);
            dpro.R += BasisFunctions.SFluxMod(ChemicalDynamicCoefficient, gf, subset, 0.0);
            dpro.R += BasisFunctions.GeoMag(
                ExtractGeomagneticParameters(subset, ChemicalDynamicCoefficient),
                geomagPrimary, geomagSecondary);
            dpro.R += BasisFunctions.UtDep(
                ExtractUTParameters(subset, ChemicalDynamicCoefficient),
                utDependent);

            dpro.ZetaR = subset.Beta[0, ChemicalDynamicReferenceHeight - bl];
            dpro.HR = subset.Beta[0, ChemicalDynamicScaleHeight - bl];
        }

        /// <summary>
        /// Computes density profile parameters for Argon (Ar)
        /// </summary>
        private static void ComputeArParameters(DensityParameters dpro, double[] gf,
                                                TemperatureProfile tpro,
                                                double[] geomagPrimary, double[,] geomagSecondary,
                                                double[] utDependent)
        {
            var subset = Initialization.AR;
            int bl = subset.Bl;

            // Reference parameters
            dpro.LnPhiF = Constants.LnVmr[7];
            dpro.LnDRef = tpro.LnDTotF + dpro.LnPhiF;
            dpro.ZRef = Constants.ZetaF;
            dpro.ZMin = -1.0;
            dpro.ZHyd = Constants.ZetaF;

            // Turbopause and scale heights
            dpro.ZetaM = subset.Beta[0, TurbopauseHeight - bl];
            dpro.HML = subset.Beta[0, LowerScaleHeight - bl];
            dpro.HMU = subset.Beta[0, UpperScaleHeight - bl];

            // Chemical/dynamical term (includes geomagnetic and UT effects)
            dpro.R = DotProduct(subset.Beta, 0, Constants.Mbf, ChemicalDynamicCoefficient, gf, 0, Constants.Mbf);
            dpro.R += BasisFunctions.GeoMag(
                ExtractGeomagneticParameters(subset, ChemicalDynamicCoefficient),
                geomagPrimary, geomagSecondary);
            dpro.R += BasisFunctions.UtDep(
                ExtractUTParameters(subset, ChemicalDynamicCoefficient),
                utDependent);

            dpro.ZetaR = subset.Beta[0, ChemicalDynamicReferenceHeight - bl];
            dpro.HR = subset.Beta[0, ChemicalDynamicScaleHeight - bl];
        }

        /// <summary>
        /// Computes density profile parameters for Atomic Nitrogen (N)
        /// </summary>
        private static void ComputeN1Parameters(DensityParameters dpro, double[] gf,
                                                TemperatureProfile tpro,
                                                double[] geomagPrimary, double[,] geomagSecondary,
                                                double[] utDependent)
        {
            var subset = Initialization.N1;
            int bl = subset.Bl;

            // Reference parameters (includes solar flux, geomagnetic, and UT effects)
            dpro.LnPhiF = 0.0;
            dpro.LnDRef = DotProduct(subset.Beta, 0, Constants.Mbf, LogReferenceDensity, gf, 0, Constants.Mbf);
            dpro.LnDRef += BasisFunctions.SFluxMod(LogReferenceDensity, gf, subset, 0.0);
            dpro.LnDRef += BasisFunctions.GeoMag(
                ExtractGeomagneticParameters(subset, LogReferenceDensity),
                geomagPrimary, geomagSecondary);
            dpro.LnDRef += BasisFunctions.UtDep(
                ExtractUTParameters(subset, LogReferenceDensity),
                utDependent);

            dpro.ZRef = Constants.ZetaB;  // Different from other species!
            dpro.ZMin = 90.0;              // Different from other species!
            dpro.ZHyd = Constants.ZetaF;

            // Turbopause and scale heights
            dpro.ZetaM = subset.Beta[0, TurbopauseHeight - bl];
            dpro.HML = subset.Beta[0, LowerScaleHeight - bl];
            dpro.HMU = subset.Beta[0, UpperScaleHeight - bl];

            // Chapman term (from Beta, not DotProduct!)
            dpro.C = subset.Beta[0, ChapmanCoefficient - bl];
            dpro.ZetaC = subset.Beta[0, ChapmanReferenceHeight - bl];
            dpro.HC = subset.Beta[0, ChapmanScaleHeight - bl];

            // Chemical/dynamical term
            dpro.R = DotProduct(subset.Beta, 0, Constants.Mbf, ChemicalDynamicCoefficient, gf, 0, Constants.Mbf);
            dpro.ZetaR = subset.Beta[0, ChemicalDynamicReferenceHeight - bl];
            dpro.HR = subset.Beta[0, ChemicalDynamicScaleHeight - bl];
        }

        /// <summary>
        /// Computes density profile parameters for Anomalous Oxygen (O*)
        /// Legacy MSISE-00 formulation - simpler than other species
        /// </summary>
        private static void ComputeOAParameters(DensityParameters dpro, double[] gf,
                                                TemperatureProfile tpro,
                                                double[] geomagPrimary, double[,] geomagSecondary,
                                                double[] utDependent)
        {
            var subset = Initialization.OA;
            int bl = subset.Bl;

            // Reference parameters (includes geomagnetic effects only)
            dpro.LnDRef = DotProduct(subset.Beta, 0, Constants.Mbf, LogReferenceDensity, gf, 0, Constants.Mbf);
            dpro.LnDRef += BasisFunctions.GeoMag(
                ExtractGeomagneticParameters(subset, LogReferenceDensity),
                geomagPrimary, geomagSecondary);

            dpro.ZRef = Constants.ZetaRefOA;  // Different from other species!
            dpro.ZMin = 120.0;                 // Different from other species!
            dpro.ZHyd = 0.0;                   // Different from other species!

            // Chapman term (from Beta, not DotProduct!)
            dpro.C = subset.Beta[0, ChapmanCoefficient - bl];
            dpro.ZetaC = subset.Beta[0, ChapmanReferenceHeight - bl];
            dpro.HC = subset.Beta[0, ChapmanScaleHeight - bl];

            // No turbopause or R term - returns here in original
            // Legacy formulation uses simpler exponential profile
        }

        /// <summary>
        /// Computes density profile parameters for Nitric Oxide (NO)
        /// </summary>
        private static void ComputeNOParameters(DensityParameters dpro, double[] gf,
                                                TemperatureProfile tpro,
                                                double[] geomagPrimary, double[,] geomagSecondary,
                                                double[] utDependent)
        {
            var subset = Initialization.NO;
            int bl = subset.Bl;

            // Skip if parameters are not defined
            if (subset.Beta[0, LogReferenceDensity - bl] == 0.0)
            {
                dpro.LnDRef = 0.0;
                return;
            }

            // Reference parameters (includes geomagnetic effects only - no SFlux or UT!)
            dpro.LnPhiF = 0.0;
            dpro.LnDRef = DotProduct(subset.Beta, 0, Constants.Mbf, LogReferenceDensity, gf, 0, Constants.Mbf);
            dpro.LnDRef += BasisFunctions.GeoMag(
                ExtractGeomagneticParameters(subset, LogReferenceDensity),
                geomagPrimary, geomagSecondary);

            dpro.ZRef = Constants.ZetaRefNO;  // Different from other species!
            dpro.ZMin = 72.5;                  // Different from other species!
            dpro.ZHyd = Constants.ZetaRefNO;   // Different from other species!

            // Turbopause and scale heights (from DotProduct, not Beta!)
            dpro.ZetaM = DotProduct(subset.Beta, 0, Constants.Mbf, TurbopauseHeight, gf, 0, Constants.Mbf);
            dpro.HML = DotProduct(subset.Beta, 0, Constants.Mbf, LowerScaleHeight, gf, 0, Constants.Mbf);
            dpro.HMU = DotProduct(subset.Beta, 0, Constants.Mbf, UpperScaleHeight, gf, 0, Constants.Mbf);

            // Chapman term (DotProduct AND geomagnetic!)
            dpro.C = DotProduct(subset.Beta, 0, Constants.Mbf, ChapmanCoefficient, gf, 0, Constants.Mbf);
            dpro.C += BasisFunctions.GeoMag(
                ExtractGeomagneticParameters(subset, ChapmanCoefficient),
                geomagPrimary, geomagSecondary);
            dpro.ZetaC = DotProduct(subset.Beta, 0, Constants.Mbf, ChapmanReferenceHeight, gf, 0, Constants.Mbf);
            dpro.HC = DotProduct(subset.Beta, 0, Constants.Mbf, ChapmanScaleHeight, gf, 0, Constants.Mbf);

            // Chemical/dynamical term (from DotProduct, not Beta!)
            dpro.R = DotProduct(subset.Beta, 0, Constants.Mbf, ChemicalDynamicCoefficient, gf, 0, Constants.Mbf);
            dpro.ZetaR = DotProduct(subset.Beta, 0, Constants.Mbf, ChemicalDynamicReferenceHeight, gf, 0, Constants.Mbf);
            dpro.HR = DotProduct(subset.Beta, 0, Constants.Mbf, ChemicalDynamicScaleHeight, gf, 0, Constants.Mbf);

            // Unconstrained splines (with geomagnetic effects)
            for (int splineIndex = 0; splineIndex < Constants.NsplNO; splineIndex++)
            {
                int columnIndex = splineIndex + 10;  // Spline coefficients start at column 10
                dpro.Cf[splineIndex] = DotProduct(subset.Beta, 0, Constants.Mbf, columnIndex, gf, 0, Constants.Mbf);
                dpro.Cf[splineIndex] += BasisFunctions.GeoMag(
                    ExtractGeomagneticParameters(subset, columnIndex),
                    geomagPrimary, geomagSecondary);
            }
        }

        // ==================================================================================================
        // HELPER METHODS
        // ==================================================================================================

        /// <summary>
        /// Extracts UT-dependent parameters for a given column from a basis subset
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
        /// Computes dot product of two arrays with specified ranges
        /// </summary>
        private static double DotProduct(double[,] beta, int betaRow1, int betaRow2, int betaCol,
                                        double[] gf, int gfStart, int gfEnd)
        {
            double sum = 0.0;
            for (int j = betaRow1; j <= betaRow2; j++)
            {
                if (j >= gfStart && j <= gfEnd)
                {
                    sum += beta[j, betaCol] * gf[j];
                }
            }
            return sum;
        }

        // ==================================================================================================
        // DFN_X: Compute number density at given altitude
        // ==================================================================================================

        /// <summary>
        /// Compute number density at specified altitude given profile parameters.
        /// Uses hydrostatic equilibrium with temperature-dependent scale heights.
        /// </summary>
        /// <param name="z">Altitude (km)</param>
        /// <param name="tnz">Temperature at altitude (K)</param>
        /// <param name="lndtotz">Natural log of total number density at altitude</param>
        /// <param name="Vz">Indefinite integral V(z) at altitude</param>
        /// <param name="Wz">Indefinite integral W(z) at altitude</param>
        /// <param name="HRfact">Scale height correction factor</param>
        /// <param name="tpro">Temperature profile parameters</param>
        /// <param name="dpro">Density profile parameters</param>
        /// <returns>Number density (m^-3) or DMissing if below minimum altitude</returns>
        public static double DfnX(double z, double tnz, double lndtotz,
                                 double Vz, double Wz, double HRfact,
                                 TemperatureProfile tpro, DensityParameters dpro)
        {
            // Below minimum height of profile
            if (z < dpro.ZMin)
            {
                return Constants.DMissing;
            }

            // Special case: Anomalous Oxygen (legacy MSISE-00 formulation)
            if (dpro.ISpec == 9)
            {
                return ComputeAnomalousOxygenDensity(z, dpro);
            }

            // Special case: Nitric Oxide - skip if parameters not defined
            if (dpro.ISpec == 10 && dpro.LnDRef == 0.0)
            {
                return Constants.DMissing;
            }

            // Compute Chapman and logistic corrections
            double ccor = ComputeDensityCorrection(z, HRfact, dpro);

            // Below height where hydrostatic terms are needed
            if (z < dpro.ZHyd)
            {
                return ComputeDensityBelowHydrostaticRegion(z, lndtotz, ccor, dpro);
            }

            // Compute density using hydrostatic equilibrium
            return ComputeDensityWithHydrostaticEquilibrium(z, tnz, Vz, Wz, ccor, dpro);
        }

        /// <summary>
        /// Computes density for anomalous oxygen using legacy MSISE-00 formulation
        /// </summary>
        private static double ComputeAnomalousOxygenDensity(double altitude, DensityParameters dpro)
        {
            double lnDensity = dpro.LnDRef
                             - (altitude - dpro.ZRef) / Constants.HOA
                             - dpro.C * Math.Exp(-(altitude - dpro.ZetaC) / dpro.HC);

            return Math.Exp(lnDensity);
        }

        /// <summary>
        /// Computes Chapman and logistic density corrections based on species type
        /// </summary>
        private static double ComputeDensityCorrection(double altitude, double scaleHeightFactor,
                                                       DensityParameters dpro)
        {
            switch (dpro.ISpec)
            {
                case 2:  // N2
                case 3:  // O2
                case 5:  // He
                case 7:  // Ar
                    // Major species: logistic correction only
                    return ComputeLogisticCorrection(altitude, scaleHeightFactor, dpro);

                case 4:  // O
                case 6:  // H
                case 8:  // N
                case 10: // NO
                    // Minor species: Chapman and logistic corrections
                    return ComputeChapmanCorrection(altitude, dpro)
                         + ComputeLogisticCorrection(altitude, scaleHeightFactor, dpro);

                default:
                    return 0.0;
            }
        }

        /// <summary>
        /// Computes Chapman profile correction term
        /// </summary>
        private static double ComputeChapmanCorrection(double altitude, DensityParameters dpro)
        {
            return -dpro.C * Math.Exp(-(altitude - dpro.ZetaC) / dpro.HC);
        }

        /// <summary>
        /// Computes logistic function correction term
        /// </summary>
        private static double ComputeLogisticCorrection(double altitude, double scaleHeightFactor,
                                                        DensityParameters dpro)
        {
            return dpro.R * (1.0 + Math.Tanh((altitude - dpro.ZetaR) / (scaleHeightFactor * dpro.HR)));
        }

        /// <summary>
        /// Computes density below the hydrostatic region using mixing ratios or splines
        /// </summary>
        private static double ComputeDensityBelowHydrostaticRegion(double altitude, double lnTotalDensity,
                                                                    double correction, DensityParameters dpro)
        {
            switch (dpro.ISpec)
            {
                case 2:  // N2
                case 3:  // O2
                case 5:  // He
                case 7:  // Ar
                    // Major species: use mixing ratio
                    return Math.Exp(lnTotalDensity + dpro.LnPhiF + correction);

                case 4:  // O
                    // Atomic oxygen: use spline interpolation
                    return EvaluateOxygenSpline(altitude, dpro);

                case 10: // NO
                    // Nitric oxide: use spline interpolation
                    return EvaluateNitricOxideSpline(altitude, dpro);

                default:
                    // H, N, O*: shouldn't reach here (they use hydrostatic above ZMin)
                    return Math.Exp(lnTotalDensity + correction);
            }
        }

        /// <summary>
        /// Evaluates cubic B-spline for atomic oxygen density
        /// </summary>
        private static double EvaluateOxygenSpline(double altitude, DensityParameters dpro)
        {
            const int cubicOrder = 4;
            Utilities.BSpline(altitude, Constants.NodesO1, Constants.NdO1, cubicOrder,
                            Initialization.EtaO1, out double[,] splineValues, out int intervalIndex);

            // Compute dot product: dpro.Cf[iz-3:iz] · Sz[-3:0, 4]
            // Sz[-3:0, 4] maps to Sz[2:5, 2] in C# (spline indices -3,-2,-1,0 → array indices 2,3,4,5)
            double lnDensity = 0.0;
            for (int j = 0; j <= 3; j++)
            {
                int coefficientIndex = intervalIndex - 3 + j;
                int splineArrayIndex = j + 2;  // Maps Fortran -3+j to C# array index
                lnDensity += dpro.Cf[coefficientIndex] * splineValues[splineArrayIndex, 2];  // Order 4 → index 2
            }

            return Math.Exp(lnDensity);
        }

        /// <summary>
        /// Evaluates cubic B-spline for nitric oxide density
        /// </summary>
        private static double EvaluateNitricOxideSpline(double altitude, DensityParameters dpro)
        {
            const int cubicOrder = 4;
            Utilities.BSpline(altitude, Constants.NodesNO, Constants.NdNO, cubicOrder,
                            Initialization.EtaNO, out double[,] splineValues, out int intervalIndex);

            // Compute dot product: dpro.Cf[iz-3:iz] · Sz[-3:0, 4]
            double lnDensity = 0.0;
            for (int j = 0; j <= 3; j++)
            {
                int coefficientIndex = intervalIndex - 3 + j;
                int splineArrayIndex = j + 2;
                lnDensity += dpro.Cf[coefficientIndex] * splineValues[splineArrayIndex, 2];
            }

            return Math.Exp(lnDensity);
        }

        /// <summary>
        /// Computes density using full hydrostatic equilibrium
        /// </summary>
        private static double ComputeDensityWithHydrostaticEquilibrium(double altitude,
                                                                       double temperatureAtAltitude,
                                                                       double Vz, double Wz,
                                                                       double correction,
                                                                       DensityParameters dpro)
        {
            // Calculate effective mass at altitude
            double effectiveMass = ComputePiecewiseMass(altitude, dpro.ZetaMi, dpro.Mi, dpro.AMi);

            // Calculate hydrostatic integral
            double hydrostaticIntegral = effectiveMass * Vz - dpro.IzRef;

            // Apply piecewise corrections if within node range
            if (altitude > dpro.ZetaMi[0] && altitude < dpro.ZetaMi[4])
            {
                int nodeIndex = FindNodeIndex(altitude, dpro.ZetaMi);
                hydrostaticIntegral -= (dpro.AMi[nodeIndex] * Wz + dpro.XMi[nodeIndex]);
            }
            else if (altitude >= dpro.ZetaMi[4])
            {
                hydrostaticIntegral -= dpro.XMi[4];
            }

            // Compute log density at altitude
            double lnDensity = dpro.LnDRef - hydrostaticIntegral * Constants.G0DivKB + correction;

            // Apply ideal gas law: n(z) = n_ref * (T_ref / T(z)) * exp(...)
            double density = Math.Exp(lnDensity) * dpro.TRef / temperatureAtAltitude;

            return density;
        }

        /// <summary>
        /// Finds the node index for piecewise mass profile
        /// </summary>
        private static int FindNodeIndex(double altitude, double[] nodeHeights)
        {
            for (int nodeIndex = 0; nodeIndex <= 3; nodeIndex++)
            {
                if (altitude < nodeHeights[nodeIndex + 1])
                {
                    return nodeIndex;
                }
            }
            return 3; // Should not reach here if altitude < nodeHeights[4]
        }

        // ==================================================================================================
        // PWMP: Piecewise effective mass profile interpolation
        // ==================================================================================================

        /// <summary>
        /// Computes piecewise linear effective mass profile.
        /// Mass varies linearly between nodes to account for changing atmospheric composition.
        /// </summary>
        /// <param name="altitude">Altitude (km)</param>
        /// <param name="nodeHeights">Heights of mass profile nodes (km)</param>
        /// <param name="nodeMasses">Effective masses at nodes (amu)</param>
        /// <param name="massSlopes">Slopes of mass profile between nodes (amu/km)</param>
        /// <returns>Effective mass at altitude (amu)</returns>
        private static double ComputePiecewiseMass(double altitude, double[] nodeHeights,
                                                   double[] nodeMasses, double[] massSlopes)
        {
            // Most probable case: above highest node
            if (altitude >= nodeHeights[4])
            {
                return nodeMasses[4];
            }

            // Second most probable case: below lowest node
            if (altitude <= nodeHeights[0])
            {
                return nodeMasses[0];
            }

            // Between nodes: linear interpolation
            for (int nodeIndex = 0; nodeIndex <= 3; nodeIndex++)
            {
                if (altitude < nodeHeights[nodeIndex + 1])
                {
                    return nodeMasses[nodeIndex] + massSlopes[nodeIndex] * (altitude - nodeHeights[nodeIndex]);
                }
            }

            // Should never reach here
            throw new InvalidOperationException("Error in piecewise mass profile computation");
        }

    }
}