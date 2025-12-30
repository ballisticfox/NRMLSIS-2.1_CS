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

namespace NRLMSIS.Calculators
{
    /// <summary>
    /// Density profile parameters structure
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
    /// Density profile functions for each species
    /// </summary>
    public static class DensityProfile
    {
        // ==================================================================================================
        // DensityParameters: Compute the species density profile parameters
        // ==================================================================================================
        /// <summary>
        /// Compute the species density profile parameters
        /// </summary>
        /// <param name="ispec">Species index (2=N2, 3=O2, 4=O, 5=He, 6=H, 7=Ar, 8=N, 9=AnomalousO, 10=NO)</param>
        /// <param name="gf">Array of horizontal and temporal basis function terms [0:maxnbf-1]</param>
        /// <param name="tpro">Structure containing temperature vertical profile parameters</param>
        /// <param name="dpro">Output: density vertical profile parameters</param>
        public static void DensityParameters(int ispec, double[] gf, TemperatureProfile tpro, out DensityParameters dpro)
        {
            dpro = new DensityParameters();
            dpro.ISpec = ispec;

            int izf, i, i1, iz;
            double Cterm, Rterm0, Rterm;
            double[] bc = new double[2];
            double Wi;
            double[,] Si;
            double Mzref;

            // Helper to extract geomag parameters and basis functions
            double[] ExtractGeomagParms(BasisSubset subset, int col)
            {
                double[] parms = new double[Constants.NMag];
                for (int idx = 0; idx < Constants.NMag; idx++)
                {
                    parms[idx] = subset.Beta[Constants.CMag + idx, col - subset.Bl];
                }
                return parms;
            }

            double[] geomagBf1 = new double[13];
            Array.Copy(gf, Constants.CMag, geomagBf1, 0, 13);
            double[,] geomagBf2 = new double[7, 2];
            int gfOffset = Constants.CMag + 13;
            for (int n = 0; n <= 6; n++)
            {
                geomagBf2[n, 0] = gf[gfOffset + n];
                geomagBf2[n, 1] = gf[gfOffset + 7 + n];
            }

            double[] utdepParms;
            double[] utdepBf = new double[9];
            Array.Copy(gf, Constants.CUt, utdepBf, 0, 9);

            switch (ispec)
            {
                case 2: // Molecular Nitrogen
                    dpro.LnPhiF = Constants.LnVmr[ispec]; // ispec=2, access LnVmr[2]
                    dpro.LnDRef = tpro.LnDTotF + dpro.LnPhiF;
                    dpro.ZRef = Constants.ZetaF;
                    dpro.ZMin = -1.0;
                    dpro.ZHyd = Constants.ZetaF;
                    dpro.ZetaM = DotProduct(Initialization.N2.Beta, 0, Constants.Mbf, 1, gf, 0, Constants.Mbf);
                    dpro.HML = Initialization.N2.Beta[0, 2 - Initialization.N2.Bl];
                    dpro.HMU = Initialization.N2.Beta[0, 3 - Initialization.N2.Bl];
                    dpro.R = 0.0;
                    if (Initialization.N2RFlag)
                    {
                        dpro.R = DotProduct(Initialization.N2.Beta, 0, Constants.Mbf, 7, gf, 0, Constants.Mbf);
                    }
                    dpro.ZetaR = Initialization.N2.Beta[0, 8 - Initialization.N2.Bl];
                    dpro.HR = Initialization.N2.Beta[0, 9 - Initialization.N2.Bl];
                    break;

                case 3: // Molecular Oxygen
                    dpro.LnPhiF = Constants.LnVmr[ispec]; // ispec=3, access LnVmr[3]
                    dpro.LnDRef = tpro.LnDTotF + dpro.LnPhiF;
                    dpro.ZRef = Constants.ZetaF;
                    dpro.ZMin = -1.0;
                    dpro.ZHyd = Constants.ZetaF;
                    dpro.ZetaM = Initialization.O2.Beta[0, 1 - Initialization.O2.Bl];
                    dpro.HML = Initialization.O2.Beta[0, 2 - Initialization.O2.Bl];
                    dpro.HMU = Initialization.O2.Beta[0, 3 - Initialization.O2.Bl];
                    dpro.R = DotProduct(Initialization.O2.Beta, 0, Constants.Mbf, 7, gf, 0, Constants.Mbf);
                    dpro.R += BasisFunctions.GeoMag(ExtractGeomagParms(Initialization.O2, 7), geomagBf1, geomagBf2);
                    dpro.ZetaR = Initialization.O2.Beta[0, 8 - Initialization.O2.Bl];
                    dpro.HR = Initialization.O2.Beta[0, 9 - Initialization.O2.Bl];
                    break;

                case 4: // Atomic Oxygen
                    dpro.LnPhiF = 0.0;
                    dpro.LnDRef = DotProduct(Initialization.O1.Beta, 0, Constants.Mbf, 0, gf, 0, Constants.Mbf);
                    dpro.ZRef = Constants.ZetaRefO1;
                    dpro.ZMin = Constants.NodesO1[3];
                    dpro.ZHyd = Constants.ZetaRefO1;
                    dpro.ZetaM = Initialization.O1.Beta[0, 1 - Initialization.O1.Bl];
                    dpro.HML = Initialization.O1.Beta[0, 2 - Initialization.O1.Bl];
                    dpro.HMU = Initialization.O1.Beta[0, 3 - Initialization.O1.Bl];
                    dpro.C = DotProduct(Initialization.O1.Beta, 0, Constants.Mbf, 4, gf, 0, Constants.Mbf);
                    dpro.ZetaC = Initialization.O1.Beta[0, 5 - Initialization.O1.Bl];
                    dpro.HC = Initialization.O1.Beta[0, 6 - Initialization.O1.Bl];
                    dpro.R = DotProduct(Initialization.O1.Beta, 0, Constants.Mbf, 7, gf, 0, Constants.Mbf);
                    dpro.R += BasisFunctions.SFluxMod(7, gf, Initialization.O1, 0.0);
                    dpro.R += BasisFunctions.GeoMag(ExtractGeomagParms(Initialization.O1, 7), geomagBf1, geomagBf2);
                    utdepParms = new double[Constants.NUt];
                    for (int idx = 0; idx < Constants.NUt; idx++)
                    {
                        utdepParms[idx] = Initialization.O1.Beta[Constants.CUt + idx, 7 - Initialization.O1.Bl];
                    }
                    dpro.R += BasisFunctions.UtDep(utdepParms, utdepBf);
                    dpro.ZetaR = Initialization.O1.Beta[0, 8 - Initialization.O1.Bl];
                    dpro.HR = Initialization.O1.Beta[0, 9 - Initialization.O1.Bl];
                    // Unconstrained splines
                    for (izf = 0; izf < Constants.NsplO1; izf++)
                    {
                        dpro.Cf[izf] = DotProduct(Initialization.O1.Beta, 0, Constants.Mbf, izf + 10, gf, 0, Constants.Mbf);
                    }
                    break;

                case 5: // Helium
                    dpro.LnPhiF = Constants.LnVmr[ispec]; // ispec=5, access LnVmr[5]
                    dpro.LnDRef = tpro.LnDTotF + dpro.LnPhiF;
                    dpro.ZRef = Constants.ZetaF;
                    dpro.ZMin = -1.0;
                    dpro.ZHyd = Constants.ZetaF;
                    dpro.ZetaM = Initialization.HE.Beta[0, 1 - Initialization.HE.Bl];
                    dpro.HML = Initialization.HE.Beta[0, 2 - Initialization.HE.Bl];
                    dpro.HMU = Initialization.HE.Beta[0, 3 - Initialization.HE.Bl];
                    dpro.R = DotProduct(Initialization.HE.Beta, 0, Constants.Mbf, 7, gf, 0, Constants.Mbf);
                    dpro.R += BasisFunctions.SFluxMod(7, gf, Initialization.HE, 1.0);
                    dpro.R += BasisFunctions.GeoMag(ExtractGeomagParms(Initialization.HE, 7), geomagBf1, geomagBf2);
                    utdepParms = new double[Constants.NUt];
                    for (int idx = 0; idx < Constants.NUt; idx++)
                    {
                        utdepParms[idx] = Initialization.HE.Beta[Constants.CUt + idx, 7 - Initialization.HE.Bl];
                    }
                    dpro.R += BasisFunctions.UtDep(utdepParms, utdepBf);
                    dpro.ZetaR = Initialization.HE.Beta[0, 8 - Initialization.HE.Bl];
                    dpro.HR = Initialization.HE.Beta[0, 9 - Initialization.HE.Bl];
                    break;

                case 6: // Atomic Hydrogen
                    dpro.LnPhiF = 0.0;
                    dpro.LnDRef = DotProduct(Initialization.H1.Beta, 0, Constants.Mbf, 0, gf, 0, Constants.Mbf);
                    dpro.ZRef = Constants.ZetaA;
                    dpro.ZMin = 75.0;
                    dpro.ZHyd = Constants.ZetaF;
                    dpro.ZetaM = Initialization.H1.Beta[0, 1 - Initialization.H1.Bl];
                    dpro.HML = Initialization.H1.Beta[0, 2 - Initialization.H1.Bl];
                    dpro.HMU = Initialization.H1.Beta[0, 3 - Initialization.H1.Bl];
                    dpro.C = DotProduct(Initialization.H1.Beta, 0, Constants.Mbf, 4, gf, 0, Constants.Mbf);
                    dpro.ZetaC = DotProduct(Initialization.H1.Beta, 0, Constants.Mbf, 5, gf, 0, Constants.Mbf);
                    dpro.HC = Initialization.H1.Beta[0, 6 - Initialization.H1.Bl];
                    dpro.R = DotProduct(Initialization.H1.Beta, 0, Constants.Mbf, 7, gf, 0, Constants.Mbf);
                    dpro.R += BasisFunctions.SFluxMod(7, gf, Initialization.H1, 0.0);
                    dpro.R += BasisFunctions.GeoMag(ExtractGeomagParms(Initialization.H1, 7), geomagBf1, geomagBf2);
                    utdepParms = new double[Constants.NUt];
                    for (int idx = 0; idx < Constants.NUt; idx++)
                    {
                        utdepParms[idx] = Initialization.H1.Beta[Constants.CUt + idx, 7 - Initialization.H1.Bl];
                    }
                    dpro.R += BasisFunctions.UtDep(utdepParms, utdepBf);
                    dpro.ZetaR = Initialization.H1.Beta[0, 8 - Initialization.H1.Bl];
                    dpro.HR = Initialization.H1.Beta[0, 9 - Initialization.H1.Bl];
                    break;

                case 7: // Argon
                    dpro.LnPhiF = Constants.LnVmr[ispec]; // ispec=7, access LnVmr[7]
                    dpro.LnDRef = tpro.LnDTotF + dpro.LnPhiF;
                    dpro.ZRef = Constants.ZetaF;
                    dpro.ZMin = -1.0;
                    dpro.ZHyd = Constants.ZetaF;
                    dpro.ZetaM = Initialization.AR.Beta[0, 1 - Initialization.AR.Bl];
                    dpro.HML = Initialization.AR.Beta[0, 2 - Initialization.AR.Bl];
                    dpro.HMU = Initialization.AR.Beta[0, 3 - Initialization.AR.Bl];
                    dpro.R = DotProduct(Initialization.AR.Beta, 0, Constants.Mbf, 7, gf, 0, Constants.Mbf);
                    dpro.R += BasisFunctions.GeoMag(ExtractGeomagParms(Initialization.AR, 7), geomagBf1, geomagBf2);
                    utdepParms = new double[Constants.NUt];
                    for (int idx = 0; idx < Constants.NUt; idx++)
                    {
                        utdepParms[idx] = Initialization.AR.Beta[Constants.CUt + idx, 7 - Initialization.AR.Bl];
                    }
                    dpro.R += BasisFunctions.UtDep(utdepParms, utdepBf);
                    dpro.ZetaR = Initialization.AR.Beta[0, 8 - Initialization.AR.Bl];
                    dpro.HR = Initialization.AR.Beta[0, 9 - Initialization.AR.Bl];
                    break;

                case 8: // Atomic Nitrogen
                    dpro.LnPhiF = 0.0;
                    dpro.LnDRef = DotProduct(Initialization.N1.Beta, 0, Constants.Mbf, 0, gf, 0, Constants.Mbf);
                    dpro.LnDRef += BasisFunctions.SFluxMod(0, gf, Initialization.N1, 0.0);
                    dpro.LnDRef += BasisFunctions.GeoMag(ExtractGeomagParms(Initialization.N1, 0), geomagBf1, geomagBf2);
                    utdepParms = new double[Constants.NUt];
                    for (int idx = 0; idx < Constants.NUt; idx++)
                    {
                        utdepParms[idx] = Initialization.N1.Beta[Constants.CUt + idx, 0 - Initialization.N1.Bl];
                    }
                    dpro.LnDRef += BasisFunctions.UtDep(utdepParms, utdepBf);
                    dpro.ZRef = Constants.ZetaB;
                    dpro.ZMin = 90.0;
                    dpro.ZHyd = Constants.ZetaF;
                    dpro.ZetaM = Initialization.N1.Beta[0, 1 - Initialization.N1.Bl];
                    dpro.HML = Initialization.N1.Beta[0, 2 - Initialization.N1.Bl];
                    dpro.HMU = Initialization.N1.Beta[0, 3 - Initialization.N1.Bl];
                    dpro.C = Initialization.N1.Beta[0, 4 - Initialization.N1.Bl];
                    dpro.ZetaC = Initialization.N1.Beta[0, 5 - Initialization.N1.Bl];
                    dpro.HC = Initialization.N1.Beta[0, 6 - Initialization.N1.Bl];
                    dpro.R = DotProduct(Initialization.N1.Beta, 0, Constants.Mbf, 7, gf, 0, Constants.Mbf);
                    dpro.ZetaR = Initialization.N1.Beta[0, 8 - Initialization.N1.Bl];
                    dpro.HR = Initialization.N1.Beta[0, 9 - Initialization.N1.Bl];
                    break;

                case 9: // Anomalous Oxygen
                    dpro.LnDRef = DotProduct(Initialization.OA.Beta, 0, Constants.Mbf, 0, gf, 0, Constants.Mbf);
                    dpro.LnDRef += BasisFunctions.GeoMag(ExtractGeomagParms(Initialization.OA, 0), geomagBf1, geomagBf2);
                    dpro.ZRef = Constants.ZetaRefOA;
                    dpro.ZMin = 120.0;
                    dpro.ZHyd = 0.0;
                    dpro.C = Initialization.OA.Beta[0, 4 - Initialization.OA.Bl];
                    dpro.ZetaC = Initialization.OA.Beta[0, 5 - Initialization.OA.Bl];
                    dpro.HC = Initialization.OA.Beta[0, 6 - Initialization.OA.Bl];
                    return; // No further parameters needed for legacy anomalous oxygen profile

                case 10: // Nitric Oxide
                    // Skip if parameters are not defined
                    if (Initialization.NO.Beta[0, 0 - Initialization.NO.Bl] == 0.0)
                    {
                        dpro.LnDRef = 0.0;
                        return;
                    }
                    dpro.LnPhiF = 0.0;
                    dpro.LnDRef = DotProduct(Initialization.NO.Beta, 0, Constants.Mbf, 0, gf, 0, Constants.Mbf);
                    dpro.LnDRef += BasisFunctions.GeoMag(ExtractGeomagParms(Initialization.NO, 0), geomagBf1, geomagBf2);
                    dpro.ZRef = Constants.ZetaRefNO;
                    dpro.ZMin = 72.5; // Cut off profile below 72.5 km
                    dpro.ZHyd = Constants.ZetaRefNO;
                    dpro.ZetaM = DotProduct(Initialization.NO.Beta, 0, Constants.Mbf, 1, gf, 0, Constants.Mbf);
                    dpro.HML = DotProduct(Initialization.NO.Beta, 0, Constants.Mbf, 2, gf, 0, Constants.Mbf);
                    dpro.HMU = DotProduct(Initialization.NO.Beta, 0, Constants.Mbf, 3, gf, 0, Constants.Mbf);
                    dpro.C = DotProduct(Initialization.NO.Beta, 0, Constants.Mbf, 4, gf, 0, Constants.Mbf);
                    dpro.C += BasisFunctions.GeoMag(ExtractGeomagParms(Initialization.NO, 4), geomagBf1, geomagBf2);
                    dpro.ZetaC = DotProduct(Initialization.NO.Beta, 0, Constants.Mbf, 5, gf, 0, Constants.Mbf);
                    dpro.HC = DotProduct(Initialization.NO.Beta, 0, Constants.Mbf, 6, gf, 0, Constants.Mbf);
                    dpro.R = DotProduct(Initialization.NO.Beta, 0, Constants.Mbf, 7, gf, 0, Constants.Mbf);
                    dpro.ZetaR = DotProduct(Initialization.NO.Beta, 0, Constants.Mbf, 8, gf, 0, Constants.Mbf);
                    dpro.HR = DotProduct(Initialization.NO.Beta, 0, Constants.Mbf, 9, gf, 0, Constants.Mbf);
                    // Unconstrained splines
                    for (izf = 0; izf < Constants.NsplNO; izf++)
                    {
                        dpro.Cf[izf] = DotProduct(Initialization.NO.Beta, 0, Constants.Mbf, izf + 10, gf, 0, Constants.Mbf);
                        dpro.Cf[izf] += BasisFunctions.GeoMag(ExtractGeomagParms(Initialization.NO, izf + 10), geomagBf1, geomagBf2);
                    }
                    break;

                default:
                    throw new InvalidOperationException("Species not yet implemented");
            }

            // Compute piecewise mass profile values and integration terms
            dpro.ZetaMi[0] = dpro.ZetaM - 2.0 * dpro.HML;
            dpro.ZetaMi[1] = dpro.ZetaM - dpro.HML;
            dpro.ZetaMi[2] = dpro.ZetaM;
            dpro.ZetaMi[3] = dpro.ZetaM + dpro.HMU;
            dpro.ZetaMi[4] = dpro.ZetaM + 2.0 * dpro.HMU;
            dpro.Mi[0] = Constants.Mbar;
            dpro.Mi[4] = Constants.SpecMass[ispec]; // Access SpecMass with ispec directly (2-10)
            dpro.Mi[2] = (dpro.Mi[0] + dpro.Mi[4]) / 2.0;
            double delM = Constants.Tanh1 * (dpro.Mi[4] - dpro.Mi[0]) / 2.0;
            dpro.Mi[1] = dpro.Mi[2] - delM;
            dpro.Mi[3] = dpro.Mi[2] + delM;
            for (i = 0; i <= 3; i++)
            {
                dpro.AMi[i] = (dpro.Mi[i + 1] - dpro.Mi[i]) / (dpro.ZetaMi[i + 1] - dpro.ZetaMi[i]);
            }

            for (i = 0; i <= 4; i++)
            {
                double delz = dpro.ZetaMi[i] - Constants.ZetaB;
                if (dpro.ZetaMi[i] < Constants.ZetaB)
                {
                    Utilities.BSpline(dpro.ZetaMi[i], Constants.NodesTN, Constants.Nd + 2, 6, Initialization.EtaTN, out Si, out iz);
                    // Extract Si[-5:0, 6] -> Si[0:5, 4]
                    Wi = 0.0;
                    for (int j = 0; j <= 5; j++)
                    {
                        Wi += tpro.Gamma[iz - 5 + j] * Si[j, 4];
                    }
                    dpro.WMi[i] = Wi + tpro.CVs * delz + tpro.CWs;
                }
                else
                {
                    dpro.WMi[i] = (0.5 * delz * delz + Utilities.Dilog(tpro.B * Math.Exp(-tpro.Sigma * delz)) / tpro.SigmaSq) / tpro.Tex
                                + tpro.CVb * delz + tpro.CWb;
                }
            }

            dpro.XMi[0] = -dpro.AMi[0] * dpro.WMi[0];
            for (i = 1; i <= 3; i++)
            {
                dpro.XMi[i] = dpro.XMi[i - 1] - dpro.WMi[i] * (dpro.AMi[i] - dpro.AMi[i - 1]);
            }
            dpro.XMi[4] = dpro.XMi[3] + dpro.WMi[4] * dpro.AMi[3];

            // Calculate hydrostatic integral at reference height, and copy temperature
            if (dpro.ZRef == Constants.ZetaF)
            {
                Mzref = Constants.Mbar;
                dpro.TRef = tpro.TZetaF;
                dpro.IzRef = Constants.Mbar * tpro.VZetaF;
            }
            else if (dpro.ZRef == Constants.ZetaB)
            {
                Mzref = Pwmp(dpro.ZRef, dpro.ZetaMi, dpro.Mi, dpro.AMi);
                dpro.TRef = tpro.Tb0;
                dpro.IzRef = 0.0;
                if ((Constants.ZetaB > dpro.ZetaMi[0]) && (Constants.ZetaB < dpro.ZetaMi[4]))
                {
                    i = 0;
                    for (i1 = 1; i1 <= 3; i1++)
                    {
                        if (Constants.ZetaB < dpro.ZetaMi[i1])
                        {
                            break;
                        }
                        else
                        {
                            i = i1;
                        }
                    }
                    dpro.IzRef = dpro.IzRef - dpro.XMi[i];
                }
                else
                {
                    dpro.IzRef = dpro.IzRef - dpro.XMi[4];
                }
            }
            else if (dpro.ZRef == Constants.ZetaA)
            {
                Mzref = Pwmp(dpro.ZRef, dpro.ZetaMi, dpro.Mi, dpro.AMi);
                dpro.TRef = tpro.TZetaA;
                dpro.IzRef = Mzref * tpro.VZetaA;
                if ((Constants.ZetaA > dpro.ZetaMi[0]) && (Constants.ZetaA < dpro.ZetaMi[4]))
                {
                    i = 0;
                    for (i1 = 1; i1 <= 3; i1++)
                    {
                        if (Constants.ZetaA < dpro.ZetaMi[i1])
                        {
                            break;
                        }
                        else
                        {
                            i = i1;
                        }
                    }
                    dpro.IzRef = dpro.IzRef - (dpro.AMi[i] * tpro.WZetaA + dpro.XMi[i]);
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
            if (ispec == 4)
            {
                Cterm = dpro.C * Math.Exp(-(dpro.ZRef - dpro.ZetaC) / dpro.HC);
                Rterm0 = Math.Tanh((dpro.ZRef - dpro.ZetaR) / (Initialization.HRFactO1Ref * dpro.HR));
                Rterm = dpro.R * (1 + Rterm0);
                bc[0] = dpro.LnDRef - Cterm + Rterm - dpro.Cf[7] * Constants.C1O1Adj[0];
                bc[1] = -Mzref * Constants.G0DivKB / tpro.TZetaA
                        - tpro.DlnTdzA
                        + Cterm / dpro.HC
                        + Rterm * (1 - Rterm0) / dpro.HR * Initialization.DHRFactO1Ref
                        - dpro.Cf[7] * Constants.C1O1Adj[1];
                // Compute coefficients for constrained splines: bc * c1o1
                for (int idx = 8; idx <= 9; idx++)
                {
                    dpro.Cf[idx] = 0.0;
                    for (int row = 0; row < 2; row++)
                    {
                        dpro.Cf[idx] += bc[row] * Constants.C1O1[row, idx - 8];
                    }
                }
            }

            // C1 constraint for NO at 122.5 km
            if (ispec == 10)
            {
                Cterm = dpro.C * Math.Exp(-(dpro.ZRef - dpro.ZetaC) / dpro.HC);
                Rterm0 = Math.Tanh((dpro.ZRef - dpro.ZetaR) / (Initialization.HRFactNORef * dpro.HR));
                Rterm = dpro.R * (1 + Rterm0);
                bc[0] = dpro.LnDRef - Cterm + Rterm - dpro.Cf[7] * Constants.C1NOAdj[0];
                bc[1] = -Mzref * Constants.G0DivKB / tpro.Tb0
                        - tpro.Tgb0 / tpro.Tb0
                        + Cterm / dpro.HC
                        + Rterm * (1 - Rterm0) / dpro.HR * Initialization.DHRFactNORef
                        - dpro.Cf[7] * Constants.C1NOAdj[1];
                // Compute coefficients for constrained splines: bc * c1NO
                for (int idx = 8; idx <= 9; idx++)
                {
                    dpro.Cf[idx] = 0.0;
                    for (int row = 0; row < 2; row++)
                    {
                        dpro.Cf[idx] += bc[row] * Constants.C1NO[row, idx - 8];
                    }
                }
            }
        }

        // ==================================================================================================
        // DFNX: Compute a species density at specified geopotential height
        // ==================================================================================================
        /// <summary>
        /// Compute a species density at specified geopotential height
        /// </summary>
        /// <param name="z">Geopotential height (km)</param>
        /// <param name="tnz">Temperature at z</param>
        /// <param name="lndtotz">Total number density at z</param>
        /// <param name="Vz">First indefinite integral of 1/T at z</param>
        /// <param name="Wz">Second indefinite integral of 1/T at z</param>
        /// <param name="HRfact">Reduction factor for chemical/dynamical correction scale height</param>
        /// <param name="tpro">Temperature vertical profile parameters</param>
        /// <param name="dpro">Density vertical profile parameters</param>
        /// <returns>Species density (m^-3)</returns>
        public static double DfnX(double z, double tnz, double lndtotz, double Vz, double Wz,
                                 double HRfact, TemperatureProfile tpro, DensityParameters dpro)
        {
            int i, i1, iz;
            double Mz;
            double[,] Sz;
            double Ihyd;
            double ccor;
            double dfnx;

            // Below minimum height of profile
            if (z < dpro.ZMin)
            {
                return Constants.DMissing;
            }

            // Anomalous Oxygen (legacy MSISE-00 formulation)
            if (dpro.ISpec == 9)
            {
                dfnx = dpro.LnDRef - (z - dpro.ZRef) / Constants.HOA - dpro.C * Math.Exp(-(z - dpro.ZetaC) / dpro.HC);
                return Math.Exp(dfnx);
            }

            // Nitric Oxide: Skip if parameters are not defined
            if (dpro.ISpec == 10)
            {
                if (dpro.LnDRef == 0.0)
                {
                    return Constants.DMissing;
                }
            }

            // Chapman and logistic corrections
            switch (dpro.ISpec)
            {
                case 2:
                case 3:
                case 5:
                case 7: // For N2, O2, He, and Ar: logistic correction only
                    ccor = dpro.R * (1 + Math.Tanh((z - dpro.ZetaR) / (HRfact * dpro.HR)));
                    break;
                case 4:
                case 6:
                case 8:
                case 10: // For O, H, N, and NO: Chapman and logistic corrections
                    ccor = -dpro.C * Math.Exp(-(z - dpro.ZetaC) / dpro.HC)
                         + dpro.R * (1 + Math.Tanh((z - dpro.ZetaR) / (HRfact * dpro.HR)));
                    break;
                default:
                    ccor = 0.0;
                    break;
            }

            // Below height where hydrostatic terms are needed
            if (z < dpro.ZHyd)
            {
                switch (dpro.ISpec)
                {
                    case 2:
                    case 3:
                    case 5:
                    case 7: // For N2, O2, He, and Ar, apply mixing ratios and exit
                        return Math.Exp(lndtotz + dpro.LnPhiF + ccor);

                    case 4: // For O, evaluate splines
                        Utilities.BSpline(z, Constants.NodesO1, Constants.NdO1, 4, Initialization.EtaO1, out Sz, out iz);
                        // Fortran: dot_product(dpro%cf(iz-3:iz), Sz(-3:0,4))
                        // Sz(-3:0,4) maps to Sz[2:5, 2] (order 4 -> index 2)
                        dfnx = 0.0;
                        for (int j = 0; j <= 3; j++)
                        {
                            dfnx += dpro.Cf[iz - 3 + j] * Sz[j + 2, 2];
                        }
                        return Math.Exp(dfnx);

                    case 10: // For NO, evaluate splines
                        Utilities.BSpline(z, Constants.NodesNO, Constants.NdNO, 4, Initialization.EtaNO, out Sz, out iz);
                        // Fortran: dot_product(dpro%cf(iz-3:iz), Sz(-3:0,4))
                        // Sz(-3:0,4) maps to Sz[2:5, 2] (order 4 -> index 2)
                        dfnx = 0.0;
                        for (int j = 0; j <= 3; j++)
                        {
                            dfnx += dpro.Cf[iz - 3 + j] * Sz[j + 2, 2];
                        }
                        return Math.Exp(dfnx);
                }
            }

            // Calculate hydrostatic term and apply to reference density
            Mz = Pwmp(z, dpro.ZetaMi, dpro.Mi, dpro.AMi);
            Ihyd = Mz * Vz - dpro.IzRef;
            if ((z > dpro.ZetaMi[0]) && (z < dpro.ZetaMi[4]))
            {
                i = 0;
                for (i1 = 1; i1 <= 3; i1++)
                {
                    if (z < dpro.ZetaMi[i1])
                    {
                        break;
                    }
                    else
                    {
                        i = i1;
                    }
                }
                Ihyd = Ihyd - (dpro.AMi[i] * Wz + dpro.XMi[i]);
            }
            else if (z >= dpro.ZetaMi[4])
            {
                Ihyd = Ihyd - dpro.XMi[4];
            }

            dfnx = dpro.LnDRef - Ihyd * Constants.G0DivKB + ccor;

            // Apply ideal gas law
            dfnx = Math.Exp(dfnx) * dpro.TRef / tnz;

            return dfnx;
        }

        // ==================================================================================================
        // PWMP: Piecewise effective mass profile interpolation
        // ==================================================================================================
        /// <summary>
        /// Piecewise effective mass profile interpolation
        /// </summary>
        private static double Pwmp(double z, double[] zm, double[] m, double[] dmdz)
        {
            // Most probable case
            if (z >= zm[4])
            {
                return m[4];
            }

            // Second most probable case
            if (z <= zm[0])
            {
                return m[0];
            }

            // None of the above
            for (int inode = 0; inode <= 3; inode++)
            {
                if (z < zm[inode + 1])
                {
                    return m[inode] + dmdz[inode] * (z - zm[inode]);
                }
            }

            // If we are here this is a problem
            throw new InvalidOperationException("Error in pwmp");
        }

        // Helper method for dot product between subset Beta array and gf array
        private static double DotProduct(double[,] beta, int betaRowStart, int betaRowEnd, int betaCol,
                                        double[] gf, int gfStart, int gfEnd)
        {
            double sum = 0.0;
            int gfIdx = gfStart;

            // Determine which subset this is to get the correct Bl offset
            int blOffset = 0;
            if (beta == Initialization.TN.Beta) blOffset = Initialization.TN.Bl;
            else if (beta == Initialization.PR.Beta) blOffset = Initialization.PR.Bl;
            else if (beta == Initialization.N2.Beta) blOffset = Initialization.N2.Bl;
            else if (beta == Initialization.O2.Beta) blOffset = Initialization.O2.Bl;
            else if (beta == Initialization.O1.Beta) blOffset = Initialization.O1.Bl;
            else if (beta == Initialization.HE.Beta) blOffset = Initialization.HE.Bl;
            else if (beta == Initialization.H1.Beta) blOffset = Initialization.H1.Bl;
            else if (beta == Initialization.AR.Beta) blOffset = Initialization.AR.Bl;
            else if (beta == Initialization.N1.Beta) blOffset = Initialization.N1.Bl;
            else if (beta == Initialization.OA.Beta) blOffset = Initialization.OA.Bl;
            else if (beta == Initialization.NO.Beta) blOffset = Initialization.NO.Bl;

            for (int i = betaRowStart; i <= betaRowEnd; i++)
            {
                sum += beta[i, betaCol - blOffset] * gf[gfIdx];
                gfIdx++;
            }
            return sum;
        }
    }
}