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
// - DnParm structure contains all density profile parameters for each species
// - Array indexing: Fortran arrays with negative indices mapped to C# 0-based arrays
// - Species indices: 2=N2, 3=O2, 4=O, 5=He, 6=H, 7=Ar, 8=N, 9=Anomalous O, 10=NO
// **************************************************************************************************

using System;

namespace NRLMSIS
{
    /// <summary>
    /// Density profile parameters structure
    /// </summary>
    public class DnParm
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
        public double[] Cf { get; set; } = new double[MsisConstants.NsplO1 + 2]; // Merged spline coefficients [0:nsplO1+1]
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
    public static class MsisDfn
    {
        // ==================================================================================================
        // DFNPARM: Compute the species density profile parameters
        // ==================================================================================================
        /// <summary>
        /// Compute the species density profile parameters
        /// </summary>
        /// <param name="ispec">Species index (2=N2, 3=O2, 4=O, 5=He, 6=H, 7=Ar, 8=N, 9=AnomalousO, 10=NO)</param>
        /// <param name="gf">Array of horizontal and temporal basis function terms [0:maxnbf-1]</param>
        /// <param name="tpro">Structure containing temperature vertical profile parameters</param>
        /// <param name="dpro">Output: density vertical profile parameters</param>
        public static void DfnParm(int ispec, double[] gf, TnParm tpro, out DnParm dpro)
        {
            dpro = new DnParm();
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
                double[] parms = new double[MsisConstants.NMag];
                for (int idx = 0; idx < MsisConstants.NMag; idx++)
                {
                    parms[idx] = subset.Beta[MsisConstants.CMag + idx, col - subset.Bl];
                }
                return parms;
            }

            double[] geomagBf1 = new double[13];
            Array.Copy(gf, MsisConstants.CMag, geomagBf1, 0, 13);
            double[,] geomagBf2 = new double[7, 2];
            int gfOffset = MsisConstants.CMag + 13;
            for (int n = 0; n <= 6; n++)
            {
                geomagBf2[n, 0] = gf[gfOffset + n];
                geomagBf2[n, 1] = gf[gfOffset + 7 + n];
            }

            double[] utdepParms;
            double[] utdepBf = new double[9];
            Array.Copy(gf, MsisConstants.CUt, utdepBf, 0, 9);

            switch (ispec)
            {
                case 2: // Molecular Nitrogen
                    dpro.LnPhiF = MsisConstants.LnVmr[ispec]; // ispec=2, access LnVmr[2]
                    dpro.LnDRef = tpro.LnDTotF + dpro.LnPhiF;
                    dpro.ZRef = MsisConstants.ZetaF;
                    dpro.ZMin = -1.0;
                    dpro.ZHyd = MsisConstants.ZetaF;
                    dpro.ZetaM = DotProduct(MsisInit.N2.Beta, 0, MsisConstants.Mbf, 1, gf, 0, MsisConstants.Mbf);
                    dpro.HML = MsisInit.N2.Beta[0, 2 - MsisInit.N2.Bl];
                    dpro.HMU = MsisInit.N2.Beta[0, 3 - MsisInit.N2.Bl];
                    dpro.R = 0.0;
                    if (MsisInit.N2RFlag)
                    {
                        dpro.R = DotProduct(MsisInit.N2.Beta, 0, MsisConstants.Mbf, 7, gf, 0, MsisConstants.Mbf);
                    }
                    dpro.ZetaR = MsisInit.N2.Beta[0, 8 - MsisInit.N2.Bl];
                    dpro.HR = MsisInit.N2.Beta[0, 9 - MsisInit.N2.Bl];
                    break;

                case 3: // Molecular Oxygen
                    dpro.LnPhiF = MsisConstants.LnVmr[ispec]; // ispec=3, access LnVmr[3]
                    dpro.LnDRef = tpro.LnDTotF + dpro.LnPhiF;
                    dpro.ZRef = MsisConstants.ZetaF;
                    dpro.ZMin = -1.0;
                    dpro.ZHyd = MsisConstants.ZetaF;
                    dpro.ZetaM = MsisInit.O2.Beta[0, 1 - MsisInit.O2.Bl];
                    dpro.HML = MsisInit.O2.Beta[0, 2 - MsisInit.O2.Bl];
                    dpro.HMU = MsisInit.O2.Beta[0, 3 - MsisInit.O2.Bl];
                    dpro.R = DotProduct(MsisInit.O2.Beta, 0, MsisConstants.Mbf, 7, gf, 0, MsisConstants.Mbf);
                    dpro.R += MsisGfn.GeoMag(ExtractGeomagParms(MsisInit.O2, 7), geomagBf1, geomagBf2);
                    dpro.ZetaR = MsisInit.O2.Beta[0, 8 - MsisInit.O2.Bl];
                    dpro.HR = MsisInit.O2.Beta[0, 9 - MsisInit.O2.Bl];
                    break;

                case 4: // Atomic Oxygen
                    dpro.LnPhiF = 0.0;
                    dpro.LnDRef = DotProduct(MsisInit.O1.Beta, 0, MsisConstants.Mbf, 0, gf, 0, MsisConstants.Mbf);
                    dpro.ZRef = MsisConstants.ZetaRefO1;
                    dpro.ZMin = MsisConstants.NodesO1[3];
                    dpro.ZHyd = MsisConstants.ZetaRefO1;
                    dpro.ZetaM = MsisInit.O1.Beta[0, 1 - MsisInit.O1.Bl];
                    dpro.HML = MsisInit.O1.Beta[0, 2 - MsisInit.O1.Bl];
                    dpro.HMU = MsisInit.O1.Beta[0, 3 - MsisInit.O1.Bl];
                    dpro.C = DotProduct(MsisInit.O1.Beta, 0, MsisConstants.Mbf, 4, gf, 0, MsisConstants.Mbf);
                    dpro.ZetaC = MsisInit.O1.Beta[0, 5 - MsisInit.O1.Bl];
                    dpro.HC = MsisInit.O1.Beta[0, 6 - MsisInit.O1.Bl];
                    dpro.R = DotProduct(MsisInit.O1.Beta, 0, MsisConstants.Mbf, 7, gf, 0, MsisConstants.Mbf);
                    dpro.R += MsisGfn.SFluxMod(7, gf, MsisInit.O1, 0.0);
                    dpro.R += MsisGfn.GeoMag(ExtractGeomagParms(MsisInit.O1, 7), geomagBf1, geomagBf2);
                    utdepParms = new double[MsisConstants.NUt];
                    for (int idx = 0; idx < MsisConstants.NUt; idx++)
                    {
                        utdepParms[idx] = MsisInit.O1.Beta[MsisConstants.CUt + idx, 7 - MsisInit.O1.Bl];
                    }
                    dpro.R += MsisGfn.UtDep(utdepParms, utdepBf);
                    dpro.ZetaR = MsisInit.O1.Beta[0, 8 - MsisInit.O1.Bl];
                    dpro.HR = MsisInit.O1.Beta[0, 9 - MsisInit.O1.Bl];
                    // Unconstrained splines
                    for (izf = 0; izf < MsisConstants.NsplO1; izf++)
                    {
                        dpro.Cf[izf] = DotProduct(MsisInit.O1.Beta, 0, MsisConstants.Mbf, izf + 10, gf, 0, MsisConstants.Mbf);
                    }
                    break;

                case 5: // Helium
                    dpro.LnPhiF = MsisConstants.LnVmr[ispec]; // ispec=5, access LnVmr[5]
                    dpro.LnDRef = tpro.LnDTotF + dpro.LnPhiF;
                    dpro.ZRef = MsisConstants.ZetaF;
                    dpro.ZMin = -1.0;
                    dpro.ZHyd = MsisConstants.ZetaF;
                    dpro.ZetaM = MsisInit.HE.Beta[0, 1 - MsisInit.HE.Bl];
                    dpro.HML = MsisInit.HE.Beta[0, 2 - MsisInit.HE.Bl];
                    dpro.HMU = MsisInit.HE.Beta[0, 3 - MsisInit.HE.Bl];
                    dpro.R = DotProduct(MsisInit.HE.Beta, 0, MsisConstants.Mbf, 7, gf, 0, MsisConstants.Mbf);
                    dpro.R += MsisGfn.SFluxMod(7, gf, MsisInit.HE, 1.0);
                    dpro.R += MsisGfn.GeoMag(ExtractGeomagParms(MsisInit.HE, 7), geomagBf1, geomagBf2);
                    utdepParms = new double[MsisConstants.NUt];
                    for (int idx = 0; idx < MsisConstants.NUt; idx++)
                    {
                        utdepParms[idx] = MsisInit.HE.Beta[MsisConstants.CUt + idx, 7 - MsisInit.HE.Bl];
                    }
                    dpro.R += MsisGfn.UtDep(utdepParms, utdepBf);
                    dpro.ZetaR = MsisInit.HE.Beta[0, 8 - MsisInit.HE.Bl];
                    dpro.HR = MsisInit.HE.Beta[0, 9 - MsisInit.HE.Bl];
                    break;

                case 6: // Atomic Hydrogen
                    dpro.LnPhiF = 0.0;
                    dpro.LnDRef = DotProduct(MsisInit.H1.Beta, 0, MsisConstants.Mbf, 0, gf, 0, MsisConstants.Mbf);
                    dpro.ZRef = MsisConstants.ZetaA;
                    dpro.ZMin = 75.0;
                    dpro.ZHyd = MsisConstants.ZetaF;
                    dpro.ZetaM = MsisInit.H1.Beta[0, 1 - MsisInit.H1.Bl];
                    dpro.HML = MsisInit.H1.Beta[0, 2 - MsisInit.H1.Bl];
                    dpro.HMU = MsisInit.H1.Beta[0, 3 - MsisInit.H1.Bl];
                    dpro.C = DotProduct(MsisInit.H1.Beta, 0, MsisConstants.Mbf, 4, gf, 0, MsisConstants.Mbf);
                    dpro.ZetaC = DotProduct(MsisInit.H1.Beta, 0, MsisConstants.Mbf, 5, gf, 0, MsisConstants.Mbf);
                    dpro.HC = MsisInit.H1.Beta[0, 6 - MsisInit.H1.Bl];
                    dpro.R = DotProduct(MsisInit.H1.Beta, 0, MsisConstants.Mbf, 7, gf, 0, MsisConstants.Mbf);
                    dpro.R += MsisGfn.SFluxMod(7, gf, MsisInit.H1, 0.0);
                    dpro.R += MsisGfn.GeoMag(ExtractGeomagParms(MsisInit.H1, 7), geomagBf1, geomagBf2);
                    utdepParms = new double[MsisConstants.NUt];
                    for (int idx = 0; idx < MsisConstants.NUt; idx++)
                    {
                        utdepParms[idx] = MsisInit.H1.Beta[MsisConstants.CUt + idx, 7 - MsisInit.H1.Bl];
                    }
                    dpro.R += MsisGfn.UtDep(utdepParms, utdepBf);
                    dpro.ZetaR = MsisInit.H1.Beta[0, 8 - MsisInit.H1.Bl];
                    dpro.HR = MsisInit.H1.Beta[0, 9 - MsisInit.H1.Bl];
                    break;

                case 7: // Argon
                    dpro.LnPhiF = MsisConstants.LnVmr[ispec]; // ispec=7, access LnVmr[7]
                    dpro.LnDRef = tpro.LnDTotF + dpro.LnPhiF;
                    dpro.ZRef = MsisConstants.ZetaF;
                    dpro.ZMin = -1.0;
                    dpro.ZHyd = MsisConstants.ZetaF;
                    dpro.ZetaM = MsisInit.AR.Beta[0, 1 - MsisInit.AR.Bl];
                    dpro.HML = MsisInit.AR.Beta[0, 2 - MsisInit.AR.Bl];
                    dpro.HMU = MsisInit.AR.Beta[0, 3 - MsisInit.AR.Bl];
                    dpro.R = DotProduct(MsisInit.AR.Beta, 0, MsisConstants.Mbf, 7, gf, 0, MsisConstants.Mbf);
                    dpro.R += MsisGfn.GeoMag(ExtractGeomagParms(MsisInit.AR, 7), geomagBf1, geomagBf2);
                    utdepParms = new double[MsisConstants.NUt];
                    for (int idx = 0; idx < MsisConstants.NUt; idx++)
                    {
                        utdepParms[idx] = MsisInit.AR.Beta[MsisConstants.CUt + idx, 7 - MsisInit.AR.Bl];
                    }
                    dpro.R += MsisGfn.UtDep(utdepParms, utdepBf);
                    dpro.ZetaR = MsisInit.AR.Beta[0, 8 - MsisInit.AR.Bl];
                    dpro.HR = MsisInit.AR.Beta[0, 9 - MsisInit.AR.Bl];
                    break;

                case 8: // Atomic Nitrogen
                    dpro.LnPhiF = 0.0;
                    dpro.LnDRef = DotProduct(MsisInit.N1.Beta, 0, MsisConstants.Mbf, 0, gf, 0, MsisConstants.Mbf);
                    dpro.LnDRef += MsisGfn.SFluxMod(0, gf, MsisInit.N1, 0.0);
                    dpro.LnDRef += MsisGfn.GeoMag(ExtractGeomagParms(MsisInit.N1, 0), geomagBf1, geomagBf2);
                    utdepParms = new double[MsisConstants.NUt];
                    for (int idx = 0; idx < MsisConstants.NUt; idx++)
                    {
                        utdepParms[idx] = MsisInit.N1.Beta[MsisConstants.CUt + idx, 0 - MsisInit.N1.Bl];
                    }
                    dpro.LnDRef += MsisGfn.UtDep(utdepParms, utdepBf);
                    dpro.ZRef = MsisConstants.ZetaB;
                    dpro.ZMin = 90.0;
                    dpro.ZHyd = MsisConstants.ZetaF;
                    dpro.ZetaM = MsisInit.N1.Beta[0, 1 - MsisInit.N1.Bl];
                    dpro.HML = MsisInit.N1.Beta[0, 2 - MsisInit.N1.Bl];
                    dpro.HMU = MsisInit.N1.Beta[0, 3 - MsisInit.N1.Bl];
                    dpro.C = MsisInit.N1.Beta[0, 4 - MsisInit.N1.Bl];
                    dpro.ZetaC = MsisInit.N1.Beta[0, 5 - MsisInit.N1.Bl];
                    dpro.HC = MsisInit.N1.Beta[0, 6 - MsisInit.N1.Bl];
                    dpro.R = DotProduct(MsisInit.N1.Beta, 0, MsisConstants.Mbf, 7, gf, 0, MsisConstants.Mbf);
                    dpro.ZetaR = MsisInit.N1.Beta[0, 8 - MsisInit.N1.Bl];
                    dpro.HR = MsisInit.N1.Beta[0, 9 - MsisInit.N1.Bl];
                    break;

                case 9: // Anomalous Oxygen
                    dpro.LnDRef = DotProduct(MsisInit.OA.Beta, 0, MsisConstants.Mbf, 0, gf, 0, MsisConstants.Mbf);
                    dpro.LnDRef += MsisGfn.GeoMag(ExtractGeomagParms(MsisInit.OA, 0), geomagBf1, geomagBf2);
                    dpro.ZRef = MsisConstants.ZetaRefOA;
                    dpro.ZMin = 120.0;
                    dpro.ZHyd = 0.0;
                    dpro.C = MsisInit.OA.Beta[0, 4 - MsisInit.OA.Bl];
                    dpro.ZetaC = MsisInit.OA.Beta[0, 5 - MsisInit.OA.Bl];
                    dpro.HC = MsisInit.OA.Beta[0, 6 - MsisInit.OA.Bl];
                    return; // No further parameters needed for legacy anomalous oxygen profile

                case 10: // Nitric Oxide
                    // Skip if parameters are not defined
                    if (MsisInit.NO.Beta[0, 0 - MsisInit.NO.Bl] == 0.0)
                    {
                        dpro.LnDRef = 0.0;
                        return;
                    }
                    dpro.LnPhiF = 0.0;
                    dpro.LnDRef = DotProduct(MsisInit.NO.Beta, 0, MsisConstants.Mbf, 0, gf, 0, MsisConstants.Mbf);
                    dpro.LnDRef += MsisGfn.GeoMag(ExtractGeomagParms(MsisInit.NO, 0), geomagBf1, geomagBf2);
                    dpro.ZRef = MsisConstants.ZetaRefNO;
                    dpro.ZMin = 72.5; // Cut off profile below 72.5 km
                    dpro.ZHyd = MsisConstants.ZetaRefNO;
                    dpro.ZetaM = DotProduct(MsisInit.NO.Beta, 0, MsisConstants.Mbf, 1, gf, 0, MsisConstants.Mbf);
                    dpro.HML = DotProduct(MsisInit.NO.Beta, 0, MsisConstants.Mbf, 2, gf, 0, MsisConstants.Mbf);
                    dpro.HMU = DotProduct(MsisInit.NO.Beta, 0, MsisConstants.Mbf, 3, gf, 0, MsisConstants.Mbf);
                    dpro.C = DotProduct(MsisInit.NO.Beta, 0, MsisConstants.Mbf, 4, gf, 0, MsisConstants.Mbf);
                    dpro.C += MsisGfn.GeoMag(ExtractGeomagParms(MsisInit.NO, 4), geomagBf1, geomagBf2);
                    dpro.ZetaC = DotProduct(MsisInit.NO.Beta, 0, MsisConstants.Mbf, 5, gf, 0, MsisConstants.Mbf);
                    dpro.HC = DotProduct(MsisInit.NO.Beta, 0, MsisConstants.Mbf, 6, gf, 0, MsisConstants.Mbf);
                    dpro.R = DotProduct(MsisInit.NO.Beta, 0, MsisConstants.Mbf, 7, gf, 0, MsisConstants.Mbf);
                    dpro.ZetaR = DotProduct(MsisInit.NO.Beta, 0, MsisConstants.Mbf, 8, gf, 0, MsisConstants.Mbf);
                    dpro.HR = DotProduct(MsisInit.NO.Beta, 0, MsisConstants.Mbf, 9, gf, 0, MsisConstants.Mbf);
                    // Unconstrained splines
                    for (izf = 0; izf < MsisConstants.NsplNO; izf++)
                    {
                        dpro.Cf[izf] = DotProduct(MsisInit.NO.Beta, 0, MsisConstants.Mbf, izf + 10, gf, 0, MsisConstants.Mbf);
                        dpro.Cf[izf] += MsisGfn.GeoMag(ExtractGeomagParms(MsisInit.NO, izf + 10), geomagBf1, geomagBf2);
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
            dpro.Mi[0] = MsisConstants.Mbar;
            dpro.Mi[4] = MsisConstants.SpecMass[ispec]; // Access SpecMass with ispec directly (2-10)
            dpro.Mi[2] = (dpro.Mi[0] + dpro.Mi[4]) / 2.0;
            double delM = MsisConstants.Tanh1 * (dpro.Mi[4] - dpro.Mi[0]) / 2.0;
            dpro.Mi[1] = dpro.Mi[2] - delM;
            dpro.Mi[3] = dpro.Mi[2] + delM;
            for (i = 0; i <= 3; i++)
            {
                dpro.AMi[i] = (dpro.Mi[i + 1] - dpro.Mi[i]) / (dpro.ZetaMi[i + 1] - dpro.ZetaMi[i]);
            }

            for (i = 0; i <= 4; i++)
            {
                double delz = dpro.ZetaMi[i] - MsisConstants.ZetaB;
                if (dpro.ZetaMi[i] < MsisConstants.ZetaB)
                {
                    MsisUtils.BSpline(dpro.ZetaMi[i], MsisConstants.NodesTN, MsisConstants.Nd + 2, 6, MsisInit.EtaTN, out Si, out iz);
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
                    dpro.WMi[i] = (0.5 * delz * delz + MsisUtils.Dilog(tpro.B * Math.Exp(-tpro.Sigma * delz)) / tpro.SigmaSq) / tpro.Tex
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
            if (dpro.ZRef == MsisConstants.ZetaF)
            {
                Mzref = MsisConstants.Mbar;
                dpro.TRef = tpro.TZetaF;
                dpro.IzRef = MsisConstants.Mbar * tpro.VZetaF;
            }
            else if (dpro.ZRef == MsisConstants.ZetaB)
            {
                Mzref = Pwmp(dpro.ZRef, dpro.ZetaMi, dpro.Mi, dpro.AMi);
                dpro.TRef = tpro.Tb0;
                dpro.IzRef = 0.0;
                if ((MsisConstants.ZetaB > dpro.ZetaMi[0]) && (MsisConstants.ZetaB < dpro.ZetaMi[4]))
                {
                    i = 0;
                    for (i1 = 1; i1 <= 3; i1++)
                    {
                        if (MsisConstants.ZetaB < dpro.ZetaMi[i1])
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
            else if (dpro.ZRef == MsisConstants.ZetaA)
            {
                Mzref = Pwmp(dpro.ZRef, dpro.ZetaMi, dpro.Mi, dpro.AMi);
                dpro.TRef = tpro.TZetaA;
                dpro.IzRef = Mzref * tpro.VZetaA;
                if ((MsisConstants.ZetaA > dpro.ZetaMi[0]) && (MsisConstants.ZetaA < dpro.ZetaMi[4]))
                {
                    i = 0;
                    for (i1 = 1; i1 <= 3; i1++)
                    {
                        if (MsisConstants.ZetaA < dpro.ZetaMi[i1])
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
                Rterm0 = Math.Tanh((dpro.ZRef - dpro.ZetaR) / (MsisInit.HRFactO1Ref * dpro.HR));
                Rterm = dpro.R * (1 + Rterm0);
                bc[0] = dpro.LnDRef - Cterm + Rterm - dpro.Cf[7] * MsisConstants.C1O1Adj[0];
                bc[1] = -Mzref * MsisConstants.G0DivKB / tpro.TZetaA
                        - tpro.DlnTdzA
                        + Cterm / dpro.HC
                        + Rterm * (1 - Rterm0) / dpro.HR * MsisInit.DHRFactO1Ref
                        - dpro.Cf[7] * MsisConstants.C1O1Adj[1];
                // Compute coefficients for constrained splines: bc * c1o1
                for (int idx = 8; idx <= 9; idx++)
                {
                    dpro.Cf[idx] = 0.0;
                    for (int row = 0; row < 2; row++)
                    {
                        dpro.Cf[idx] += bc[row] * MsisConstants.C1O1[row, idx - 8];
                    }
                }
            }

            // C1 constraint for NO at 122.5 km
            if (ispec == 10)
            {
                Cterm = dpro.C * Math.Exp(-(dpro.ZRef - dpro.ZetaC) / dpro.HC);
                Rterm0 = Math.Tanh((dpro.ZRef - dpro.ZetaR) / (MsisInit.HRFactNORef * dpro.HR));
                Rterm = dpro.R * (1 + Rterm0);
                bc[0] = dpro.LnDRef - Cterm + Rterm - dpro.Cf[7] * MsisConstants.C1NOAdj[0];
                bc[1] = -Mzref * MsisConstants.G0DivKB / tpro.Tb0
                        - tpro.Tgb0 / tpro.Tb0
                        + Cterm / dpro.HC
                        + Rterm * (1 - Rterm0) / dpro.HR * MsisInit.DHRFactNORef
                        - dpro.Cf[7] * MsisConstants.C1NOAdj[1];
                // Compute coefficients for constrained splines: bc * c1NO
                for (int idx = 8; idx <= 9; idx++)
                {
                    dpro.Cf[idx] = 0.0;
                    for (int row = 0; row < 2; row++)
                    {
                        dpro.Cf[idx] += bc[row] * MsisConstants.C1NO[row, idx - 8];
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
                                 double HRfact, TnParm tpro, DnParm dpro)
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
                return MsisConstants.DMissing;
            }

            // Anomalous Oxygen (legacy MSISE-00 formulation)
            if (dpro.ISpec == 9)
            {
                dfnx = dpro.LnDRef - (z - dpro.ZRef) / MsisConstants.HOA - dpro.C * Math.Exp(-(z - dpro.ZetaC) / dpro.HC);
                return Math.Exp(dfnx);
            }

            // Nitric Oxide: Skip if parameters are not defined
            if (dpro.ISpec == 10)
            {
                if (dpro.LnDRef == 0.0)
                {
                    return MsisConstants.DMissing;
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
                        MsisUtils.BSpline(z, MsisConstants.NodesO1, MsisConstants.NdO1, 4, MsisInit.EtaO1, out Sz, out iz);
                        // Fortran: dot_product(dpro%cf(iz-3:iz), Sz(-3:0,4))
                        // Sz(-3:0,4) maps to Sz[2:5, 2] (order 4 -> index 2)
                        dfnx = 0.0;
                        for (int j = 0; j <= 3; j++)
                        {
                            dfnx += dpro.Cf[iz - 3 + j] * Sz[j + 2, 2];
                        }
                        return Math.Exp(dfnx);

                    case 10: // For NO, evaluate splines
                        MsisUtils.BSpline(z, MsisConstants.NodesNO, MsisConstants.NdNO, 4, MsisInit.EtaNO, out Sz, out iz);
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

            dfnx = dpro.LnDRef - Ihyd * MsisConstants.G0DivKB + ccor;

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
            if (beta == MsisInit.TN.Beta) blOffset = MsisInit.TN.Bl;
            else if (beta == MsisInit.PR.Beta) blOffset = MsisInit.PR.Bl;
            else if (beta == MsisInit.N2.Beta) blOffset = MsisInit.N2.Bl;
            else if (beta == MsisInit.O2.Beta) blOffset = MsisInit.O2.Bl;
            else if (beta == MsisInit.O1.Beta) blOffset = MsisInit.O1.Bl;
            else if (beta == MsisInit.HE.Beta) blOffset = MsisInit.HE.Bl;
            else if (beta == MsisInit.H1.Beta) blOffset = MsisInit.H1.Bl;
            else if (beta == MsisInit.AR.Beta) blOffset = MsisInit.AR.Bl;
            else if (beta == MsisInit.N1.Beta) blOffset = MsisInit.N1.Bl;
            else if (beta == MsisInit.OA.Beta) blOffset = MsisInit.OA.Bl;
            else if (beta == MsisInit.NO.Beta) blOffset = MsisInit.NO.Bl;
            
            for (int i = betaRowStart; i <= betaRowEnd; i++)
            {
                sum += beta[i, betaCol - blOffset] * gf[gfIdx];
                gfIdx++;
            }
            return sum;
        }
    }
}