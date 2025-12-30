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
//
// MSISCALC: Interface with re-ordered input arguments and output arrays.
//
//     PREREQUISITES:
//       Must first run MsisInit.MsisInitialize to load parameters and set switches. The 
//       MsisCalc method checks for initialization and does a default
//       initialization if necessary. This self-initialization will be removed
//       in future versions.
//
//     CALLING SEQUENCE:
//       MsisCalc.Calculate(day, utsec, z, lat, lon, sfluxavg, sflux, ap, out tn, out dn, out tex);
//  
//     INPUT VARIABLES:
//       day       Day of year (1.0 to 365.0 or 366.0)
//       utsec     Universal time (seconds)
//       z         Geodetic altitude (km) (default) or Geopotential height (km)
//       lat       Geodetic latitude (deg)
//       lon       Geodetic longitude (deg)
//       sfluxavg  81 day average, centered on input time, of F10.7 solar
//                 activity index
//       sflux     Daily F10.7 for previous day
//       ap        Geomagnetic activity index array [0:6] (7 elements):
//                   [0] Daily Ap
//                   [1] 3 hr ap index for current time
//                   [2] 3 hr ap index for 3 hrs before current time
//                   [3] 3 hr ap index for 6 hrs before current time
//                   [4] 3 hr ap index for 9 hrs before current time
//                   [5] Average of eight 3 hr ap indices from 12 to 33 hrs
//                       prior to current time
//                   [6] Average of eight 3 hr ap indices from 36 to 57 hrs
//                       prior to current time
//                 ap[1:6] are only used when switch_legacy[8] = -1.0 in MsisInitialize
//
//     NOTES ON INPUT VARIABLES: 
//       - The day-of-year dependence of the model only uses the day argument. If
//         a continuous day-of-year dependence is desired, this argument should
//         include the fractional day (e.g., day = <day of year> + utsec/86400.0)
//       - If lZAltType = true (default) in the MsisInitialize call, then z is
//         treated as geodetic altitude.
//         If lZAltType = false, then z is treated as geopotential height.
//       - F107 and F107A values are the 10.7 cm radio flux at the Sun-Earth
//         distance, not the radio flux at 1 AU. 
//
//     OUTPUT VARIABLES:
//       tn     Temperature at altitude (K)
//       dn[0]  Total mass density (kg/m3)
//       dn[1]  N2 number density (m-3)
//       dn[2]  O2 number density (m-3)
//       dn[3]  O number density (m-3)
//       dn[4]  He number density (m-3)
//       dn[5]  H number density (m-3)
//       dn[6]  Ar number density (m-3)
//       dn[7]  N number density (m-3)
//       dn[8]  Anomalous oxygen number density (m-3)
//       dn[9]  NO number density (m-3)
//       tex    Exospheric temperature (K) (optional output)
//
//     NOTES ON OUTPUT VARIABLES: 
//       - Missing density values are returned as 9.999e-38
//       - Species included in mass density calculation are set in MsisInitialize
//
// ===========================================================================

using System;
using System.Linq;

namespace NRLMSIS
{
    /// <summary>
    /// Main MSIS calculation entry point
    /// </summary>
    public static class MsisCalc
    {
        private static double lastDay = -9999.0;
        private static double lastUtsec = -9999.0;
        private static double lastLat = -9999.0;
        private static double lastLon = -9999.0;
        private static double lastZ = -9999.0;
        private static double lastSflux = -9999.0;
        private static double lastSfluxavg = -9999.0;
        private static double[] lastAp = Enumerable.Repeat(-9999.0, 7).ToArray();
        private static double[] gf = new double[MsisConstants.MaxNbf];
        private static double[,] Sz = new double[6, 5]; // [-5:0, 2:6] -> [0:5, 0:4]
        private static int iz;
        private static TnParm tpro = new TnParm();
        private static DnParm[] dpro = Enumerable.Range(0, MsisConstants.NSpec)
                                                   .Select(_ => new DnParm())
                                                   .ToArray();

        // ==================================================================================================
        // MSISCALC: The main MSIS subroutine entry point
        // ==================================================================================================
        /// <summary>
        /// Calculate MSIS atmospheric parameters
        /// </summary>
        /// <param name="day">Day of year (1.0 to 365.0 or 366.0)</param>
        /// <param name="utsec">Universal time (seconds)</param>
        /// <param name="z">Geodetic altitude (km) or geopotential height (km)</param>
        /// <param name="lat">Geodetic latitude (degrees)</param>
        /// <param name="lon">Geodetic longitude (degrees)</param>
        /// <param name="sfluxavg">81-day average F10.7</param>
        /// <param name="sflux">Daily F10.7</param>
        /// <param name="ap">Geomagnetic activity index array [0:6]</param>
        /// <param name="tn">Output: Temperature at altitude (K)</param>
        /// <param name="dn">Output: Density array [0:9] - see documentation for details</param>
        /// <param name="tex">Output: Exospheric temperature (K)</param>
        public static void Calculate(double day, double utsec, double z, double lat, double lon,
                                    double sfluxavg, double sflux, double[] ap,
                                    out double tn, out double[] dn, out double tex)
        {
            // Check if model has been initialized; if not, perform default initialization
            if (!MsisInit.InitFlag)
            {
                MsisInit.MsisInitialize();
            }

            // Initialize output arrays
            dn = new double[10];
            tex = 0.0;

            double zeta;
            double lndtotz = 0.0;
            double Vz = 0.0;
            double Wz = 0.0;
            double HRfact;
            double lnPz;
            double delz;
            int i, j, kmax;

            // Calculate geopotential height, if necessary
            if (MsisInit.ZAltFlag)
            {
                zeta = MsisUtils.Alt2Gph(lat, z);
            }
            else
            {
                zeta = z;
            }

            // If only altitude changes then update the local spline weights
            if (zeta < MsisConstants.ZetaB)
            {
                if (zeta != lastZ)
                {
                    if (zeta < MsisConstants.ZetaF)
                    {
                        kmax = 5;
                    }
                    else
                    {
                        kmax = 6;
                    }
                    MsisUtils.BSpline(zeta, MsisConstants.NodesTN, MsisConstants.Nd + 2, kmax, 
                                     MsisInit.EtaTN, out Sz, out iz);
                    lastZ = zeta;
                }
            }

            // If location, time, or solar/geomagnetic conditions change then recompute the profile parameters
            bool needsUpdate = (day != lastDay) || (utsec != lastUtsec) ||
                             (lat != lastLat) || (lon != lastLon) ||
                             (sflux != lastSflux) || (sfluxavg != lastSfluxavg) ||
                             !ap.SequenceEqual(lastAp);

            if (needsUpdate)
            {
                MsisGfn.Globe(day, utsec, lat, lon, sfluxavg, sflux, ap, out gf);
                MsisTfn.TfnParm(gf, out tpro);
                // ispec ranges 2 to 10 (N2 through NO)
                // SpecFlag[ispec-1] checks if species ispec should be calculated
                // SpecFlag[1] checks species 2 (N2), SpecFlag[9] checks species 10 (NO)
                // dpro[2] = N2, dpro[3] = O2, ..., dpro[10] = NO
                for (int ispec = 2; ispec <= 10; ispec++)
                {
                    if (MsisInit.SpecFlag[ispec - 1])
                    {
                        MsisDfn.DfnParm(ispec, gf, tpro, out dpro[ispec]);
                    }
                }
                lastDay = day;
                lastUtsec = utsec;
                lastLat = lat;
                lastLon = lon;
                lastSflux = sflux;
                lastSfluxavg = sfluxavg;
                lastAp = (double[])ap.Clone();
            }

            // Exospheric temperature
            tex = tpro.Tex;

            // Temperature at altitude
            // Fortran Sz(-3:0, 4) maps to C# Sz[2:5, 2]
            double[] wght = new double[4];
            for (int idx = 0; idx < 4; idx++)
            {
                wght[idx] = Sz[idx + 2, 2];
            }
            tn = MsisTfn.TfnX(zeta, iz, wght, tpro);

            // Temperature integration terms at altitude, total number density
            delz = zeta - MsisConstants.ZetaB;
            if (zeta < MsisConstants.ZetaF)
            {
                // Below zetaF (70 km)
                i = Math.Max(iz - 4, 0);
                if (iz < 4)
                {
                    j = -iz;
                }
                else
                {
                    j = -4;
                }
                // Fortran: dot_product(tpro%beta(i:iz), Sz(j:0,5))
                // Sz(j:0,5) where j ranges from -4 to 0, order 5
                // In C#: Sz[j+5 to 5, 3] (since order 5 -> index 3, and j+5 maps -4->1, 0->5)
                Vz = 0.0;
                int szIdx = j + 5;
                for (int idx = i; idx <= iz; idx++)
                {
                    Vz += tpro.Beta[idx] * Sz[szIdx, 3];
                    szIdx++;
                }
                Vz += tpro.CVs;
                Wz = 0.0;
                lnPz = MsisConstants.LnP0 - MsisConstants.MbarG0DivKB * (Vz - tpro.VZeta0);
                lndtotz = lnPz - Math.Log(MsisConstants.KB * tn);
            }
            else
            {
                // Above zetaF (70 km)
                if (zeta < MsisConstants.ZetaB)
                {
                    // Fortran: dot_product(tpro%beta(iz-4:iz), Sz(-4:0,5))
                    // Sz(-4:0,5) maps to Sz[1:5, 3]
                    Vz = 0.0;
                    for (int idx = 0; idx <= 4; idx++)
                    {
                        Vz += tpro.Beta[iz - 4 + idx] * Sz[idx + 1, 3];
                    }
                    Vz += tpro.CVs;

                    // Fortran: dot_product(tpro%gamma(iz-5:iz), Sz(-5:0,6))
                    // Sz(-5:0,6) maps to Sz[0:5, 4]
                    Wz = 0.0;
                    for (int idx = 0; idx <= 5; idx++)
                    {
                        Wz += tpro.Gamma[iz - 5 + idx] * Sz[idx, 4];
                    }
                    Wz += tpro.CVs * delz + tpro.CWs;
                }
                else
                {
                    // Above zetaB (>122.5 km) - Bates profile region
                    Vz = (delz + Math.Log(tn / tpro.Tex) / tpro.Sigma) / tpro.Tex + tpro.CVb;
                    Wz = (0.5 * delz * delz + MsisUtils.Dilog(tpro.B * Math.Exp(-tpro.Sigma * delz)) / tpro.SigmaSq) / tpro.Tex
                       + tpro.CVb * delz + tpro.CWb;
                }
            }

            // Species number densities at altitude
            // Fortran: ispec ranges 2 to nspec-1 (which is 2 to 10)
            // C#: ispec ranges 2 to 10
            //     dn[1] = N2 (ispec=2), dn[2] = O2 (ispec=3), ..., dn[9] = NO (ispec=10)
            //     dpro[2] = N2, dpro[3] = O2, ..., dpro[10] = NO
            //     SpecFlag[1] = N2, SpecFlag[2] = O2, ..., SpecFlag[9] = NO
            HRfact = 0.5 * (1.0 + Math.Tanh(MsisConstants.HGamma * (zeta - MsisConstants.ZetaGamma)));
            for (int ispec = 2; ispec <= 10; ispec++)
            {
                if (MsisInit.SpecFlag[ispec - 1])
                {
                    dn[ispec - 1] = MsisDfn.DfnX(zeta, tn, lndtotz, Vz, Wz, HRfact, tpro, dpro[ispec]);
                }
                else
                {
                    dn[ispec - 1] = MsisConstants.DMissing;
                }
            }

            // Mass density
            // Fortran: dn(1) = dot_product(dn, masswgt)
            // C#: dn[0] corresponds to Fortran's dn(1)
            if (MsisInit.SpecFlag[0])
            {
                dn[0] = 0.0;
                for (int idx = 1; idx < 10; idx++)
                {
                    dn[0] += dn[idx] * MsisInit.MassWgt[idx];
                }
            }
            else
            {
                dn[0] = MsisConstants.DMissing;
            }
        }

        /// <summary>
        /// Calculate MSIS atmospheric parameters (overload without tex output)
        /// </summary>
        public static void Calculate(double day, double utsec, double z, double lat, double lon,
                                    double sfluxavg, double sflux, double[] ap,
                                    out double tn, out double[] dn)
        {
            Calculate(day, utsec, z, lat, lon, sfluxavg, sflux, ap, out tn, out dn, out _);
        }
    }
}