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
// GTD8D: Legacy wrapper with input and output arguments used in NRLMSISE-00
//
//     PREREQUISITES:
//       Must first run Initialization.Initializationialize to load parameters and set switches. The
//       MsisCalc.Calculate method checks for initialization and does a default
//       initialization if necessary. This self-initialization will be removed
//       in future versions.
//
//     CALLING SEQUENCE:
//       Gtd8d.Calculate(iyd, sec, alt, glat, glong, stl, f107a, f107, ap, mass, out d, out t);
//
//     INPUT VARIABLES:
//       iyd    Year and day as YYDDD (day of year from 1 to 365 (or 366))
//                (Year is ignored in current model)
//       sec    Universal time (seconds)
//       alt    Geodetic altitude (km)
//       glat   Geodetic latitude (deg)
//       glong  Geodetic longitude (deg)
//       stl    Local solar time (Ignored; calculated from sec and glong)
//       f107a  81 day average, centered on input time, of F10.7 solar activity index
//       f107   Daily F10.7 for previous day
//       ap     Geomagnetic activity index array [0:6] (7 elements):
//                [0] Daily Ap
//                [1] 3 hr ap index for current time
//                [2] 3 hr ap index for 3 hrs before current time
//                [3] 3 hr ap index for 6 hrs before current time
//                [4] 3 hr ap index for 9 hrs before current time
//                [5] Average of eight 3 hr ap indices from 12 to 33 hrs
//                    prior to current time
//                [6] Average of eight 3 hr ap indices from 36 to 57 hrs
//                    prior to current time
//              ap[1:6] are only used when switch_legacy[8] = -1.0 in Initializationialize
//       mass   Mass number (Ignored in 2.1)
//
//     NOTES ON INPUT VARIABLES:
//       - If lZAltType = false in the Initializationialize call, then the alt input
//         argument is treated as geopotential height.
//       - The stl input argument is ignored in NRLMSIS 2.1. Instead, local time
//         is computed from universal time and longitude.
//       - F107 and F107A values are the 10.7 cm radio flux at the Sun-Earth
//         distance, not the radio flux at 1 AU.
//       - The mass input argument is ignored in NRLMSIS 2.1; species to be
//         calculated are set in Initializationialize.
//
//     OUTPUT VARIABLES:
//       d[0]  He number density (cm-3)
//       d[1]  O number density (cm-3)
//       d[2]  N2 number density (cm-3)
//       d[3]  O2 number density (cm-3)
//       d[4]  Ar number density (cm-3)
//       d[5]  Total mass density (g/cm3)
//       d[6]  H number density (cm-3)
//       d[7]  N number density (cm-3)
//       d[8]  Anomalous oxygen number density (cm-3)
//       d[9]  NO number density (cm-3)
//       t[0]  Exospheric temperature (K)
//       t[1]  Temperature at altitude (K)
//
//     NOTES ON OUTPUT VARIABLES:
//       - Missing density values are returned as 9.999e-38
//       - Species included in mass density calculation are set in Initializationialize
//
// ===========================================================================

using System;

namespace NRLMSIS
{
    /// <summary>
    /// Legacy NRLMSISE-00 compatible interface wrapper
    /// </summary>
    public static class LegacyInterface
    {
        // ==================================================================================================
        // Legacy wrapper
        // ==================================================================================================
        /// <summary>
        /// Legacy NRLMSISE-00 style interface for NRLMSIS 2.1
        /// </summary>
        /// <param name="iyd">Year and day as YYDDD (day of year from 1 to 365 or 366)</param>
        /// <param name="sec">Universal time (seconds)</param>
        /// <param name="alt">Geodetic altitude (km)</param>
        /// <param name="glat">Geodetic latitude (degrees)</param>
        /// <param name="glong">Geodetic longitude (degrees)</param>
        /// <param name="stl">Local solar time (ignored)</param>
        /// <param name="f107a">81-day average F10.7</param>
        /// <param name="f107">Daily F10.7</param>
        /// <param name="ap">Geomagnetic activity index array [0:6]</param>
        /// <param name="mass">Mass number (ignored)</param>
        /// <param name="d">Output: Density array [0:9] in legacy format (cm-3 and g/cm3)</param>
        /// <param name="t">Output: Temperature array [0:1] - [0]=exospheric, [1]=at altitude</param>
        public static void Calculate(int iyd, float sec, float alt, float glat, float glong, float stl,
                                    float f107a, float f107, float[] ap, int mass,
                                    out float[] d, out float[] t)
        {
            // Convert the legacy input arguments to the new interface values and precision
            double xday = iyd % 1000; // Extract day of year from YYDDD format
            double xutsec = sec;
            double xalt = alt;
            double xlat = glat;
            double xlon = glong;
            double xsfluxavg = f107a;
            double xsflux = f107;
            double[] xap = new double[7];
            for (int i = 0; i < 7; i++)
            {
                xap[i] = ap[i];
            }

            // Call the new subroutine
            double xtn, xtex;
            double[] xdn;
            MSISCalculator.Calculate(xday, xutsec, xalt, xlat, xlon, xsfluxavg, xsflux, xap,
                             out xtn, out xdn, out xtex);

            // Initialize output arrays
            d = new float[10];
            t = new float[2];

            // Convert the output arguments to the legacy format (MKS to CGS, re-order species)
            t[0] = (float)xtex;  // Exospheric temperature
            t[1] = (float)xtn;   // Temperature at altitude

            // Convert densities from m-3 to cm-3 (divide by 1e6)
            // Convert mass density from kg/m3 to g/cm3 (divide by 1e3)
            // Process all densities except mass density first
            for (int i = 1; i < 10; i++)
            {
                if (xdn[i] != Constants.DMissing)
                {
                    xdn[i] = xdn[i] * 1e-6; // m-3 to cm-3
                }
            }
            // Special handling for mass density (index 0)
            if (xdn[0] != Constants.DMissing)
            {
                xdn[0] = xdn[0] * 1e-3; // kg/m3 to g/cm3
            }

            // Re-order species to legacy format
            // MSIS 2.1 internal order: [0]=mass, [1]=N2, [2]=O2, [3]=O, [4]=He, [5]=H, [6]=Ar, [7]=N, [8]=AnomalousO, [9]=NO
            // Legacy output order: [0]=He, [1]=O, [2]=N2, [3]=O2, [4]=Ar, [5]=mass, [6]=H, [7]=N, [8]=AnomalousO, [9]=NO
            d[0] = (float)xdn[4];  // He (internal index 4 -> legacy index 0)
            d[1] = (float)xdn[3];  // O  (internal index 3 -> legacy index 1)
            d[2] = (float)xdn[1];  // N2 (internal index 1 -> legacy index 2)
            d[3] = (float)xdn[2];  // O2 (internal index 2 -> legacy index 3)
            d[4] = (float)xdn[6];  // Ar (internal index 6 -> legacy index 4)
            d[5] = (float)xdn[0];  // Mass density (internal index 0 -> legacy index 5)
            d[6] = (float)xdn[5];  // H  (internal index 5 -> legacy index 6)
            d[7] = (float)xdn[7];  // N  (internal index 7 -> legacy index 7)
            d[8] = (float)xdn[8];  // Anomalous O (internal index 8 -> legacy index 8)
            d[9] = (float)xdn[9];  // NO (internal index 9 -> legacy index 9)
        }
    }
}