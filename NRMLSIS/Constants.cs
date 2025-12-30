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
// MSIS_Constants: Contains constants and hardwired parameters
// **************************************************************************************************

using System;

namespace NRLMSIS
{
    /// <summary>
    /// Constants and hardwired parameters for NRLMSIS 2.1
    /// </summary>
    public static class Constants
    {
        // Floating Point Precision
        // Note: C# uses double (64-bit) by default. Use float for single precision if needed.
#if SINGLE_PRECISION
        public const float DMissing = 9.999e-38f;
#else
        public const double DMissing = 9.999e-38;
#endif

        // Trigonometric constants
        public const double Pi = 3.1415926535897932384626433832795;
        public const double Deg2Rad = Pi / 180.0;
        public const double Doy2Rad = 2.0 * Pi / 365.0;
        public const double Lst2Rad = Pi / 12.0;
        public static readonly double Tanh1 = Math.Tanh(1.0);

        // Thermodynamic constants
        // Boltzmann constant (CODATA 2018) (J/kg)
        public const double KB = 1.380649e-23;
        
        // Avogadro constant (CODATA 2018)
        public const double NA = 6.02214076e23;
        
        // Reference gravity (CIMO Guide 2014) (m/s^2)
        public const double G0 = 9.80665;

        // Species molecular masses (kg/molecule) (CIPM 2007)
        // Index: 0=Mass density (dummy), 1=N2, 2=O2, 3=O, 4=He, 5=H, 6=Ar, 7=N, 8=Anomalous O, 9=NO
        // NOTE: Array is 0-based but logically represents Fortran's 1-based indexing
        // Access with species number directly: SpecMass[2] for N2, SpecMass[3] for O2, etc.
        public static readonly double[] SpecMass = new double[]
        {
            0.0,                        // Index 0: Mass density (dummy value)
            0.0,                        // Index 1: Unused (to maintain Fortran indexing)
            28.0134 / (1.0e3 * NA),     // Index 2: N2
            31.9988 / (1.0e3 * NA),     // Index 3: O2
            31.9988 / 2.0 / (1.0e3 * NA), // Index 4: O
            4.0 / (1.0e3 * NA),         // Index 5: He
            1.0 / (1.0e3 * NA),         // Index 6: H
            39.948 / (1.0e3 * NA),      // Index 7: Ar
            28.0134 / 2.0 / (1.0e3 * NA), // Index 8: N
            31.9988 / 2.0 / (1.0e3 * NA), // Index 9: Anomalous O
            (28.0134 + 31.9988) / 2.0 / (1.0e3 * NA) // Index 10: NO
        };

        // Dry air mean mass in fully mixed atmosphere (CIPM 2007)
        public const double Mbar = 28.96546 / (1.0e3 * NA); // kg/molecule

        // Dry air log volume mixing ratios (CIPM 2007)
        // NOTE: Array is 0-based but logically represents Fortran's 1-based indexing
        public static readonly double[] LnVmr = new double[]
        {
            Math.Log(1.0),        // Index 0: Mass density (dummy value)
            Math.Log(1.0),        // Index 1: Unused
            Math.Log(0.780848),   // Index 2: N2
            Math.Log(0.209390),   // Index 3: O2
            Math.Log(1.0),        // Index 4: O (dummy value)
            Math.Log(0.0000052),  // Index 5: He
            Math.Log(1.0),        // Index 6: H (dummy value)
            Math.Log(0.009332),   // Index 7: Ar
            Math.Log(1.0),        // Index 8: N (dummy value)
            Math.Log(1.0),        // Index 9: Anomalous O (dummy value)
            Math.Log(1.0)         // Index 10: NO (dummy value)
        };

        // Natural log of global average surface pressure (Pa)
        public const double LnP0 = 11.515614;

        // Derived constants
        public const double G0DivKB = G0 / KB * 1.0e3;           // K/(kg km)
        public const double MbarG0DivKB = Mbar * G0 / KB * 1.0e3; // K/km

        // Vertical profile parameters
        public const int NSpec = 11;  // Number of species including temperature
        public const int Nd = 27;     // Number of temperature profile nodes
        public const int P = 4;       // Spline order
        public const int Nl = Nd - P; // Last temperature profile level index
        public const int Nls = 9;     // Last parameter index for each species (excluding O, NO splines)
        
        public const double BwAlt = 122.5;      // Reference geopotential height for Bates Profile
        public const double ZetaF = 70.0;       // Fully mixed below this, uses constant mixing ratios
        public const double ZetaB = BwAlt;      // Bates Profile above this altitude
        public const double ZetaA = 85.0;       // Default reference height for active minor species
        public const double ZetaGamma = 100.0;  // Reference height of tanh taper
        public const double HGamma = 1.0 / 30.0; // Inverse scale height of tanh taper

        // Nodes for temperature profile splines (0-based indexing in C#)
        public static readonly double[] NodesTN = new double[]
        {
            -15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0,
            55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 92.5, 102.5, 112.5, 122.5, 132.5, 142.5,
            152.5, 162.5, 172.5
        };

        public const int IzFmx = 13;  // Fully mixed below this spline index
        public const int IzFx = 14;   // Spline index at zetaF
        public const int IzAx = 17;   // Spline index at zetaA
        public const int ItEx = Nl;   // Index of Bates exospheric temperature
        public const int ItGb0 = Nl - 1; // Index of Bates temperature gradient at lower boundary
        public const int ItB0 = Nl - 2;  // Index of Bates temperature at lower boundary

        // O1 Spline parameters
        public const int NdO1 = 13;
        public const int NsplO1 = NdO1 - 5; // Number of unconstrained spline parameters for O1
        
        public static readonly double[] NodesO1 = new double[]
        {
            35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 92.5, 102.5, 112.5
        };
        
        public const double ZetaRefO1 = ZetaA; // Joining height for O1 splines

        // NO Spline parameters
        public const int NdNO = 13;
        public const int NsplNO = NdNO - 5; // Number of unconstrained spline parameters for NO
        
        public static readonly double[] NodesNO = new double[]
        {
            47.5, 55.0, 62.5, 70.0, 77.5, 85.0, 92.5, 100.0, 107.5, 115.0, 122.5, 130.0, 137.5, 145.0
        };
        
        public const double ZetaRefNO = ZetaB; // Joining height for NO splines

        // C2 Continuity matrix for temperature (3x3)
        // NOTE: Fortran uses column-major order; this is the correct transpose for C# row-major
        public static readonly double[,] C2Tn = new double[,]
        {
            {  1.0,   1.0,   1.0 },
            { -10.0,  0.0,  10.0 },
            {  33.333333333333336, -16.666666666666668,  33.333333333333336 }
        };

        // C1 Continuity for O1 (2x2)
        // NOTE: Fortran uses column-major order; this is the correct transpose for C# row-major
        public static readonly double[,] C1O1 = new double[,]
        {
            {  1.75,              -1.624999900076852 },
            { -2.916666573405061,  21.458332647194382 }
        };

        public static readonly double[] C1O1Adj = new double[]
        {
            0.257142857142857, -0.102857142686844
        };

        // C1 Continuity for NO (2x2)
        // NOTE: Fortran uses column-major order; this is the correct transpose for C# row-major
        public static readonly double[,] C1NO = new double[,]
        {
            {  1.5,  0.0 },
            { -3.75, 15.0 }
        };

        public static readonly double[] C1NOAdj = new double[]
        {
            0.166666666666667, -0.066666666666667
        };

        // Anomalous Oxygen parameters (legacy profile from NRLMSISE-00)
        public const double ZetaRefOA = ZetaB;
        public const double TOA = 4000.0; // Temperature of anomalous oxygen density (K)
        public static readonly double HOA = (KB * TOA) / ((16.0 / (1.0e3 * NA)) * G0) * 1.0e-3; // Scale height (km)

        // Horizontal and time-dependent basis function (gfn) parameters
        public const int MaxNbf = 512;   // Number of basis functions to be allocated
        public const int MaxN = 6;       // Maximum latitude (Legendre) spectral degree
        public const int MaxL = 3;       // Maximum local time (tidal) spectral order
        public const int MaxM = 2;       // Maximum longitude (stationary planetary wave) order
        public const int MaxS = 2;       // Maximum day of year (intra-annual) Fourier order
        public const int AMaxN = 6;      // Maximum Legendre degree in time independent terms
        public const int AMaxS = 2;      // Maximum intra-annual order in zonal mean terms
        public const int TMaxL = 3;      // Maximum tidal order used
        public const int TMaxN = 6;      // Maximum Legendre degree coupled with tides
        public const int TMaxS = 2;      // Maximum intra-annual order coupled with tides
        public const int PMaxM = 2;      // Maximum stationary planetary wave order used
        public const int PMaxN = 6;      // Maximum Legendre degree coupled with SPW
        public const int PMaxS = 2;      // Maximum intra-annual order coupled with SPW
        public const int NSfx = 5;       // Number of linear solar flux terms
        public const int NSfxMod = 5;    // Number of nonlinear modulating solar flux terms
        public const int NMag = 54;      // Number of terms in geomagnetic parameterization
        public const int NUt = 12;       // Number of terms in UT parameterization
        
        public const int CTimeInd = 0;
        public const int CIntAnn = CTimeInd + (AMaxN + 1);
        public const int CTide = CIntAnn + ((AMaxN + 1) * 2 * AMaxS);
        public const int CSpw = CTide + (4 * TMaxS + 2) * (TMaxL * (TMaxN + 1) - (TMaxL * (TMaxL + 1)) / 2);
        public const int CSfx = CSpw + (4 * PMaxS + 2) * (PMaxM * (PMaxN + 1) - (PMaxM * (PMaxM + 1)) / 2);
        public const int CExtra = CSfx + NSfx;
        public const int Mbf = 383;
        public const int CNonLin = Mbf + 1;
        public const int CSfxMod = CNonLin;
        public const int CMag = CSfxMod + NSfxMod;
        public const int CUt = CMag + NMag;

        // Weights for calculation log pressure spline coefficients from temperature coefficients
        public static readonly double[] Gwht = new double[]
        {
            5.0 / 24.0, 55.0 / 24.0, 55.0 / 24.0, 5.0 / 24.0
        };

        // Constants for analytical integration by parts of hydrostatic piecewise effective mass profile
        public static readonly double[] WBeta;
        public static readonly double[] WGamma;

        // Non-zero bspline values at zetaB (5th and 6th order)
        public static readonly double[] S5ZetaB = new double[]
        {
            0.041666666666667, 0.458333333333333, 0.458333333333333, 0.041666666666667
        };

        public static readonly double[] S6ZetaB = new double[]
        {
            0.008771929824561, 0.216228070175439, 0.550000000000000, 0.216666666666667, 0.008333333333333
        };

        // Weights for calculating temperature gradient at zetaA
        public static readonly double[] WghtAxdz = new double[]
        {
            -0.102857142857, 0.0495238095238, 0.053333333333
        };

        // Non-zero bspline values at zetaA (4th, 5th and 6th order)
        public static readonly double[] S4ZetaA = new double[]
        {
            0.257142857142857, 0.653968253968254, 0.088888888888889
        };

        public static readonly double[] S5ZetaA = new double[]
        {
            0.085714285714286, 0.587590187590188, 0.313020313020313, 0.013675213675214
        };

        public static readonly double[] S6ZetaA = new double[]
        {
            0.023376623376623, 0.378732378732379, 0.500743700743701, 0.095538448479625, 0.001608848667672
        };

        // Non-zero bspline values at zetaF (4th and 5th order)
        public static readonly double[] S4ZetaF = new double[]
        {
            0.166666666666667, 0.666666666666667, 0.166666666666667
        };

        public static readonly double[] S5ZetaF = new double[]
        {
            0.041666666666667, 0.458333333333333, 0.458333333333333, 0.041666666666667
        };

        // Non-zero bspline values at zeta=0 (5th order)
        public static readonly double[] S5Zeta0 = new double[]
        {
            0.458333333333333, 0.458333333333333, 0.041666666666667
        };

        // Static constructor to initialize computed arrays
        static Constants()
        {
            // Initialize WBeta: (nodesTN(4:nd) - nodesTN(0:nl)) / 4.0
            WBeta = new double[Nl + 1];
            for (int i = 0; i <= Nl; i++)
            {
                WBeta[i] = (NodesTN[i + 4] - NodesTN[i]) / 4.0;
            }

            // Initialize WGamma: (nodesTN(5:nd+1) - nodesTN(0:nl)) / 5.0
            WGamma = new double[Nl + 1];
            for (int i = 0; i <= Nl; i++)
            {
                WGamma[i] = (NodesTN[i + 5] - NodesTN[i]) / 5.0;
            }
        }
    }

    // References:
    // CODATA Internationally recommended 2018 values of the fundamental physical constants.
    //   https://pml.nist.gov/cuu/Constants/; https://pml.nist.gov/cuu/pdf/wallet_2018.pdf
    // Picard, A., Davis, R. S., Glaeser, M., and Fujii, K. (2007). Revised formula for the density of
    //   air (CIPM 2007). Metrologia 45, 149–155. doi:10.1088/0026-1394/45/2/004
    // World Meteorological Organization (2014). WMO guide to meteorological instruments and methods of observation
    //   (the CIMO Guide). Part I, Chapter 12. https://www.wmo.int/pages/prog/www/IMOP/CIMO-Guide.html
}