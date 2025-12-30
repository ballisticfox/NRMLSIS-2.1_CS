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
// MSISINIT: Initialization of MSIS parameters, switches, and options.
//
//     PREREQUISITES:
//       MSIS binary parameter file (msis21.parm)
//
//     CALLING SEQUENCE:
//       MsisInit([OPTIONAL ARGUMENTS])
//  
//     OPTIONAL ARGUMENTS:
//       parmPath        File path pointing to the MSIS parameter file.
//                         Default: "" (current directory)
//       parmFile        Name of MSIS parameter file.
//                         Default: "msis21.parm"
//       switchGfn       Boolean array of 512 switches for individual terms. For
//                         advanced users.
//                         Default values: true (all switches on)
//       switchLegacy    Float array (1:25) of legacy switches that
//                         control groups of terms:
//                            1 - F10.7
//                            2 - Time independent
//                            3 - Symmetrical annual
//                            4 - Symmetrical semiannual
//                            5 - Asymmetrical annual
//                            6 - Asymmetrical semiannual
//                            7 - Diurnal
//                            8 - Semidiurnal
//                            9 - Geomagnetic activity:
//                                  1.0 = Daily Ap mode
//                                 -1.0 = Storm-time ap mode
//                           10 - All UT/long effects
//                           11 - Longitudinal
//                           12 - UT and mixed UT/long
//                           13 - Mixed Ap/UT/long
//                           14 - Terdiurnal
//                           15-25 - Not used in NRLMSIS 2.1
//                         For all switches:
//                           0.0 = Off
//                           1.0 = On
//                           2.0 = Main effects off, cross terms on
//                         Default values: 1.0
//       lZAltType       Boolean flag for altitude input type:
//                         true = Geodetic altitude (km)
//                         false = Geopotential height (km)
//                         Default: true (Geodetic altitude)
//       lSpecSelect     Boolean array (1:10) flagging which densities to 
//                         calculate.
//                         true = Calculate, false = Do not calculate
//                            1 - Mass density
//                            2 - N2
//                            3 - O2
//                            4 - O
//                            5 - He
//                            6 - H
//                            7 - Ar
//                            8 - N
//                            9 - Anomalous O
//                           10 - NO
//                         Default values: true
//       lMassInclude    Boolean array (1:10) flagging which species to include
//                         in mass density calculation. Same ordering as 
//                         lSpecSelect.
//                         Default values: true
//       lN2Msis00       Boolean flag for retrieving NRLMSISE-00 upper
//                         thermospheric N2 variation. See paper for details.
//                           false: Thermospheric N2 determined entirely by
//                             temperature profile and the constant mixing ratio
//                             of N2 in the lower atmosphere. 
//                           true: Upper thermospheric N2 relaxes to NRLMSISE-00
//                             values.
//                         Default: false
//
// ===========================================================================

using System;
using System.IO;
using System.Linq;

namespace NRLMSIS
{
    /// <summary>
    /// Parameter subset structure for MSIS model
    /// </summary>
    public class BasisSubset
    {
        public string Name { get; set; } = "";
        public int Bl { get; set; }
        public int Nl { get; set; }
        public double[,] Beta { get; set; } = new double[0, 0];
        public bool[,] Active { get; set; } = new bool[0, 0];
        public int[,] Fitb { get; set; } = new int[0, 0];
    }

    /// <summary>
    /// Initialization and parameter management for NRLMSIS 2.1
    /// </summary>
    public static class MsisInit
    {
        // Model flags
        public static bool InitFlag { get; private set; } = false;
        public static bool HaveParmSpace { get; private set; } = false;
        public static bool ZAltFlag { get; set; } = true;
        // SpecFlag and MassFlag: 0-based arrays [0:9] representing species 1-10
        // SpecFlag[0] = species 1 (mass density), SpecFlag[1] = species 2 (N2), ..., SpecFlag[9] = species 10 (NO)
        public static bool[] SpecFlag { get; set; } = Enumerable.Repeat(true, MsisConstants.NSpec - 1).ToArray();
        public static bool[] MassFlag { get; set; } = Enumerable.Repeat(true, MsisConstants.NSpec - 1).ToArray();
        public static bool N2RFlag { get; set; } = false;
        public static bool[] Zsfx { get; set; } = new bool[MsisConstants.Mbf + 1];
        public static bool[] Tsfx { get; set; } = new bool[MsisConstants.Mbf + 1];
        public static bool[] Psfx { get; set; } = new bool[MsisConstants.Mbf + 1];
        public static bool[] Smod { get; set; } = new bool[MsisConstants.Nl + 1];
        public static bool[] Swg { get; set; } = Enumerable.Repeat(true, MsisConstants.MaxNbf).ToArray();
        public static double[] MassWgt { get; set; } = new double[MsisConstants.NSpec - 1];
        public static float[] SwLeg { get; set; } = Enumerable.Repeat(1.0f, 25).ToArray();
        public static float[] Swc { get; set; } = new float[25];
        public static float[] Sav { get; set; } = new float[25];

        // Model parameter arrays (initialized with empty instances, populated by InitParmSpace)
        public static BasisSubset TN { get; private set; } = new BasisSubset();
        public static BasisSubset PR { get; private set; } = new BasisSubset();
        public static BasisSubset N2 { get; private set; } = new BasisSubset();
        public static BasisSubset O2 { get; private set; } = new BasisSubset();
        public static BasisSubset O1 { get; private set; } = new BasisSubset();
        public static BasisSubset HE { get; private set; } = new BasisSubset();
        public static BasisSubset H1 { get; private set; } = new BasisSubset();
        public static BasisSubset AR { get; private set; } = new BasisSubset();
        public static BasisSubset N1 { get; private set; } = new BasisSubset();
        public static BasisSubset OA { get; private set; } = new BasisSubset();
        public static BasisSubset NO { get; private set; } = new BasisSubset();
        public static int NVertParm { get; private set; }

        // Reciprocal node difference arrays (constant values needed for B-spline calculations)
        public static double[,] EtaTN { get; set; } = new double[31, 5]; // [0:30, 2:6] -> [0:30, 0:4]
        public static double[,] EtaO1 { get; set; } = new double[31, 5];
        public static double[,] EtaNO { get; set; } = new double[31, 5];

        // C1 constraint terms for O and NO related to the tapered logistic correction
        public static double HRFactO1Ref { get; private set; }
        public static double DHRFactO1Ref { get; private set; }
        public static double HRFactNORef { get; private set; }
        public static double DHRFactNORef { get; private set; }

        // ==================================================================================================
        // MSISINIT: Entry point for initializing model and loading parameters
        // ==================================================================================================
        /// <summary>
        /// Initialize MSIS model and load parameters
        /// </summary>
        public static void MsisInitialize(
            string parmPath = "",
            string parmFile = "msis21.parm",
            bool[]? switchGfn = null,
            float[]? switchLegacy = null,
            bool? lZAltType = null,
            bool[]? lSpecSelect = null,
            bool[]? lMassInclude = null,
            bool? lN2Msis00 = null)
        {
            string parmPath1 = parmPath;
            string parmFile1 = parmFile;

            // Initialize model parameter space
            if (!HaveParmSpace) InitParmSpace();

            // Load parameter set
            LoadParmSet(Path.Combine(parmPath1, parmFile1));

            // Set switches
            Swg = Enumerable.Repeat(true, MsisConstants.MaxNbf).ToArray();
            SwLeg = Enumerable.Repeat(1.0f, 25).ToArray();

            if (switchGfn != null)
            {
                Swg = (bool[])switchGfn.Clone();
            }
            else if (switchLegacy != null)
            {
                SwLeg = (float[])switchLegacy.Clone();
                TSelec(SwLeg);
            }

            // Input altitude type flag
            ZAltFlag = lZAltType ?? true;

            // Species flags for number density and mass density
            if (lSpecSelect != null)
            {
                SpecFlag = (bool[])lSpecSelect.Clone();
            }
            else
            {
                SpecFlag = Enumerable.Repeat(true, MsisConstants.NSpec - 1).ToArray();
            }

            if (SpecFlag[0])
            {
                if (lMassInclude != null)
                {
                    MassFlag = (bool[])lMassInclude.Clone();
                }
                else
                {
                    MassFlag = Enumerable.Repeat(true, MsisConstants.NSpec - 1).ToArray();
                }
            }
            else
            {
                MassFlag = Enumerable.Repeat(false, MsisConstants.NSpec - 1).ToArray();
            }

            // Where massflag is true, specflag must also be true
            for (int i = 0; i < MassFlag.Length; i++)
            {
                if (MassFlag[i]) SpecFlag[i] = true;
            }

            // Calculate mass weights
            // Fortran: masswgt(1:nspec-1) where(massflag) masswgt = 1.0
            //          masswgt(1) = 0.0 (mass density)
            //          masswgt = masswgt * specmass
            //          masswgt(10) = 0.0 (NO)
            // C#: MassWgt[0:9] representing species 1-10
            //     MassWgt[0] = species 1 (mass), MassWgt[1] = species 2 (N2), etc.
            MassWgt = new double[MsisConstants.NSpec - 1];
            for (int i = 0; i < MassFlag.Length; i++)
            {
                MassWgt[i] = MassFlag[i] ? 1.0 : 0.0;
            }
            MassWgt[0] = 0.0; // Mass density itself doesn't contribute

            // Multiply by species masses
            // MassWgt and MassFlag are [0:9] representing species [1:10]
            // SpecMass is [0:10] indexed by species number
            // MassWgt[0] = species 1, multiply by SpecMass[1] (but it's 0 anyway)
            // MassWgt[1] = species 2 (N2), multiply by SpecMass[2]
            // MassWgt[i] = species i+1, multiply by SpecMass[i+1]
            for (int i = 0; i < MassWgt.Length; i++)
            {
                MassWgt[i] *= MsisConstants.SpecMass[i + 1];
            }
            MassWgt[9] = 0.0; // NO (species 10) doesn't contribute to mass density

            // Flag for retrieving NRLMSISE-00 thermospheric N2 variations
            N2RFlag = lN2Msis00 ?? false;

            // Set model initialization flag
            InitFlag = true;
        }

        // ==================================================================================================
        // INITPARMSPACE: Initialize and allocate the model parameter space
        // ==================================================================================================
        private static void InitParmSpace()
        {
            // Vertical parameter counter (number of vertical parameters in the parameter file)
            NVertParm = 0;

            // Model formulation parameter subsets
            TN = InitSubset(0, MsisConstants.Nl, MsisConstants.MaxNbf, "TN");
            PR = InitSubset(0, MsisConstants.Nl, MsisConstants.MaxNbf, "PR");
            N2 = InitSubset(0, MsisConstants.Nls, MsisConstants.MaxNbf, "N2");
            O2 = InitSubset(0, MsisConstants.Nls, MsisConstants.MaxNbf, "O2");
            O1 = InitSubset(0, MsisConstants.Nls + MsisConstants.NsplO1, MsisConstants.MaxNbf, "O1");
            HE = InitSubset(0, MsisConstants.Nls, MsisConstants.MaxNbf, "HE");
            H1 = InitSubset(0, MsisConstants.Nls, MsisConstants.MaxNbf, "H1");
            AR = InitSubset(0, MsisConstants.Nls, MsisConstants.MaxNbf, "AR");
            N1 = InitSubset(0, MsisConstants.Nls, MsisConstants.MaxNbf, "N1");
            OA = InitSubset(0, MsisConstants.Nls, MsisConstants.MaxNbf, "OA");
            NO = InitSubset(0, MsisConstants.Nls + MsisConstants.NsplNO, MsisConstants.MaxNbf, "NO");

            // Add the surface pressure parameter to the vertical parameter counter
            NVertParm = NVertParm + 1;

            // Set solar flux modulation flags
            Zsfx = new bool[MsisConstants.Mbf + 1];
            Tsfx = new bool[MsisConstants.Mbf + 1];
            Psfx = new bool[MsisConstants.Mbf + 1];

            // F1, solar flux modulation of the zonal mean asymmetric annual terms
            Zsfx[9] = Zsfx[10] = true;    // Pl(1,0) annual terms
            Zsfx[13] = Zsfx[14] = true;   // Pl(3,0) annual terms
            Zsfx[17] = Zsfx[18] = true;   // Pl(5,0) annual terms

            // F2, solar flux modulation of the tides
            for (int i = MsisConstants.CTide; i < MsisConstants.CSpw; i++)
            {
                Tsfx[i] = true;
            }

            // F3, solar flux modulation of stationary planetary wave 1
            for (int i = MsisConstants.CSpw; i < MsisConstants.CSpw + 60; i++)
            {
                Psfx[i] = true;
            }

            // Calculate reciprocal node difference arrays
            for (int k = 2; k <= 6; k++)
            {
                for (int j = 0; j <= MsisConstants.Nl; j++)
                {
                    EtaTN[j, k - 2] = 1.0 / (MsisConstants.NodesTN[j + k - 1] - MsisConstants.NodesTN[j]);
                }
            }

            for (int k = 2; k <= 4; k++)
            {
                for (int j = 0; j <= MsisConstants.NdO1 - k + 1; j++)
                {
                    EtaO1[j, k - 2] = 1.0 / (MsisConstants.NodesO1[j + k - 1] - MsisConstants.NodesO1[j]);
                }
                for (int j = 0; j <= MsisConstants.NdNO - k + 1; j++)
                {
                    EtaNO[j, k - 2] = 1.0 / (MsisConstants.NodesNO[j + k - 1] - MsisConstants.NodesNO[j]);
                }
            }

            // Calculate C1 constraint terms for O and NO related to the tapered logistic correction
            double gammaTerm0 = Math.Tanh((MsisConstants.ZetaRefO1 - MsisConstants.ZetaGamma) * MsisConstants.HGamma);
            HRFactO1Ref = 0.5 * (1.0 + gammaTerm0);
            DHRFactO1Ref = (1.0 - (MsisConstants.ZetaRefO1 - MsisConstants.ZetaGamma) * (1.0 - gammaTerm0) * MsisConstants.HGamma) / HRFactO1Ref;

            gammaTerm0 = Math.Tanh((MsisConstants.ZetaRefNO - MsisConstants.ZetaGamma) * MsisConstants.HGamma);
            HRFactNORef = 0.5 * (1.0 + gammaTerm0);
            DHRFactNORef = (1.0 - (MsisConstants.ZetaRefNO - MsisConstants.ZetaGamma) * (1.0 - gammaTerm0) * MsisConstants.HGamma) / HRFactNORef;

            // Set parameter space initialization flag
            HaveParmSpace = true;
        }

        // --------------------------------------------------------------------------------------------------
        // INITSUBSET: Initialize and allocate a parameter subset
        // --------------------------------------------------------------------------------------------------
        private static BasisSubset InitSubset(int bl, int nl, int maxnbf, string name)
        {
            var subset = new BasisSubset
            {
                Name = name,
                Bl = bl,
                Nl = nl,
                Beta = new double[maxnbf, nl - bl + 1],      // [0:maxnbf-1, bl:nl]
                Active = new bool[maxnbf, nl - bl + 1],
                Fitb = new int[maxnbf, nl - bl + 1]
            };

            // Increment vertical parameter counter except for pressure
            if (name != "PR")
            {
                NVertParm = NVertParm + nl - bl + 1;
            }

            return subset;
        }

        // ==================================================================================================
        // LOADPARMSET: Read in a parameter file
        // ==================================================================================================
        private static void LoadParmSet(string name)
        {
            // Check if file exists
            if (!File.Exists(name))
            {
                // Try with full path
                string currentDir = Directory.GetCurrentDirectory();
                string fullPath = Path.GetFullPath(name);

                throw new FileNotFoundException(
                    $"MSIS parameter set '{name}' not found.\n" +
                    $"Current directory: {currentDir}\n" +
                    $"Looking for: {fullPath}\n" +
                    $"Please ensure the file exists and has read permissions.");
            }

            // Read in parameter values into temporary double-precision array
            double[,] parmIn = new double[MsisConstants.MaxNbf, NVertParm];

            try
            {
                using (BinaryReader reader = new BinaryReader(File.Open(name, FileMode.Open, FileAccess.Read)))
                {
                    // Read parameters in column-major order (Fortran style)
                    for (int col = 0; col < NVertParm; col++)
                    {
                        for (int row = 0; row < MsisConstants.MaxNbf; row++)
                        {
                            parmIn[row, col] = reader.ReadDouble();
                        }
                    }
                }
            }
            catch (UnauthorizedAccessException ex)
            {
                throw new UnauthorizedAccessException(
                    $"Access denied to '{name}'. Please check file permissions.\n" +
                    $"On macOS/Linux, ensure the file has read permissions (chmod +r {name})", ex);
            }
            catch (IOException ex)
            {
                throw new IOException($"Error reading parameter file '{name}': {ex.Message}", ex);
            }

            // Transfer parameters to structures
            int i0 = 0;
            int i1 = TN.Nl - TN.Bl;
            CopyParameters(parmIn, TN.Beta, i0, i1, TN.Bl, TN.Nl);

            i0 = i1 + 1;
            i1 = i0;
            for (int j = 0; j < MsisConstants.MaxNbf; j++)
            {
                PR.Beta[j, PR.Bl] = parmIn[j, i0];
            }

            i0 = i1 + 1;
            i1 = i0 + N2.Nl - N2.Bl;
            CopyParameters(parmIn, N2.Beta, i0, i1, N2.Bl, N2.Nl);

            i0 = i1 + 1;
            i1 = i0 + O2.Nl - O2.Bl;
            CopyParameters(parmIn, O2.Beta, i0, i1, O2.Bl, O2.Nl);

            i0 = i1 + 1;
            i1 = i0 + O1.Nl - O1.Bl;
            CopyParameters(parmIn, O1.Beta, i0, i1, O1.Bl, O1.Nl);

            i0 = i1 + 1;
            i1 = i0 + HE.Nl - HE.Bl;
            CopyParameters(parmIn, HE.Beta, i0, i1, HE.Bl, HE.Nl);

            i0 = i1 + 1;
            i1 = i0 + H1.Nl - H1.Bl;
            CopyParameters(parmIn, H1.Beta, i0, i1, H1.Bl, H1.Nl);

            i0 = i1 + 1;
            i1 = i0 + AR.Nl - AR.Bl;
            CopyParameters(parmIn, AR.Beta, i0, i1, AR.Bl, AR.Nl);

            i0 = i1 + 1;
            i1 = i0 + N1.Nl - N1.Bl;
            CopyParameters(parmIn, N1.Beta, i0, i1, N1.Bl, N1.Nl);

            i0 = i1 + 1;
            i1 = i0 + OA.Nl - OA.Bl;
            CopyParameters(parmIn, OA.Beta, i0, i1, OA.Bl, OA.Nl);

            i0 = i1 + 1;
            i1 = i0 + NO.Nl - NO.Bl;
            CopyParameters(parmIn, NO.Beta, i0, i1, NO.Bl, NO.Nl);

            // Set solar flux modulation flags; if on for a given vertical parameter, then sfluxmod is called by tfnparm
            Smod = new bool[MsisConstants.Nl + 1];
            for (int iz = 0; iz <= MsisConstants.Nl; iz++)
            {
                if (TN.Beta[MsisConstants.CSfxMod + 0, iz - TN.Bl] != 0 ||
                    TN.Beta[MsisConstants.CSfxMod + 1, iz - TN.Bl] != 0 ||
                    TN.Beta[MsisConstants.CSfxMod + 2, iz - TN.Bl] != 0)
                {
                    Smod[iz] = true;
                }
            }

            // Compute log pressure spline coefficients from temperature spline coefficients
            PressParm();
        }

        // Helper method to copy parameters from input array to subset array
        private static void CopyParameters(double[,] source, double[,] dest, int srcColStart, int srcColEnd, int destBl, int destNl)
        {
            for (int j = 0; j < MsisConstants.MaxNbf; j++)
            {
                for (int col = 0; col <= srcColEnd - srcColStart; col++)
                {
                    dest[j, col] = source[j, srcColStart + col];
                }
            }
        }

        // ==================================================================================================
        // PRESSPARM: Compute log pressure spline coefficients from temperature spline coefficients
        // ==================================================================================================
        private static void PressParm()
        {
            // Integrate pressure on nodes up to the last fully mixed level
            for (int j = 0; j <= MsisConstants.Mbf; j++)
            {
                double lnz = 0.0;
                for (int b = 0; b <= 3; b++)
                {
                    lnz = lnz + TN.Beta[j, b] * MsisConstants.Gwht[b] * MsisConstants.MbarG0DivKB;
                }
                PR.Beta[j, 1] = -lnz;

                for (int iz = 1; iz <= MsisConstants.IzFmx; iz++)
                {
                    lnz = 0.0;
                    for (int b = 0; b <= 3; b++)
                    {
                        lnz = lnz + TN.Beta[j, iz + b] * MsisConstants.Gwht[b] * MsisConstants.MbarG0DivKB;
                    }
                    PR.Beta[j, iz + 1] = PR.Beta[j, iz] - lnz;
                }
            }
        }

        // ==================================================================================================
        // TSELEC: Legacy switches and mapping to new switches
        // ==================================================================================================
        private static void TSelec(float[] sv)
        {
            // Set cross-terms flags
            // Note: sv is 0-based in C#, but logically represents Fortran's 1:25 array
            for (int i = 0; i < 25; i++)
            {
                Sav[i] = sv[i];
                SwLeg[i] = sv[i] % 2.0f;
                if (Math.Abs(sv[i]) == 1.0f || Math.Abs(sv[i]) == 2.0f)
                {
                    Swc[i] = 1.0f;
                }
                else
                {
                    Swc[i] = 0.0f;
                }
            }

            // Main effects
            // Note: In the comments below, switch indices refer to the logical 1-based Fortran indices
            Swg[0] = true; // Global term must be on
            SetRange(Swg, MsisConstants.CSfx, MsisConstants.CSfx + MsisConstants.NSfx - 1, SwLeg[0] == 1.0f); // Solar flux (switch 1)
            Swg[310] = (SwLeg[0] == 1.0f); // Solar flux (truncated quadratic F10.7a function) (switch 1)
            SetRange(Swg, 1, 6, SwLeg[1] == 1.0f); // Time independent (switch 2)
            SetRange(Swg, 304, 305, SwLeg[1] == 1.0f); // Time independent (extra, F10.7a modulated terms) (switch 2)
            SetRange(Swg, 311, 312, SwLeg[1] == 1.0f); // Time independent (extra, truncated quadratic F10.7a modulated terms) (switch 2)
            SetRange(Swg, 313, 314, SwLeg[1] == 1.0f); // Time independent (extra, dF10.7 modulated terms) (switch 2)
            SetIndices(Swg, new[] { 7, 8, 11, 12, 15, 16, 19, 20 }, SwLeg[2] == 1.0f); // Symmetric annual (switch 3)
            SetRange(Swg, 306, 307, SwLeg[2] == 1.0f); // Global AO (extra, solar-flux modulated terms) (switch 3)
            SetIndices(Swg, new[] { 21, 22, 25, 26, 29, 30, 33, 34 }, SwLeg[3] == 1.0f); // Symmetric semiannual (switch 4)
            SetRange(Swg, 308, 309, SwLeg[3] == 1.0f); // Global SAO (extra, solar-flux modulated terms) (switch 4)
            SetIndices(Swg, new[] { 9, 10, 13, 14, 17, 18 }, SwLeg[4] == 1.0f); // Asymmetric annual (switch 5)
            SetIndices(Swg, new[] { 23, 24, 27, 28, 31, 32 }, SwLeg[5] == 1.0f); // Asymmetric semiannual (switch 6)
            SetRange(Swg, 35, 94, SwLeg[6] == 1.0f); // Diurnal (switch 7)
            SetRange(Swg, 300, 303, SwLeg[6] == 1.0f); // Solar zenith angle (switch 7)
            SetRange(Swg, 95, 144, SwLeg[7] == 1.0f); // Semidiurnal (switch 8)
            SetRange(Swg, 145, 184, SwLeg[13] == 1.0f); // Terdiurnal (switch 14)

            // Geomagnetic activity mode master switch
            Swg[MsisConstants.CMag] = Swg[MsisConstants.CMag + 1] = false;
            if ((SwLeg[8] > 0) || (SwLeg[12] == 1)) // switch 9 or switch 13
            {
                Swg[MsisConstants.CMag] = Swg[MsisConstants.CMag + 1] = true; // Daily mode master switch
            }
            if (SwLeg[8] < 0) // switch 9
            {
                Swg[MsisConstants.CMag] = false;
                Swg[MsisConstants.CMag + 1] = true; // Storm-time mode master switch
            }

            SetRange(Swg, MsisConstants.CMag + 2, MsisConstants.CMag + 12, SwLeg[8] == 1.0f); // Daily geomagnetic activity terms (switch 9)
            SetRange(Swg, MsisConstants.CMag + 28, MsisConstants.CMag + 40, SwLeg[8] == -1.0f); // Storm-time geomagnetic activity terms (switch 9)
            SetRange(Swg, MsisConstants.CSpw, MsisConstants.CSfx - 1, (SwLeg[10] == 1.0f) && (SwLeg[9] == 1.0f)); // Longitudinal (switches 11 and 10)
            SetRange(Swg, MsisConstants.CUt, MsisConstants.CUt + MsisConstants.NUt - 1, (SwLeg[11] == 1.0f) && (SwLeg[9] == 1.0f)); // UT/Lon (switches 12 and 10)
            SetRange(Swg, MsisConstants.CMag + 13, MsisConstants.CMag + 25, (SwLeg[12] == 1.0f) && (SwLeg[9] == 1.0f)); // Mixed UT/Lon/Geomag (Daily mode terms) (switches 13 and 10)
            SetRange(Swg, MsisConstants.CMag + 41, MsisConstants.CMag + 53, (SwLeg[12] == 1.0f) && (SwLeg[9] == 1.0f)); // Mixed UT/Lon/Geomag (Storm-time mode terms) (switches 13 and 10)

            // Cross terms
            SetRange(Swg, MsisConstants.CSfxMod, MsisConstants.CSfxMod + MsisConstants.NSfxMod - 1, Swc[0] == 1.0f); // Solar activity modulation
            if (Swc[0] == 0)
            {
                SetRange(Swg, 302, 303, false); // Solar zenith angle
                SetRange(Swg, 304, 305, false); // Time independent
                SetRange(Swg, 306, 307, false); // Global AO
                SetRange(Swg, 308, 309, false); // Global SAO
                SetRange(Swg, 311, 314, false); // Time independent
                Swg[447] = false; // UT/Lon
                Swg[454] = false; // UT/Lon
            }
            if (Swc[1] == 0) // Time independent (latitude terms)
            {
                SetRange(Swg, 9, 20, false); // AO
                SetRange(Swg, 23, 34, false); // SAO
                SetRange(Swg, 35, 184, false); // All tides
                SetRange(Swg, 185, 294, false); // All SPW
                SetRange(Swg, 392, 414, false); // Daily geomagnetic activity
                SetRange(Swg, 420, 442, false); // Storm-time geomagnetic activity
                SetRange(Swg, 449, 453, false); // UT/Lon
            }
            if (Swc[2] == 0) // Symmetric annual
            {
                SetRange(Swg, 201, 204, false); // SPW1 (2,1)
                SetRange(Swg, 209, 212, false); // SPW1 (4,1)
                SetRange(Swg, 217, 220, false); // SPW1 (6,1)
                SetRange(Swg, 255, 258, false); // SPW2 (2,2)
                SetRange(Swg, 263, 266, false); // SPW2 (4,2)
                SetRange(Swg, 271, 274, false); // SPW2 (6,2)
                SetRange(Swg, 306, 307, false); // Global AO solar flux modulation
            }
            if (Swc[3] == 0) // Symmetric semiannual
            {
                SetRange(Swg, 225, 228, false); // SPW1 (2,1)
                SetRange(Swg, 233, 236, false); // SPW1 (4,1)
                SetRange(Swg, 241, 244, false); // SPW1 (6,1)
                SetRange(Swg, 275, 278, false); // SPW2 (2,2)
                SetRange(Swg, 283, 286, false); // SPW2 (4,2)
                SetRange(Swg, 291, 294, false); // SPW2 (6,2)
                SetRange(Swg, 308, 309, false); // Global SAO solar flux modulation
            }
            if (Swc[4] == 0) // Asymmetric annual
            {
                SetRange(Swg, 47, 50, false); // Diurnal (1,1)
                SetRange(Swg, 51, 54, false); // Diurnal (2,1)
                SetRange(Swg, 55, 58, false); // Diurnal (3,1)
                SetRange(Swg, 59, 62, false); // Diurnal (4,1)
                SetRange(Swg, 63, 66, false); // Diurnal (5,1)
                SetRange(Swg, 67, 70, false); // Diurnal (6,1)
                SetRange(Swg, 105, 108, false); // Semidiurnal (2,2)
                SetRange(Swg, 109, 112, false); // Semidiurnal (3,2)
                SetRange(Swg, 113, 116, false); // Semidiurnal (4,2)
                SetRange(Swg, 117, 120, false); // Semidiurnal (5,2)
                SetRange(Swg, 121, 124, false); // Semidiurnal (6,2)
                SetRange(Swg, 153, 156, false); // Terdiurnal (3,3)
                SetRange(Swg, 157, 160, false); // Terdiurnal (4,3)
                SetRange(Swg, 161, 164, false); // Terdiurnal (5,3)
                SetRange(Swg, 165, 168, false); // Terdiurnal (6,3)
                SetRange(Swg, 197, 200, false); // SPW1 (1,1)
                SetRange(Swg, 205, 208, false); // SPW1 (3,1)
                SetRange(Swg, 213, 216, false); // SPW1 (5,1)
                SetRange(Swg, 259, 262, false); // SPW2 (3,2)
                SetRange(Swg, 267, 270, false); // SPW2 (5,2)
                SetRange(Swg, 394, 397, false); // Geomag (Daily mode terms)
                SetRange(Swg, 407, 410, false); // Mixed UT/Lon/Geomag (Daily mode terms)
                SetRange(Swg, 422, 425, false); // Geomag (Storm-time mode terms)
                SetRange(Swg, 435, 438, false); // Mixed UT/Lon/Geomag (Storm-time mode terms)
                Swg[446] = false; // UT/Lon
            }
            if (Swc[5] == 0) // Asymmetric semiannual
            {
                SetRange(Swg, 221, 224, false); // SPW1 (1,1)
                SetRange(Swg, 229, 232, false); // SPW1 (3,1)
                SetRange(Swg, 237, 240, false); // SPW1 (5,1)
                SetRange(Swg, 279, 282, false); // SPW2 (3,2)
                SetRange(Swg, 287, 290, false); // SPW2 (5,2)
            }
            if (Swc[6] == 0) // Diurnal
            {
                SetRange(Swg, 398, 401, false); // Geomag (Daily mode terms)
                SetRange(Swg, 426, 429, false); // Geomag (Storm-time mode terms)
            }
            if (Swc[10] == 0) // Longitude
            {
                SetRange(Swg, 402, 410, false); // Mixed UT/Lon/Geomag (Daily mode terms)
                SetRange(Swg, 430, 438, false); // Mixed UT/Lon/Geomag (Storm-time mode terms)
                SetRange(Swg, 452, 453, false); // UT/Lon
            }
            if (Swc[11] == 0) // UT/Lon
            {
                SetRange(Swg, 411, 414, false); // Mixed UT/Lon/Geomag (Daily mode terms)
                SetRange(Swg, 439, 440, false); // Mixed UT/Lon/Geomag (Storm-time mode terms)
            }
        }

        // Helper methods for setting ranges of boolean arrays
        private static void SetRange(bool[] array, int start, int end, bool value)
        {
            for (int i = start; i <= end; i++)
            {
                array[i] = value;
            }
        }

        private static void SetIndices(bool[] array, int[] indices, bool value)
        {
            foreach (int i in indices)
            {
                array[i] = value;
            }
        }

        // ==================================================================================================
        // TRETRV: Legacy routine for retrieving switch settings
        // ==================================================================================================
        /// <summary>
        /// Retrieve legacy switch settings
        /// </summary>
        public static float[] TRetrv()
        {
            return (float[])Sav.Clone();
        }
    }
}