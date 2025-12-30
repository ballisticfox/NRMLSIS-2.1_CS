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
// REFACTORED MSISCALC: Improved version with thread-safety, validation, and readability
//
// Key improvements over original:
// - Thread-safe: Uses Context instead of static state
// - Input validation: Clear error messages for invalid inputs
// - Readable: Helper methods and named constants instead of magic numbers
// - Maintainable: Separated concerns and better documentation
//
// ===========================================================================

using System;

namespace NRLMSIS.Calculators
{
    /// <summary>
    /// Main MSIS calculation entry point with improved thread-safety and validation.
    /// This is a refactored version of the original MsisCalc that uses context objects
    /// for thread-safe caching and provides better error handling.
    /// </summary>
    public static class MSISCalculator
    {
        // ==================================================================================================
        // CALCULATE: Main calculation method with context for thread-safety
        // ==================================================================================================
        /// <summary>
        /// Calculate MSIS atmospheric parameters using a provided context for caching.
        /// This overload is thread-safe and recommended for multi-threaded applications.
        /// </summary>
        /// <param name="day">Day of year (1.0 to 366.0, fractional days allowed)</param>
        /// <param name="utsec">Universal time in seconds (0 to 86400)</param>
        /// <param name="z">Altitude in km (geodetic or geopotential based on initialization)</param>
        /// <param name="lat">Geodetic latitude in degrees (-90 to 90)</param>
        /// <param name="lon">Geodetic longitude in degrees (-180 to 360)</param>
        /// <param name="sfluxavg">81-day average F10.7 solar flux</param>
        /// <param name="sflux">Daily F10.7 solar flux for previous day</param>
        /// <param name="ap">Geomagnetic activity index array [0:6] (7 elements)</param>
        /// <param name="context">Computation context for caching (create one per thread)</param>
        /// <param name="tn">Output: Temperature at altitude (K)</param>
        /// <param name="dn">Output: Density array [0:9] in SI units (see remarks)</param>
        /// <param name="tex">Output: Exospheric temperature (K)</param>
        /// <remarks>
        /// Output density array (dn):
        ///   [0] Total mass density (kg/m³)
        ///   [1] N2 number density (m⁻³)
        ///   [2] O2 number density (m⁻³)
        ///   [3] O number density (m⁻³)
        ///   [4] He number density (m⁻³)
        ///   [5] H number density (m⁻³)
        ///   [6] Ar number density (m⁻³)
        ///   [7] N number density (m⁻³)
        ///   [8] Anomalous O number density (m⁻³)
        ///   [9] NO number density (m⁻³)
        ///
        /// Missing values are returned as 9.999e-38.
        /// Species calculated are determined by Initialization configuration.
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when any input parameter is outside valid range.
        /// </exception>
        /// <exception cref="ArgumentException">
        /// Thrown when ap array has incorrect length.
        /// </exception>
        /// <exception cref="InvalidOperationException">
        /// Thrown when MSIS has not been initialized.
        /// </exception>
        public static void Calculate(
            double day, double utsec, double z, double lat, double lon,
            double sfluxavg, double sflux, double[] ap,
            Context context,
            out double tn, out double[] dn, out double tex)
        {
            // Validate initialization
            InputValidator.ValidateInitialization();

            // Validate all input parameters
            InputValidator.ValidateInputs(day, utsec, z, lat, lon, sfluxavg, sflux, ap);

            // Initialize output
            dn = new double[10];
            tex = 0.0;
            tn = 0.0;

            // Convert altitude to geopotential height if needed
            double geopotentialHeight = Initialization.ZAltFlag
                ? Utilities.Alt2Gph(lat, z)
                : z;

            // Check if we need to recalculate basis functions
            // (they depend on location, time, and space weather)
            bool needsProfileUpdate = context.IsLocationTimeChanged(
                day, utsec, lat, lon, sflux, sfluxavg, ap);

            if (needsProfileUpdate)
            {
                // Calculate horizontal and temporal basis functions
                BasisFunctions.Globe(day, utsec, lat, lon, sfluxavg, sflux, ap,
                            out double[] gf);

                // Copy to context cache
                Array.Copy(gf, context.BasisFunctions, gf.Length);

                // Calculate temperature profile parameters
                TemperatureFunction.TemperatureParameters(context.BasisFunctions, out var tpro);

                // Copy temperature profile to context
                CopyTemperatureProfile(tpro, context.TemperatureProfile);

                // Calculate density profile parameters for each species
                for (int ispec = 2; ispec <= 10; ispec++)
                {
                    if (Initialization.SpecFlag[ispec - 1])
                    {
                        DensityProfile.DensityParameters(ispec, context.BasisFunctions,
                                       context.TemperatureProfile,
                                       out var dpro);
                        CopyDensityProfile(dpro, context.DensityProfiles[ispec]);
                    }
                }

                // Update cache markers
                context.UpdateLocationTimeCache(day, utsec, lat, lon, sflux, sfluxavg, ap);
            }

            // Check if we need to recalculate spline weights
            // (they only depend on altitude)
            if (geopotentialHeight < AltitudeRegimes.ZetaB)
            {
                bool needsSplineUpdate = context.IsAltitudeChanged(geopotentialHeight);

                if (needsSplineUpdate)
                {
                    // Determine maximum spline order needed
                    int maxOrder = geopotentialHeight < AltitudeRegimes.ZetaF ? 5 : 6;

                    // Calculate B-spline weights
                    Utilities.BSpline(
                        geopotentialHeight,
                        Constants.NodesTN,
                        Constants.Nd + 2,
                        maxOrder,
                        Initialization.EtaTN,
                        out double[,] splineWeights,
                        out int splineIndex);

                    // Copy to context cache
                    Array.Copy(splineWeights, context.SplineWeights, splineWeights.Length);
                    context.SplineIndex = splineIndex;

                    // Update altitude cache marker
                    context.UpdateAltitudeCache(geopotentialHeight);
                }
            }

            // Calculate temperature at altitude
            tn = CalculateTemperature(geopotentialHeight, context);

            // Get exospheric temperature
            tex = context.TemperatureProfile.Tex;

            // Calculate densities for all species
            CalculateDensities(geopotentialHeight, tn, context, dn);
        }

        // ==================================================================================================
        // CALCULATE: Overload without context (creates temporary context)
        // ==================================================================================================
        /// <summary>
        /// Calculate MSIS atmospheric parameters without explicit context.
        /// A temporary context is created for this calculation.
        /// For multi-threaded use or repeated calculations, use the overload with Context.
        /// </summary>
        public static void Calculate(
            double day, double utsec, double z, double lat, double lon,
            double sfluxavg, double sflux, double[] ap,
            out double tn, out double[] dn, out double tex)
        {
            using var context = new Context();
            Calculate(day, utsec, z, lat, lon, sfluxavg, sflux, ap,
                     context, out tn, out dn, out tex);
        }

        // ==================================================================================================
        // CALCULATE: Overload without tex output
        // ==================================================================================================
        /// <summary>
        /// Calculate MSIS atmospheric parameters without exospheric temperature output.
        /// </summary>
        public static void Calculate(
            double day, double utsec, double z, double lat, double lon,
            double sfluxavg, double sflux, double[] ap,
            Context context,
            out double tn, out double[] dn)
        {
            Calculate(day, utsec, z, lat, lon, sfluxavg, sflux, ap,
                     context, out tn, out dn, out _);
        }

        // ==================================================================================================
        // Private helper methods
        // ==================================================================================================

        /// <summary>
        /// Calculates temperature at the specified altitude using cached profile parameters.
        /// </summary>
        internal static double CalculateTemperature(double altitude, Context context)
        {
            // Extract spline weights for temperature calculation (order 4 = cubic)
            // Fortran: Sz(-3:0, 4) maps to C# Sz[2:5, 2]
            Span<double> wght = context.WeightBuffer;
            for (int idx = 0; idx < 4; idx++)
            {
                wght[idx] = context.SplineWeights[idx + 2, 2]; // Order 4 -> index 2
            }

            // Calculate temperature at altitude
            return TemperatureFunction.TfnX(altitude, context.SplineIndex, wght.ToArray(),
                               context.TemperatureProfile);
        }

        /// <summary>
        /// Calculates densities for all enabled species at the specified altitude.
        /// </summary>
        internal static void CalculateDensities(
            double altitude,
            double temperature,
            Context context,
            double[] densities)
        {
            var tpro = context.TemperatureProfile;

            // Calculate temperature integration terms
            double delz = altitude - AltitudeRegimes.ZetaB;
            double lndtotz, Vz, Wz;

            if (altitude < AltitudeRegimes.ZetaF)
            {
                // Below fully mixed region (< 70 km)
                CalculateIntegralsLowAltitude(altitude, temperature, context,
                                             out lndtotz, out Vz, out Wz);
            }
            else if (altitude < AltitudeRegimes.ZetaB)
            {
                // Transition region (70-122.5 km)
                CalculateIntegralsMidAltitude(altitude, temperature, context,
                                             out lndtotz, out Vz, out Wz);
            }
            else
            {
                // Upper thermosphere (> 122.5 km) - Bates profile
                CalculateIntegralsHighAltitude(altitude, temperature, delz, tpro,
                                              out lndtotz, out Vz, out Wz);
            }

            // Calculate chemical/dynamical taper factor
            double taperFactor = 0.5 * (1.0 + Math.Tanh(
                AltitudeRegimes.HGamma * (altitude - AltitudeRegimes.ZetaGamma)));

            // Calculate number densities for each species (2-10)
            for (int ispec = 2; ispec <= 10; ispec++)
            {
                if (Initialization.SpecFlag[ispec - 1])
                {
                    densities[ispec - 1] = DensityProfile.DfnX(
                        altitude, temperature, lndtotz, Vz, Wz, taperFactor,
                        tpro, context.DensityProfiles[ispec]);
                }
                else
                {
                    densities[ispec - 1] = Constants.DMissing;
                }
            }

            // Calculate total mass density if enabled
            if (Initialization.SpecFlag[0])
            {
                densities[0] = 0.0;
                for (int idx = 1; idx < 10; idx++)
                {
                    if (densities[idx] != Constants.DMissing)
                    {
                        densities[0] += densities[idx] * Initialization.MassWgt[idx];
                    }
                }
            }
            else
            {
                densities[0] = Constants.DMissing;
            }
        }

        /// <summary>
        /// Calculates integration terms for low altitude (below 70 km).
        /// </summary>
        private static void CalculateIntegralsLowAltitude(
            double altitude,
            double temperature,
            Context context,
            out double lnTotalDensity,
            out double firstIntegral,
            out double secondIntegral)
        {
            var tpro = context.TemperatureProfile;
            int iz = context.SplineIndex;

            // Calculate first integral (V) using 5th order splines
            int i = Math.Max(iz - 4, 0);
            int j = iz < 4 ? -iz : -4;

            firstIntegral = 0.0;
            int szIdx = j + 5; // Map Fortran index to C#
            for (int idx = i; idx <= iz; idx++)
            {
                firstIntegral += tpro.Beta[idx] * context.SplineWeights[szIdx, 3]; // Order 5 -> index 3
                szIdx++;
            }
            firstIntegral += tpro.CVs;

            // Second integral not needed for low altitude
            secondIntegral = 0.0;

            // Calculate log pressure and total number density
            double lnPressure = Constants.LnP0 -
                Constants.MbarG0DivKB * (firstIntegral - tpro.VZeta0);
            lnTotalDensity = lnPressure - Math.Log(Constants.KB * temperature);
        }

        /// <summary>
        /// Calculates integration terms for mid altitude (70-122.5 km).
        /// </summary>
        private static void CalculateIntegralsMidAltitude(
            double altitude,
            double temperature,
            Context context,
            out double lnTotalDensity,
            out double firstIntegral,
            out double secondIntegral)
        {
            var tpro = context.TemperatureProfile;
            int iz = context.SplineIndex;

            // Calculate first integral (V) using 5th order splines
            firstIntegral = 0.0;
            for (int idx = 0; idx <= 4; idx++)
            {
                firstIntegral += tpro.Beta[iz - 4 + idx] * context.SplineWeights[idx + 1, 3];
            }
            firstIntegral += tpro.CVs;

            // Calculate second integral (W) using 6th order splines
            secondIntegral = 0.0;
            for (int idx = 0; idx <= 5; idx++)
            {
                secondIntegral += tpro.Gamma[iz - 5 + idx] * context.SplineWeights[idx, 4];
            }
            double delz = altitude - AltitudeRegimes.ZetaB;
            secondIntegral += tpro.CVs * delz + tpro.CWs;

            // Log total density not used in this region
            lnTotalDensity = 0.0;
        }

        /// <summary>
        /// Calculates integration terms for high altitude (above 122.5 km).
        /// </summary>
        private static void CalculateIntegralsHighAltitude(
            double altitude,
            double temperature,
            double delz,
            TemperatureProfile tpro,
            out double lnTotalDensity,
            out double firstIntegral,
            out double secondIntegral)
        {
            // Bates profile integrals
            firstIntegral = (delz + Math.Log(temperature / tpro.Tex) / tpro.Sigma) / tpro.Tex
                          + tpro.CVb;

            secondIntegral = (0.5 * delz * delz +
                            Utilities.Dilog(tpro.B * Math.Exp(-tpro.Sigma * delz)) / tpro.SigmaSq)
                           / tpro.Tex + tpro.CVb * delz + tpro.CWb;

            // Log total density not used in this region
            lnTotalDensity = 0.0;
        }

        /// <summary>
        /// Deep copies temperature profile parameters.
        /// </summary>
        private static void CopyTemperatureProfile(TemperatureProfile source, TemperatureProfile dest)
        {
            Array.Copy(source.Cf, dest.Cf, source.Cf.Length);
            Array.Copy(source.Beta, dest.Beta, source.Beta.Length);
            Array.Copy(source.Gamma, dest.Gamma, source.Gamma.Length);

            dest.TZetaF = source.TZetaF;
            dest.TZetaA = source.TZetaA;
            dest.DlnTdzA = source.DlnTdzA;
            dest.LnDTotF = source.LnDTotF;
            dest.Tex = source.Tex;
            dest.Tgb0 = source.Tgb0;
            dest.Tb0 = source.Tb0;
            dest.Sigma = source.Sigma;
            dest.SigmaSq = source.SigmaSq;
            dest.B = source.B;
            dest.CVs = source.CVs;
            dest.CVb = source.CVb;
            dest.CWs = source.CWs;
            dest.CWb = source.CWb;
            dest.VZetaF = source.VZetaF;
            dest.VZetaA = source.VZetaA;
            dest.WZetaA = source.WZetaA;
            dest.VZeta0 = source.VZeta0;
        }

        /// <summary>
        /// Deep copies density profile parameters.
        /// </summary>
        private static void CopyDensityProfile(DensityParameters source, DensityParameters dest)
        {
            if (source.Cf != null && dest.Cf != null)
            {
                Array.Copy(source.Cf, dest.Cf, Math.Min(source.Cf.Length, dest.Cf.Length));
            }
            if (source.Mi != null && dest.Mi != null)
            {
                Array.Copy(source.Mi, dest.Mi, source.Mi.Length);
            }
            if (source.ZetaMi != null && dest.ZetaMi != null)
            {
                Array.Copy(source.ZetaMi, dest.ZetaMi, source.ZetaMi.Length);
            }
            if (source.AMi != null && dest.AMi != null)
            {
                Array.Copy(source.AMi, dest.AMi, source.AMi.Length);
            }
            if (source.WMi != null && dest.WMi != null)
            {
                Array.Copy(source.WMi, dest.WMi, source.WMi.Length);
            }
            if (source.XMi != null && dest.XMi != null)
            {
                Array.Copy(source.XMi, dest.XMi, source.XMi.Length);
            }

            dest.LnPhiF = source.LnPhiF;
            dest.LnDRef = source.LnDRef;
            dest.ZetaM = source.ZetaM;
            dest.HML = source.HML;
            dest.HMU = source.HMU;
            dest.C = source.C;
            dest.ZetaC = source.ZetaC;
            dest.HC = source.HC;
            dest.R = source.R;
            dest.ZetaR = source.ZetaR;
            dest.HR = source.HR;
            dest.ZRef = source.ZRef;
            dest.IzRef = source.IzRef;
            dest.TRef = source.TRef;
            dest.ZMin = source.ZMin;
            dest.ZHyd = source.ZHyd;
            dest.ISpec = source.ISpec;
        }
    }
}