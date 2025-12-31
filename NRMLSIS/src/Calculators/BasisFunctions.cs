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
// MSIS_Gfn: Contains subroutines to calculate global (horizontal and time-dependent) model
//           basis functions
// **************************************************************************************************

using System;
using NRLMSIS.Infrastructure;

namespace NRLMSIS.Calculators
{
    /// <summary>
    /// Global basis functions for NRLMSIS 2.1.
    /// Computes horizontal and time-dependent basis functions using:
    /// - Associated Legendre polynomials for latitude dependence
    /// - Fourier harmonics for temporal and longitudinal variations
    /// - Solar flux modulation terms
    /// - Geomagnetic activity effects
    /// </summary>
    public static class BasisFunctions
    {
        // Module-level cache for computed values
        private static double[,] plg = new double[Constants.MaxN + 1, Constants.MaxN + 1]; // Associated Legendre polynomials
        private static double[] cdoy = new double[2]; // Cosine of day-of-year harmonics [annual, semiannual]
        private static double[] sdoy = new double[2]; // Sine of day-of-year harmonics [annual, semiannual]
        private static double[] clst = new double[3]; // Cosine of local solar time harmonics [diurnal, semidiurnal, terdiurnal]
        private static double[] slst = new double[3]; // Sine of local solar time harmonics [diurnal, semidiurnal, terdiurnal]
        private static double[] clon = new double[2]; // Cosine of longitude harmonics [wave-1, wave-2]
        private static double[] slon = new double[2]; // Sine of longitude harmonics [wave-1, wave-2]

        // Reference values for solar flux
        private static double sfluxAvgRef = 150.0; // Reference F10.7 value (NRLMSISE-00 standard)
        private static double sfluxAvgQuadCutoff = 150.0; // Cutoff for truncated quadratic F10.7a function

        // Cache validation - tracks last computed values
        private static double lastLat = -999.9;
        private static double lastDoy = -999.9;
        private static double lastLst = -999.9;
        private static double lastLon = -999.9;

        // ==================================================================================================
        // GLOBE: Calculate horizontal and time-dependent basis functions
        //        (Same purpose as NRLMSISE-00 "GLOBE7" subroutine)
        // ==================================================================================================
        /// <summary>
        /// Calculate horizontal and time-dependent basis functions.
        /// Computes all spatial and temporal variations used in the MSIS model.
        /// </summary>
        /// <param name="doy">Day of year (1-366)</param>
        /// <param name="utsec">Universal time in seconds (0-86400)</param>
        /// <param name="lat">Geodetic latitude in degrees (-90 to 90)</param>
        /// <param name="lon">Geodetic longitude in degrees (-180 to 360)</param>
        /// <param name="sfluxavg">81-day average F10.7 solar flux</param>
        /// <param name="sflux">Daily F10.7 solar flux</param>
        /// <param name="ap">Ap geomagnetic activity index history array [0:6] (7 elements)</param>
        /// <param name="bf">Output: basis function terms [0:maxnbf-1]</param>
        public static void Globe(double doy, double utsec, double lat, double lon,
                                double sfluxavg, double sflux, double[] ap, out double[] bf)
        {
            bf = new double[Constants.MaxNbf];

            // Calculate local solar time
            double lst = (utsec / 3600.0 + lon / 15.0) % 24.0;
            if (lst < 0) lst += 24.0;

            // Update cached trigonometric values (only if inputs changed)
            UpdateLegendrePolynomials(lat);
            UpdateDayOfYearHarmonics(doy);
            UpdateLocalTimeHarmonics(lst);
            UpdateLongitudeHarmonics(lon);

            //---------------------------------------------
            // LINEAR TERMS
            //---------------------------------------------
            int basisIndex = 0;

            basisIndex = BuildTimeIndependentTerms(bf, basisIndex);
            basisIndex = BuildIntraAnnualTerms(bf, basisIndex);
            basisIndex = BuildMigratingTideTerms(bf, basisIndex);
            basisIndex = BuildStationaryPlanetaryWaveTerms(bf, basisIndex);
            basisIndex = BuildSolarFluxTerms(bf, basisIndex, sfluxavg, sflux);
            basisIndex = BuildExtraLinearTerms(bf, basisIndex, sfluxavg, sflux, doy, lst, lat, lon);

            //---------------------------------------------
            // NONLINEAR TERMS
            //---------------------------------------------
            basisIndex = Constants.CNonLin;

            basisIndex = BuildSolarFluxModulationTerms(bf, basisIndex, sfluxavg, sflux);
            basisIndex = BuildGeomagneticTerms(bf, basisIndex, ap, doy, lst, lat, lon, utsec);
            BuildUniversalTimeTerms(bf, Constants.CUt, utsec, doy, sfluxavg, lon);

            //---------------------------------------------
            // Apply Switches
            //---------------------------------------------
            ApplySwitches(bf);
        }

        // ==================================================================================================
        // CACHED VALUE UPDATES - Only recompute when inputs change
        // ==================================================================================================

        /// <summary>
        /// Updates Associated Legendre polynomials if latitude has changed.
        /// Computes P_n^m(cos(colatitude)) for spherical harmonic expansions.
        /// </summary>
        private static void UpdateLegendrePolynomials(double lat)
        {
            if (lat == lastLat)
                return;

            // Note: Legendre polynomials defined in colatitude, so sin(lat) = cos(colat)
            double cosColatitude = Math.Sin(lat * Constants.Deg2Rad);  // clat in original
            double sinColatitude = Math.Cos(lat * Constants.Deg2Rad);  // slat in original
            double cosColat2 = cosColatitude * cosColatitude;
            double cosColat4 = cosColat2 * cosColat2;
            double sinColat2 = sinColatitude * sinColatitude;

            // P_n^0 - Zonal harmonics (m=0)
            plg[0, 0] = 1.0;
            plg[1, 0] = cosColatitude;
            plg[2, 0] = 0.5 * (3.0 * cosColat2 - 1.0);
            plg[3, 0] = 0.5 * (5.0 * cosColatitude * cosColat2 - 3.0 * cosColatitude);
            plg[4, 0] = (35.0 * cosColat4 - 30.0 * cosColat2 + 3.0) / 8.0;
            plg[5, 0] = (63.0 * cosColat2 * cosColat2 * cosColatitude - 70.0 * cosColat2 * cosColatitude + 15.0 * cosColatitude) / 8.0;
            plg[6, 0] = (11.0 * cosColatitude * plg[5, 0] - 5.0 * plg[4, 0]) / 6.0;

            // P_n^1 - Sectoral harmonics (m=1)
            plg[1, 1] = sinColatitude;
            plg[2, 1] = 3.0 * cosColatitude * sinColatitude;
            plg[3, 1] = 1.5 * (5.0 * cosColat2 - 1.0) * sinColatitude;
            plg[4, 1] = 2.5 * (7.0 * cosColat2 * cosColatitude - 3.0 * cosColatitude) * sinColatitude;
            plg[5, 1] = 1.875 * (21.0 * cosColat4 - 14.0 * cosColat2 + 1.0) * sinColatitude;
            plg[6, 1] = (11.0 * cosColatitude * plg[5, 1] - 6.0 * plg[4, 1]) / 5.0;

            // P_n^2 - Tesseral harmonics (m=2)
            plg[2, 2] = 3.0 * sinColat2;
            plg[3, 2] = 15.0 * sinColat2 * cosColatitude;
            plg[4, 2] = 7.5 * (7.0 * cosColat2 - 1.0) * sinColat2;
            plg[5, 2] = 3.0 * cosColatitude * plg[4, 2] - 2.0 * plg[3, 2];
            plg[6, 2] = (11.0 * cosColatitude * plg[5, 2] - 7.0 * plg[4, 2]) / 4.0;

            // P_n^3 - Higher order tesseral harmonics (m=3)
            plg[3, 3] = 15.0 * sinColat2 * sinColatitude;
            plg[4, 3] = 105.0 * sinColat2 * sinColatitude * cosColatitude;
            plg[5, 3] = (9.0 * cosColatitude * plg[4, 3] - 7.0 * plg[3, 3]) / 2.0;
            plg[6, 3] = (11.0 * cosColatitude * plg[5, 3] - 8.0 * plg[4, 3]) / 3.0;

            lastLat = lat;
        }

        /// <summary>
        /// Updates day-of-year Fourier harmonics if day has changed.
        /// Computes annual and semiannual variations.
        /// </summary>
        private static void UpdateDayOfYearHarmonics(double doy)
        {
            if (doy == lastDoy)
                return;

            cdoy[0] = Math.Cos(Constants.Doy2Rad * doy);        // Annual (365.25 day period)
            sdoy[0] = Math.Sin(Constants.Doy2Rad * doy);
            cdoy[1] = Math.Cos(Constants.Doy2Rad * doy * 2.0);  // Semiannual (182.625 day period)
            sdoy[1] = Math.Sin(Constants.Doy2Rad * doy * 2.0);

            lastDoy = doy;
        }

        /// <summary>
        /// Updates local solar time Fourier harmonics if LST has changed.
        /// Computes diurnal, semidiurnal, and terdiurnal variations.
        /// </summary>
        private static void UpdateLocalTimeHarmonics(double lst)
        {
            if (lst == lastLst)
                return;

            clst[0] = Math.Cos(Constants.Lst2Rad * lst);        // Diurnal (24 hour period)
            slst[0] = Math.Sin(Constants.Lst2Rad * lst);
            clst[1] = Math.Cos(Constants.Lst2Rad * lst * 2.0);  // Semidiurnal (12 hour period)
            slst[1] = Math.Sin(Constants.Lst2Rad * lst * 2.0);
            clst[2] = Math.Cos(Constants.Lst2Rad * lst * 3.0);  // Terdiurnal (8 hour period)
            slst[2] = Math.Sin(Constants.Lst2Rad * lst * 3.0);

            lastLst = lst;
        }

        /// <summary>
        /// Updates longitude Fourier harmonics if longitude has changed.
        /// Computes wave-1 and wave-2 planetary wave patterns.
        /// </summary>
        private static void UpdateLongitudeHarmonics(double lon)
        {
            if (lon == lastLon)
                return;

            clon[0] = Math.Cos(Constants.Deg2Rad * lon);        // Wave-1 (360° wavelength)
            slon[0] = Math.Sin(Constants.Deg2Rad * lon);
            clon[1] = Math.Cos(Constants.Deg2Rad * lon * 2.0);  // Wave-2 (180° wavelength)
            slon[1] = Math.Sin(Constants.Deg2Rad * lon * 2.0);

            lastLon = lon;
        }

        // ==================================================================================================
        // BASIS FUNCTION BUILDERS - LINEAR TERMS
        // ==================================================================================================

        /// <summary>
        /// Builds time-independent terms (pure latitude dependence).
        /// Uses zonal Legendre polynomials P_n^0.
        /// </summary>
        private static int BuildTimeIndependentTerms(double[] bf, int startIndex)
        {
            int basisIndex = startIndex;

            if (basisIndex != Constants.CTimeInd)
                throw new InvalidOperationException("Time-independent terms index mismatch");

            for (int latitudeDegree = 0; latitudeDegree <= Constants.AMaxN; latitudeDegree++)
            {
                bf[basisIndex] = plg[latitudeDegree, 0];
                basisIndex++;
            }

            return basisIndex;
        }

        /// <summary>
        /// Builds intra-annual terms (seasonal variations).
        /// Combines zonal Legendre polynomials with annual/semiannual harmonics.
        /// </summary>
        private static int BuildIntraAnnualTerms(double[] bf, int startIndex)
        {
            int basisIndex = startIndex;

            if (basisIndex != Constants.CIntAnn)
                throw new InvalidOperationException("Intra-annual terms index mismatch");

            for (int seasonalHarmonic = 1; seasonalHarmonic <= Constants.AMaxS; seasonalHarmonic++)
            {
                double cosDayOfYear = cdoy[seasonalHarmonic - 1];
                double sinDayOfYear = sdoy[seasonalHarmonic - 1];

                for (int latitudeDegree = 0; latitudeDegree <= Constants.AMaxN; latitudeDegree++)
                {
                    double legendrePolynomial = plg[latitudeDegree, 0];
                    bf[basisIndex] = legendrePolynomial * cosDayOfYear;
                    bf[basisIndex + 1] = legendrePolynomial * sinDayOfYear;
                    basisIndex += 2;
                }
            }

            return basisIndex;
        }

        /// <summary>
        /// Builds migrating tide terms (local time dependence).
        /// Includes tides and their seasonal modulation.
        /// </summary>
        private static int BuildMigratingTideTerms(double[] bf, int startIndex)
        {
            int basisIndex = startIndex;

            if (basisIndex != Constants.CTide)
                throw new InvalidOperationException("Migrating tide terms index mismatch");

            for (int localTimeHarmonic = 1; localTimeHarmonic <= Constants.TMaxL; localTimeHarmonic++)
            {
                double cosLocalSolarTime = clst[localTimeHarmonic - 1];
                double sinLocalSolarTime = slst[localTimeHarmonic - 1];

                // Base tidal components
                for (int latitudeDegree = localTimeHarmonic; latitudeDegree <= Constants.TMaxN; latitudeDegree++)
                {
                    double legendrePolynomial = plg[latitudeDegree, localTimeHarmonic];
                    bf[basisIndex] = legendrePolynomial * cosLocalSolarTime;
                    bf[basisIndex + 1] = legendrePolynomial * sinLocalSolarTime;
                    basisIndex += 2;
                }

                // Seasonal modulation of tides
                for (int seasonalHarmonic = 1; seasonalHarmonic <= Constants.TMaxS; seasonalHarmonic++)
                {
                    double cosDayOfYear = cdoy[seasonalHarmonic - 1];
                    double sinDayOfYear = sdoy[seasonalHarmonic - 1];

                    for (int latitudeDegree = localTimeHarmonic; latitudeDegree <= Constants.TMaxN; latitudeDegree++)
                    {
                        double legendrePolynomial = plg[latitudeDegree, localTimeHarmonic];
                        bf[basisIndex] = legendrePolynomial * cosLocalSolarTime * cosDayOfYear;
                        bf[basisIndex + 1] = legendrePolynomial * sinLocalSolarTime * cosDayOfYear;
                        bf[basisIndex + 2] = legendrePolynomial * cosLocalSolarTime * sinDayOfYear;
                        bf[basisIndex + 3] = legendrePolynomial * sinLocalSolarTime * sinDayOfYear;
                        basisIndex += 4;
                    }
                }
            }

            return basisIndex;
        }

        /// <summary>
        /// Builds stationary planetary wave terms (longitude dependence).
        /// Includes waves and their seasonal modulation.
        /// </summary>
        private static int BuildStationaryPlanetaryWaveTerms(double[] bf, int startIndex)
        {
            int basisIndex = startIndex;

            if (basisIndex != Constants.CSpw)
                throw new InvalidOperationException("Stationary planetary wave terms index mismatch");

            for (int longitudeWaveNumber = 1; longitudeWaveNumber <= Constants.PMaxM; longitudeWaveNumber++)
            {
                double cosLongitude = clon[longitudeWaveNumber - 1];
                double sinLongitude = slon[longitudeWaveNumber - 1];

                // Base planetary wave components
                for (int latitudeDegree = longitudeWaveNumber; latitudeDegree <= Constants.PMaxN; latitudeDegree++)
                {
                    double legendrePolynomial = plg[latitudeDegree, longitudeWaveNumber];
                    bf[basisIndex] = legendrePolynomial * cosLongitude;
                    bf[basisIndex + 1] = legendrePolynomial * sinLongitude;
                    basisIndex += 2;
                }

                // Seasonal modulation of planetary waves
                for (int seasonalHarmonic = 1; seasonalHarmonic <= Constants.PMaxS; seasonalHarmonic++)
                {
                    double cosDayOfYear = cdoy[seasonalHarmonic - 1];
                    double sinDayOfYear = sdoy[seasonalHarmonic - 1];

                    for (int latitudeDegree = longitudeWaveNumber; latitudeDegree <= Constants.PMaxN; latitudeDegree++)
                    {
                        double legendrePolynomial = plg[latitudeDegree, longitudeWaveNumber];
                        bf[basisIndex] = legendrePolynomial * cosLongitude * cosDayOfYear;
                        bf[basisIndex + 1] = legendrePolynomial * sinLongitude * cosDayOfYear;
                        bf[basisIndex + 2] = legendrePolynomial * cosLongitude * sinDayOfYear;
                        bf[basisIndex + 3] = legendrePolynomial * sinLongitude * sinDayOfYear;
                        basisIndex += 4;
                    }
                }
            }

            return basisIndex;
        }

        /// <summary>
        /// Builds solar flux modulation terms.
        /// Linear and quadratic dependence on F10.7 variations.
        /// </summary>
        private static int BuildSolarFluxTerms(double[] bf, int startIndex, double sfluxavg, double sflux)
        {
            int basisIndex = startIndex;

            if (basisIndex != Constants.CSfx)
                throw new InvalidOperationException("Solar flux terms index mismatch");

            double avgFluxDeviation = sfluxavg - sfluxAvgRef;
            double dailyFluxDeviation = sflux - sfluxavg;

            bf[basisIndex] = avgFluxDeviation;                          // Linear in F10.7 average
            bf[basisIndex + 1] = avgFluxDeviation * avgFluxDeviation;   // Quadratic in F10.7 average
            bf[basisIndex + 2] = dailyFluxDeviation;                    // Linear in F10.7 daily variation
            bf[basisIndex + 3] = dailyFluxDeviation * dailyFluxDeviation; // Quadratic in F10.7 daily variation
            bf[basisIndex + 4] = dailyFluxDeviation * avgFluxDeviation;  // Cross term

            basisIndex += Constants.NSfx;

            return basisIndex;
        }

        /// <summary>
        /// Builds additional linear terms including solar zenith angle and
        /// solar flux modulation of various components.
        /// This is the "CExtra" section which includes coupled/modulation terms.
        /// </summary>
        private static int BuildExtraLinearTerms(double[] bf, int startIndex,
                                                 double sfluxavg, double sflux,
                                                 double doy, double lst, double lat, double lon)
        {
            int basisIndex = startIndex;

            if (basisIndex != Constants.CExtra)
                throw new InvalidOperationException("Extra linear terms index mismatch");

            double avgFluxDeviation = sfluxavg - sfluxAvgRef;
            double dailyFluxDeviation = sflux - sfluxavg;
            double solarZenithAngle = SolZen(doy, lst, lat, lon);

            // Solar zenith angle logistic functions
            bf[basisIndex] = -0.5 * Math.Tanh((solarZenithAngle - 98.0) / 6.0);     // For O, H
            bf[basisIndex + 1] = -0.5 * Math.Tanh((solarZenithAngle - 101.5) / 20.0); // For NO

            // Solar flux modulation of SZA terms
            bf[basisIndex + 2] = avgFluxDeviation * bf[basisIndex];
            bf[basisIndex + 3] = avgFluxDeviation * bf[basisIndex + 1];

            // Solar flux modulation of Legendre polynomials
            bf[basisIndex + 4] = avgFluxDeviation * plg[2, 0];  // P(2,0) modulation
            bf[basisIndex + 5] = avgFluxDeviation * plg[4, 0];  // P(4,0) modulation

            // Solar flux modulation of annual/semiannual global variations
            bf[basisIndex + 6] = avgFluxDeviation * plg[0, 0] * cdoy[0];  // Global annual oscillation (AO)
            bf[basisIndex + 7] = avgFluxDeviation * plg[0, 0] * sdoy[0];
            bf[basisIndex + 8] = avgFluxDeviation * plg[0, 0] * cdoy[1];  // Global semiannual oscillation (SAO)
            bf[basisIndex + 9] = avgFluxDeviation * plg[0, 0] * sdoy[1];

            // Truncated quadratic F10.7a term
            if (sfluxavg <= sfluxAvgQuadCutoff)
            {
                bf[basisIndex + 10] = avgFluxDeviation * avgFluxDeviation;
            }
            else
            {
                double cutoffDeviation = sfluxAvgQuadCutoff - sfluxAvgRef;
                bf[basisIndex + 10] = cutoffDeviation * (2.0 * avgFluxDeviation - cutoffDeviation);
            }

            // Legendre polynomial modulation of truncated quadratic term
            bf[basisIndex + 11] = bf[basisIndex + 10] * plg[2, 0];  // P(2,0) modulation
            bf[basisIndex + 12] = bf[basisIndex + 10] * plg[4, 0];  // P(4,0) modulation

            // Daily flux variation modulation
            bf[basisIndex + 13] = dailyFluxDeviation * plg[2, 0];   // P(2,0) modulation of df
            bf[basisIndex + 14] = dailyFluxDeviation * plg[4, 0];   // P(4,0) modulation of df

            // Total: 15 terms in CExtra section
            return basisIndex + 15;
        }

        // ==================================================================================================
        // BASIS FUNCTION BUILDERS - NONLINEAR TERMS
        // ==================================================================================================

        /// <summary>
        /// Builds nonlinear solar flux modulation terms.
        /// This is the "CSfxMod" section (start of nonlinear terms).
        /// </summary>
        private static int BuildSolarFluxModulationTerms(double[] bf, int startIndex,
                                                         double sfluxavg, double sflux)
        {
            int basisIndex = startIndex;

            if (basisIndex != Constants.CSfxMod)
                throw new InvalidOperationException("Solar flux modulation terms index mismatch");

            double avgFluxDeviation = sfluxavg - sfluxAvgRef;
            double dailyFluxDeviation = sflux - sfluxavg;

            bf[basisIndex] = avgFluxDeviation;                              // Linear
            bf[basisIndex + 1] = avgFluxDeviation * avgFluxDeviation;       // Quadratic
            bf[basisIndex + 2] = dailyFluxDeviation;                        // Linear daily
            bf[basisIndex + 3] = dailyFluxDeviation * dailyFluxDeviation;   // Quadratic daily
            bf[basisIndex + 4] = dailyFluxDeviation * avgFluxDeviation;     // Cross term

            basisIndex += Constants.NSfxMod;
            return basisIndex;
        }

        /// <summary>
        /// Builds geomagnetic activity terms for legacy compatibility.
        /// </summary>
        private static int BuildGeomagneticTerms(double[] bf, int startIndex,
                                                double[] ap, double doy, double lst,
                                                double lat, double lon, double utsec)
        {
            int basisIndex = startIndex;

            if (basisIndex != Constants.CMag)
                throw new InvalidOperationException("Geomagnetic terms index mismatch");

            // Ap index history (7 values)
            for (int i = 0; i < 7; i++)
            {
                bf[basisIndex + i] = ap[i] - 4.0;
            }

            // Additional geomagnetic basis functions
            bf[basisIndex + 8] = Constants.Doy2Rad * doy;
            bf[basisIndex + 9] = Constants.Lst2Rad * lst;
            bf[basisIndex + 10] = Constants.Deg2Rad * lon;
            bf[basisIndex + 11] = Constants.Lst2Rad * utsec / 3600.0;
            bf[basisIndex + 12] = Math.Abs(lat);

            basisIndex += 13;

            // Legendre polynomials for m=0,1
            for (int longitudeWaveNumber = 0; longitudeWaveNumber <= 1; longitudeWaveNumber++)
            {
                for (int latitudeDegree = 0; latitudeDegree <= Constants.AMaxN; latitudeDegree++)
                {
                    bf[basisIndex] = plg[latitudeDegree, longitudeWaveNumber];
                    basisIndex++;
                }
            }

            return basisIndex;
        }

        /// <summary>
        /// Builds universal time dependent terms for high-latitude variations.
        /// </summary>
        private static void BuildUniversalTimeTerms(double[] bf, int startIndex,
                                                   double utsec, double doy,
                                                   double sfluxavg, double lon)
        {
            int basisIndex = startIndex;

            if (basisIndex != Constants.CUt)
                throw new InvalidOperationException("UT-dependent terms index mismatch");

            bf[basisIndex] = Constants.Lst2Rad * utsec / 3600.0;
            bf[basisIndex + 1] = Constants.Doy2Rad * doy;
            bf[basisIndex + 2] = sfluxavg - sfluxAvgRef;  // avgFluxDeviation
            bf[basisIndex + 3] = Constants.Deg2Rad * lon;
            bf[basisIndex + 4] = plg[1, 0];
            bf[basisIndex + 5] = plg[3, 0];
            bf[basisIndex + 6] = plg[5, 0];
            bf[basisIndex + 7] = plg[3, 2];
            bf[basisIndex + 8] = plg[5, 2];
        }

        /// <summary>
        /// Applies switches to enable/disable basis function terms.
        /// </summary>
        private static void ApplySwitches(double[] bf)
        {
            for (int i = 0; i <= Constants.Mbf; i++)
            {
                if (!Initialization.Swg[i])
                    bf[i] = 0.0;
            }
        }

        // ==================================================================================================
        // HELPER FUNCTIONS
        // ==================================================================================================

        /// <summary>
        /// Calculate solar zenith angle (adapted from IRI subroutine)
        /// </summary>
        private static double SolZen(double ddd, double lst, double lat, double lon)
        {
            const double Humr = Constants.Pi / 12.0;
            double[] p = { 0.017203534, 0.034407068, 0.051610602, 0.068814136, 0.103221204 };

            double wlon = 360.0 - lon;
            double teqnx = ddd + (lst + wlon / 15.0) / 24.0 + 0.9369;
            teqnx = ddd + 0.9369;

            // Solar declination
            double dec = 23.256 * Math.Sin(p[0] * (teqnx - 82.242)) + 0.381 * Math.Sin(p[1] * (teqnx - 44.855))
                       + 0.167 * Math.Sin(p[2] * (teqnx - 23.355)) - 0.013 * Math.Sin(p[3] * (teqnx + 11.97))
                       + 0.011 * Math.Sin(p[4] * (teqnx - 10.410)) + 0.339137;
            dec = dec * Constants.Deg2Rad;

            // Equation of time
            double tf = teqnx - 0.5;
            double teqt = -7.38 * Math.Sin(p[0] * (tf - 4.0)) - 9.87 * Math.Sin(p[1] * (tf + 9.0))
                        + 0.27 * Math.Sin(p[2] * (tf - 53.0)) - 0.2 * Math.Cos(p[3] * (tf - 17.0));

            double phi = Humr * (lst - 12.0) + teqt * Constants.Deg2Rad / 4.0;
            double rlat = lat * Constants.Deg2Rad;

            // Cosine of solar zenith angle
            double cosx = Math.Sin(rlat) * Math.Sin(dec) + Math.Cos(rlat) * Math.Cos(dec) * Math.Cos(phi);
            if (Math.Abs(cosx) > 1.0) cosx = Math.Sign(cosx) * 1.0;

            return Math.Acos(cosx) / Constants.Deg2Rad;
        }

        // ==================================================================================================
        // PUBLIC INTERFACE METHODS (PRESERVED FROM ORIGINAL FOR COMPATIBILITY)
        // ==================================================================================================

        /// <summary>
        /// Legacy nonlinear modulation of intra-annual, tide, and SPW terms
        /// </summary>
        public static double SFluxMod(int iz, double[] gf, BasisSubset parmset, double dffact)
        {
            double f1, f2, f3, sum;

            // Intra-annual modulation factor
            if (Initialization.Swg[Constants.CSfxMod])
            {
                f1 = parmset.Beta[Constants.CSfxMod, iz] * gf[Constants.CSfxMod]
                   + (parmset.Beta[Constants.CSfx + 2, iz] * gf[Constants.CSfxMod + 2]
                    + parmset.Beta[Constants.CSfx + 3, iz] * gf[Constants.CSfxMod + 3]) * dffact;
            }
            else
            {
                f1 = 0.0;
            }

            // Migrating tide (local time) modulation factor
            if (Initialization.Swg[Constants.CSfxMod + 1])
            {
                f2 = parmset.Beta[Constants.CSfxMod + 1, iz] * gf[Constants.CSfxMod]
                   + (parmset.Beta[Constants.CSfx + 2, iz] * gf[Constants.CSfxMod + 2]
                    + parmset.Beta[Constants.CSfx + 3, iz] * gf[Constants.CSfxMod + 3]) * dffact;
            }
            else
            {
                f2 = 0.0;
            }

            // SPW (longitude) modulation factor
            if (Initialization.Swg[Constants.CSfxMod + 2])
            {
                f3 = parmset.Beta[Constants.CSfxMod + 2, iz] * gf[Constants.CSfxMod];
            }
            else
            {
                f3 = 0.0;
            }

            sum = 0.0;
            for (int j = 0; j <= Constants.Mbf; j++)
            {
                // Apply intra-annual modulation
                if (Initialization.Zsfx[j])
                {
                    sum += parmset.Beta[j, iz] * gf[j] * f1;
                    continue;
                }
                // Apply migrating tide modulation
                if (Initialization.Tsfx[j])
                {
                    sum += parmset.Beta[j, iz] * gf[j] * f2;
                    continue;
                }
                // Apply SPW modulation
                if (Initialization.Psfx[j])
                {
                    sum += parmset.Beta[j, iz] * gf[j] * f3;
                    continue;
                }
            }

            return sum;
        }

        /// <summary>
        /// Legacy nonlinear ap dependence (daily ap mode and ap history mode), including mixed
        /// ap/UT/Longitude terms.
        /// </summary>
        public static double GeoMag(double[] p0, double[] bf, double[,] plg)
        {
            // Return zero if both master switches are off
            if (!(Initialization.Swg[Constants.CMag] || Initialization.Swg[Constants.CMag + 1]))
            {
                return 0.0;
            }

            // Copy parameters
            double[] p = (double[])p0.Clone();
            bool[] swg1 = new bool[Constants.NMag];
            Array.Copy(Initialization.Swg, Constants.CMag, swg1, 0, Constants.NMag);

            double geomag;

            // Calculate function
            if (swg1[0] == swg1[1])
            {
                // Daily Ap mode
                if (p[1] == 0) // If k00s is zero, then cannot compute function
                {
                    return 0.0;
                }

                // Apply switches
                for (int i = 2; i <= 25; i++)
                {
                    if (!swg1[i]) p[i] = 0.0;
                }
                p[8] = p0[8]; // Need doy phase term

                double delA = G0Fn(bf[0], p[0], p[1]);
                geomag = (p[2] * plg[0, 0] + p[3] * plg[2, 0] + p[4] * plg[4, 0]                               // time independent
                    + (p[5] * plg[1, 0] + p[6] * plg[3, 0] + p[7] * plg[5, 0]) * Math.Cos(bf[8] - p[8])        // doy modulation
                    + (p[9] * plg[1, 1] + p[10] * plg[3, 1] + p[11] * plg[5, 1]) * Math.Cos(bf[9] - p[12])     // local time modulation
                    + (1.0 + p[13] * plg[1, 0]) *
                      (p[14] * plg[2, 1] + p[15] * plg[4, 1] + p[16] * plg[6, 1]) * Math.Cos(bf[10] - p[17])   // longitude effect
                    + (p[18] * plg[1, 1] + p[19] * plg[3, 1] + p[20] * plg[5, 1]) * Math.Cos(bf[10] - p[21]) *
                      Math.Cos(bf[8] - p[8])                                                                    // longitude with doy modulation
                    + (p[22] * plg[1, 0] + p[23] * plg[3, 0] + p[24] * plg[5, 0]) * Math.Cos(bf[11] - p[25]))  // universal time
                    * delA;
            }
            else
            {
                // 3-hour ap history mode
                if (p[28] == 0) // If beta00 is zero, then cannot compute function
                {
                    return 0.0;
                }

                // Apply switches
                for (int i = 30; i < Constants.NMag; i++)
                {
                    if (!swg1[i]) p[i] = 0.0;
                }
                p[36] = p0[36]; // Need doy phase term

                double gbeta = p[28] / (1 + p[29] * (45.0 - bf[12]));
                double ex = Math.Exp(-10800.0 * gbeta);
                double sumex = 1 + (1 - Math.Pow(ex, 19.0)) * Math.Pow(ex, 0.5) / (1 - ex);
                double[] G = new double[7]; // G[1:6], index 0 unused
                for (int i = 1; i <= 6; i++)
                {
                    G[i] = G0Fn(bf[i], p[26], p[27]);
                }
                double delA = (G[1]
                            + (G[2] * ex + G[3] * ex * ex + G[4] * Math.Pow(ex, 3.0)
                             + (G[5] * Math.Pow(ex, 4.0) + G[6] * Math.Pow(ex, 12.0)) * (1 - Math.Pow(ex, 8.0)) / (1 - ex))) / sumex;

                geomag = (p[30] * plg[0, 0] + p[31] * plg[2, 0] + p[32] * plg[4, 0]                              // time independent
                    + (p[33] * plg[1, 0] + p[34] * plg[3, 0] + p[35] * plg[5, 0]) * Math.Cos(bf[8] - p[36])      // doy modulation
                    + (p[37] * plg[1, 1] + p[38] * plg[3, 1] + p[39] * plg[5, 1]) * Math.Cos(bf[9] - p[40])      // local time modulation
                    + (1.0 + p[41] * plg[1, 0]) *
                      (p[42] * plg[2, 1] + p[43] * plg[4, 1] + p[44] * plg[6, 1]) * Math.Cos(bf[10] - p[45])     // longitude effect
                    + (p[46] * plg[1, 1] + p[47] * plg[3, 1] + p[48] * plg[5, 1]) * Math.Cos(bf[10] - p[49]) *
                      Math.Cos(bf[8] - p[36])                                                                     // longitude with doy modulation
                    + (p[50] * plg[1, 0] + p[51] * plg[3, 0] + p[52] * plg[5, 0]) * Math.Cos(bf[11] - p[53]))    // universal time
                    * delA;
            }

            return geomag;
        }

        /// <summary>
        /// Helper function for GeoMag
        /// </summary>
        private static double G0Fn(double a, double k00r, double k00s)
        {
            return a + (k00r - 1.0) * (a + (Math.Exp(-a * k00s) - 1.0) / k00s);
        }

        /// <summary>
        /// Legacy nonlinear UT dependence
        /// </summary>
        public static double UtDep(double[] p0, double[] bf)
        {
            // Copy parameters
            double[] p = (double[])p0.Clone();
            bool[] swg1 = new bool[Constants.NUt];
            Array.Copy(Initialization.Swg, Constants.CUt, swg1, 0, Constants.NUt);

            // Apply switches
            for (int i = 3; i < Constants.NUt; i++)
            {
                if (!swg1[i]) p[i] = 0.0;
            }

            // Calculate function
            double utdep = Math.Cos(bf[0] - p[0]) *
                          (1 + p[3] * bf[4] * Math.Cos(bf[1] - p[1])) *
                          (1 + p[4] * bf[2]) * (1 + p[5] * bf[4]) *
                          (p[6] * bf[4] + p[7] * bf[5] + p[8] * bf[6]) +
                          Math.Cos(bf[0] - p[2] + 2 * bf[3]) * (p[9] * bf[7] + p[10] * bf[8]) * (1 + p[11] * bf[2]);

            return utdep;
        }
    }
}