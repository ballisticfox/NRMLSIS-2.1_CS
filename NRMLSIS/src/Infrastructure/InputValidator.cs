// #######################################################################
// MSIS® (NRL-SOF-014-1) SOFTWARE
// NRLMSIS® empirical atmospheric model software. Use is governed by the
// Open Source Academic Research License Agreement contained in the file
// nrlmsis2.1_license.txt, which is part of this software package. BY
// USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
// CONDITIONS OF THE LICENSE.
// #######################################################################

using System;

namespace NRLMSIS
{
    /// <summary>
    /// Validates input parameters for MSIS atmospheric model calculations.
    /// Provides clear error messages when inputs are out of valid ranges.
    /// </summary>
    public static class InputValidator
    {
        // Valid ranges for input parameters
        private const double MinDayOfYear = 1.0;
        private const double MaxDayOfYear = 366.0;

        private const double MinUniversalTime = 0.0;
        private const double MaxUniversalTime = 86400.0;

        private const double MinAltitude = -5.0;   // Slightly below sea level
        private const double MaxAltitude = 1000.0; // Upper limit of model validity

        private const double MinLatitude = -90.0;
        private const double MaxLatitude = 90.0;

        private const double MinLongitude = -180.0;
        private const double MaxLongitude = 360.0;  // Allow both -180:180 and 0:360

        private const double MinF107 = 0.0;
        private const double MaxF107 = 500.0;  // Extremely high, but physically possible

        private const int RequiredApLength = 7;

        /// <summary>
        /// Validates all input parameters for a MSIS calculation.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when any parameter is outside its valid range.
        /// </exception>
        /// <exception cref="ArgumentException">
        /// Thrown when array parameters have incorrect length.
        /// </exception>
        public static void ValidateInputs(
            double day,
            double utsec,
            double altitude,
            double latitude,
            double longitude,
            double sfluxavg,
            double sflux,
            double[] ap)
        {
            ValidateDayOfYear(day);
            ValidateUniversalTime(utsec);
            ValidateAltitude(altitude);
            ValidateLatitude(latitude);
            ValidateLongitude(longitude);
            ValidateF107Average(sfluxavg);
            ValidateF107Daily(sflux);
            ValidateApArray(ap);
        }

        /// <summary>
        /// Validates day of year parameter.
        /// </summary>
        /// <param name="day">Day of year (1-366)</param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when day is not in range [1, 366].
        /// </exception>
        public static void ValidateDayOfYear(double day)
        {
            if (double.IsNaN(day))
            {
                throw new ArgumentException(
                    "Day of year cannot be NaN",
                    nameof(day));
            }

            if (double.IsInfinity(day))
            {
                throw new ArgumentException(
                    "Day of year cannot be infinite",
                    nameof(day));
            }

            if (day < MinDayOfYear || day > MaxDayOfYear)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(day),
                    day,
                    $"Day of year must be between {MinDayOfYear} and {MaxDayOfYear}. " +
                    $"For example, January 1 = 1, December 31 = 365 (or 366 for leap years). " +
                    $"Fractional days are allowed (e.g., 172.5 for noon on day 172).");
            }
        }

        /// <summary>
        /// Validates universal time parameter.
        /// </summary>
        /// <param name="utsec">Universal time in seconds (0-86400)</param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when utsec is not in range [0, 86400).
        /// </exception>
        public static void ValidateUniversalTime(double utsec)
        {
            if (double.IsNaN(utsec))
            {
                throw new ArgumentException(
                    "Universal time cannot be NaN",
                    nameof(utsec));
            }

            if (double.IsInfinity(utsec))
            {
                throw new ArgumentException(
                    "Universal time cannot be infinite",
                    nameof(utsec));
            }

            if (utsec < MinUniversalTime || utsec >= MaxUniversalTime)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(utsec),
                    utsec,
                    $"Universal time must be in range [{MinUniversalTime}, {MaxUniversalTime}) seconds. " +
                    $"Note: 86400 seconds = 24 hours. For times >= 86400, consider incrementing the day.");
            }
        }

        /// <summary>
        /// Validates altitude parameter.
        /// </summary>
        /// <param name="altitude">Altitude in kilometers</param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when altitude is not in range [-5, 1000] km.
        /// </exception>
        public static void ValidateAltitude(double altitude)
        {
            if (double.IsNaN(altitude))
            {
                throw new ArgumentException(
                    "Altitude cannot be NaN",
                    nameof(altitude));
            }

            if (double.IsInfinity(altitude))
            {
                throw new ArgumentException(
                    "Altitude cannot be infinite",
                    nameof(altitude));
            }

            if (altitude < MinAltitude || altitude > MaxAltitude)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(altitude),
                    altitude,
                    $"Altitude must be between {MinAltitude} and {MaxAltitude} km. " +
                    $"The MSIS model is valid from the Earth's surface up to ~1000 km. " +
                    $"For altitudes > 1000 km, results may be less accurate.");
            }

            // Warning for altitudes above typical thermosphere
            if (altitude > 600)
            {
                // Could log a warning here if logging is available
                // For now, just allow it to proceed
            }
        }

        /// <summary>
        /// Validates latitude parameter.
        /// </summary>
        /// <param name="latitude">Latitude in degrees (-90 to 90)</param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when latitude is not in range [-90, 90] degrees.
        /// </exception>
        public static void ValidateLatitude(double latitude)
        {
            if (double.IsNaN(latitude))
            {
                throw new ArgumentException(
                    "Latitude cannot be NaN",
                    nameof(latitude));
            }

            if (double.IsInfinity(latitude))
            {
                throw new ArgumentException(
                    "Latitude cannot be infinite",
                    nameof(latitude));
            }

            if (latitude < MinLatitude || latitude > MaxLatitude)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(latitude),
                    latitude,
                    $"Latitude must be between {MinLatitude} and {MaxLatitude} degrees. " +
                    $"Positive values indicate North, negative values indicate South. " +
                    $"For example: 45.0 = 45°N, -30.0 = 30°S");
            }
        }

        /// <summary>
        /// Validates longitude parameter.
        /// </summary>
        /// <param name="longitude">Longitude in degrees</param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when longitude is not in valid range.
        /// </exception>
        public static void ValidateLongitude(double longitude)
        {
            if (double.IsNaN(longitude))
            {
                throw new ArgumentException(
                    "Longitude cannot be NaN",
                    nameof(longitude));
            }

            if (double.IsInfinity(longitude))
            {
                throw new ArgumentException(
                    "Longitude cannot be infinite",
                    nameof(longitude));
            }

            // Accept both -180:180 and 0:360 conventions
            if (longitude < MinLongitude || longitude > MaxLongitude)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(longitude),
                    longitude,
                    $"Longitude must be between {MinLongitude} and {MaxLongitude} degrees. " +
                    $"You can use either -180:180 or 0:360 convention. " +
                    $"Positive values indicate East, negative values indicate West. " +
                    $"For example: -122.0 = 122°W, 238.0 = 122°W (equivalent)");
            }
        }

        /// <summary>
        /// Validates F10.7 average flux parameter.
        /// </summary>
        /// <param name="f107avg">81-day average F10.7 solar flux</param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when F10.7 average is not in valid range.
        /// </exception>
        public static void ValidateF107Average(double f107avg)
        {
            if (double.IsNaN(f107avg))
            {
                throw new ArgumentException(
                    "F10.7 average cannot be NaN",
                    nameof(f107avg));
            }

            if (double.IsInfinity(f107avg))
            {
                throw new ArgumentException(
                    "F10.7 average cannot be infinite",
                    nameof(f107avg));
            }

            if (f107avg < MinF107 || f107avg > MaxF107)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(f107avg),
                    f107avg,
                    $"F10.7 average must be between {MinF107} and {MaxF107} SFU. " +
                    $"Solar flux units (SFU): 1 SFU = 10^-22 W/(m^2 Hz). " +
                    $"Typical range during solar cycle: 70-250 SFU. " +
                    $"Note: Use the observed flux at Earth distance, not adjusted to 1 AU.");
            }
        }

        /// <summary>
        /// Validates F10.7 daily flux parameter.
        /// </summary>
        /// <param name="f107">Daily F10.7 solar flux</param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when F10.7 daily is not in valid range.
        /// </exception>
        public static void ValidateF107Daily(double f107)
        {
            if (double.IsNaN(f107))
            {
                throw new ArgumentException(
                    "F10.7 daily cannot be NaN",
                    nameof(f107));
            }

            if (double.IsInfinity(f107))
            {
                throw new ArgumentException(
                    "F10.7 daily cannot be infinite",
                    nameof(f107));
            }

            if (f107 < MinF107 || f107 > MaxF107)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(f107),
                    f107,
                    $"F10.7 daily must be between {MinF107} and {MaxF107} SFU. " +
                    $"This should be the observed flux for the previous day. " +
                    $"Note: Use the observed flux at Earth distance, not adjusted to 1 AU.");
            }
        }

        /// <summary>
        /// Validates Ap geomagnetic activity index array.
        /// </summary>
        /// <param name="ap">Array of 7 Ap indices</param>
        /// <exception cref="ArgumentNullException">
        /// Thrown when ap array is null.
        /// </exception>
        /// <exception cref="ArgumentException">
        /// Thrown when ap array doesn't have exactly 7 elements.
        /// </exception>
        public static void ValidateApArray(double[] ap)
        {
            if (ap == null)
            {
                throw new ArgumentNullException(
                    nameof(ap),
                    "Ap array cannot be null. It must contain 7 geomagnetic activity indices.");
            }

            if (ap.Length != RequiredApLength)
            {
                throw new ArgumentException(
                    $"Ap array must have exactly {RequiredApLength} elements, but has {ap.Length}. " +
                    $"The array should contain:\n" +
                    $"  [0] Daily Ap\n" +
                    $"  [1] 3-hour ap index for current time\n" +
                    $"  [2] 3-hour ap index for 3 hours before current time\n" +
                    $"  [3] 3-hour ap index for 6 hours before current time\n" +
                    $"  [4] 3-hour ap index for 9 hours before current time\n" +
                    $"  [5] Average of eight 3-hour ap indices from 12-33 hours prior\n" +
                    $"  [6] Average of eight 3-hour ap indices from 36-57 hours prior",
                    nameof(ap));
            }

            // Validate each Ap value
            for (int i = 0; i < ap.Length; i++)
            {
                if (double.IsNaN(ap[i]))
                {
                    throw new ArgumentException(
                        $"Ap[{i}] cannot be NaN",
                        nameof(ap));
                }

                if (double.IsInfinity(ap[i]))
                {
                    throw new ArgumentException(
                        $"Ap[{i}] cannot be infinite",
                        nameof(ap));
                }

                if (ap[i] < 0 || ap[i] > 400)
                {
                    throw new ArgumentOutOfRangeException(
                        nameof(ap),
                        ap[i],
                        $"Ap[{i}] must be between 0 and 400. " +
                        $"Typical range: 0-100, with extreme geomagnetic storms reaching ~400. " +
                        $"Current value: {ap[i]}");
                }
            }
        }

        /// <summary>
        /// Validates that MSIS has been initialized before use.
        /// </summary>
        /// <exception cref="InvalidOperationException">
        /// Thrown when MSIS has not been initialized.
        /// </exception>
        public static void ValidateInitialization()
        {
            if (!Initialization.InitFlag)
            {
                throw new InvalidOperationException(
                    "MSIS model has not been initialized. " +
                    "Call Initialization.Initialize() before performing calculations. " +
                    "Example: Initialization.Initialize(parmPath: \"./\", parmFile: \"msis21.parm\");");
            }
        }
    }
}