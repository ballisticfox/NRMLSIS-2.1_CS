// #######################################################################
// MSIS® (NRL-SOF-014-1) SOFTWARE
// NRLMSIS® empirical atmospheric model software. Use is governed by the
// Open Source Academic Research License Agreement contained in the file
// nrlmsis2.1_license.txt, which is part of this software package. BY
// USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
// CONDITIONS OF THE LICENSE.
// #######################################################################

using System;

namespace NRLMSIS.Core
{
    /// <summary>
    /// Represents the spatio-temporal location for an atmospheric calculation.
    /// This is an immutable value type that groups related location parameters.
    /// </summary>
    /// <remarks>
    /// This struct is designed to be lightweight and efficient:
    /// - Stack-allocated (no heap allocation)
    /// - Immutable after construction (thread-safe)
    /// - Value semantics (can be used as dictionary key)
    /// 
    /// For real-time physics applications, creating these structs has negligible overhead
    /// (~40 bytes on stack, no GC pressure).
    /// </remarks>
    public readonly struct AtmosphericLocation : IEquatable<AtmosphericLocation>
    {
        /// <summary>
        /// Gets the day of year (1.0 to 366.0, fractional days allowed).
        /// Example: 172.5 represents noon on the 172nd day of the year.
        /// </summary>
        public double DayOfYear { get; init; }

        /// <summary>
        /// Gets the universal time in seconds (0 to 86400).
        /// </summary>
        public double UniversalTimeSeconds { get; init; }

        /// <summary>
        /// Gets the geodetic latitude in degrees (-90 to 90).
        /// Positive values indicate North, negative values indicate South.
        /// </summary>
        public double LatitudeDegrees { get; init; }

        /// <summary>
        /// Gets the geodetic longitude in degrees (-180 to 360).
        /// Both -180:180 and 0:360 conventions are supported.
        /// </summary>
        public double LongitudeDegrees { get; init; }

        /// <summary>
        /// Gets the altitude in kilometers.
        /// The interpretation (geodetic vs geopotential) depends on configuration.
        /// </summary>
        public double AltitudeKm { get; init; }

        /// <summary>
        /// Creates a new atmospheric location.
        /// </summary>
        /// <param name="dayOfYear">Day of year (1-366)</param>
        /// <param name="universalTimeSeconds">Universal time in seconds (0-86400)</param>
        /// <param name="latitudeDegrees">Latitude in degrees (-90 to 90)</param>
        /// <param name="longitudeDegrees">Longitude in degrees (-180 to 360)</param>
        /// <param name="altitudeKm">Altitude in kilometers</param>
        public AtmosphericLocation(
            double dayOfYear,
            double universalTimeSeconds,
            double latitudeDegrees,
            double longitudeDegrees,
            double altitudeKm)
        {
            DayOfYear = dayOfYear;
            UniversalTimeSeconds = universalTimeSeconds;
            LatitudeDegrees = latitudeDegrees;
            LongitudeDegrees = longitudeDegrees;
            AltitudeKm = altitudeKm;
        }

        /// <summary>
        /// Determines whether this location equals another location.
        /// Two locations are equal if all their components are equal.
        /// </summary>
        public bool Equals(AtmosphericLocation other) =>
            DayOfYear == other.DayOfYear &&
            UniversalTimeSeconds == other.UniversalTimeSeconds &&
            LatitudeDegrees == other.LatitudeDegrees &&
            LongitudeDegrees == other.LongitudeDegrees &&
            AltitudeKm == other.AltitudeKm;

        /// <summary>
        /// Determines whether this location equals another location (excluding altitude).
        /// Useful for checking if basis functions need recalculation.
        /// Basis functions depend on day, time, latitude, and longitude, but not altitude.
        /// </summary>
        public bool EqualsExcludingAltitude(AtmosphericLocation other) =>
            DayOfYear == other.DayOfYear &&
            UniversalTimeSeconds == other.UniversalTimeSeconds &&
            LatitudeDegrees == other.LatitudeDegrees &&
            LongitudeDegrees == other.LongitudeDegrees;

        /// <summary>
        /// Determines whether this location equals an object.
        /// </summary>
        public override bool Equals(object obj) =>
            obj is AtmosphericLocation location && Equals(location);

        /// <summary>
        /// Gets the hash code for this location.
        /// </summary>
        public override int GetHashCode() =>
            HashCode.Combine(DayOfYear, UniversalTimeSeconds, LatitudeDegrees, 
                           LongitudeDegrees, AltitudeKm);

        /// <summary>
        /// Determines whether two locations are equal.
        /// </summary>
        public static bool operator ==(AtmosphericLocation left, AtmosphericLocation right) =>
            left.Equals(right);

        /// <summary>
        /// Determines whether two locations are not equal.
        /// </summary>
        public static bool operator !=(AtmosphericLocation left, AtmosphericLocation right) =>
            !left.Equals(right);

        /// <summary>
        /// Returns a string representation of this location.
        /// </summary>
        public override string ToString() =>
            $"Location(Day={DayOfYear:F1}, UT={UniversalTimeSeconds:F0}s, " +
            $"Lat={LatitudeDegrees:F2}°, Lon={LongitudeDegrees:F2}°, Alt={AltitudeKm:F1}km)";
    }
}
