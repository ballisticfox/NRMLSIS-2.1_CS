// #######################################################################
// MSIS® (NRL-SOF-014-1) SOFTWARE
// NRLMSIS® empirical atmospheric model software. Use is governed by the
// Open Source Academic Research License Agreement contained in the file
// nrlmsis2.1_license.txt, which is part of this software package. BY
// USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
// CONDITIONS OF THE LICENSE.
// #######################################################################

using System;
using System.Linq;

namespace NRLMSIS.Core
{
    /// <summary>
    /// Solar and geomagnetic activity indices affecting the atmosphere.
    /// This is an immutable value type that groups related space weather parameters.
    /// </summary>
    /// <remarks>
    /// Space weather conditions drive thermospheric density variations:
    /// - F10.7 solar flux represents solar EUV/UV radiation
    /// - Ap indices represent geomagnetic activity from solar wind
    /// 
    /// This struct is stack-allocated (~64 bytes) with no GC impact.
    /// The ApIndices array is defensively copied to ensure immutability.
    /// 
    /// For real-time physics: Cache this struct and update only when space weather changes
    /// (typically once per day or less).
    /// </remarks>
    public readonly struct SpaceWeather : IEquatable<SpaceWeather>
    {
        /// <summary>
        /// Gets the 81-day average F10.7 solar flux (solar flux units).
        /// This represents solar EUV activity averaged over ~3 solar rotations.
        /// Typical range: 70-250 SFU (1 SFU = 10^-22 W/(m^2 Hz)).
        /// </summary>
        public double F107Average { get; init; }

        /// <summary>
        /// Gets the daily F10.7 solar flux for the previous day (solar flux units).
        /// This captures day-to-day solar activity variations.
        /// </summary>
        public double F107Daily { get; init; }

        /// <summary>
        /// Gets the geomagnetic activity index history array (7 elements).
        /// [0] Daily Ap
        /// [1] 3-hour ap index for current time
        /// [2] 3-hour ap index for 3 hours before current time
        /// [3] 3-hour ap index for 6 hours before current time
        /// [4] 3-hour ap index for 9 hours before current time
        /// [5] Average of eight 3-hour ap indices from 12-33 hours prior
        /// [6] Average of eight 3-hour ap indices from 36-57 hours prior
        /// </summary>
        public double[] ApIndices { get; init; }

        /// <summary>
        /// Creates new space weather conditions.
        /// </summary>
        /// <param name="f107Average">81-day average F10.7</param>
        /// <param name="f107Daily">Daily F10.7 for previous day</param>
        /// <param name="apIndices">Ap geomagnetic activity history (7 elements)</param>
        /// <exception cref="ArgumentNullException">Thrown if apIndices is null</exception>
        /// <exception cref="ArgumentException">Thrown if apIndices doesn't have 7 elements</exception>
        public SpaceWeather(double f107Average, double f107Daily, double[] apIndices)
        {
            if (apIndices == null)
                throw new ArgumentNullException(nameof(apIndices));
            if (apIndices.Length != 7)
                throw new ArgumentException("ApIndices must have exactly 7 elements", nameof(apIndices));

            F107Average = f107Average;
            F107Daily = f107Daily;
            // Create defensive copy to ensure immutability
            ApIndices = (double[])apIndices.Clone();
        }

        /// <summary>
        /// Determines whether this space weather equals another.
        /// </summary>
        public bool Equals(SpaceWeather other)
        {
            if (F107Average != other.F107Average || F107Daily != other.F107Daily)
                return false;

            if (ApIndices == null && other.ApIndices == null)
                return true;
            if (ApIndices == null || other.ApIndices == null)
                return false;

            return ApIndices.SequenceEqual(other.ApIndices);
        }

        /// <summary>
        /// Determines whether this space weather equals an object.
        /// </summary>
        public override bool Equals(object obj) =>
            obj is SpaceWeather weather && Equals(weather);

        /// <summary>
        /// Gets the hash code for this space weather.
        /// </summary>
        public override int GetHashCode()
        {
            var hash = new HashCode();
            hash.Add(F107Average);
            hash.Add(F107Daily);
            if (ApIndices != null)
            {
                foreach (var ap in ApIndices)
                    hash.Add(ap);
            }
            return hash.ToHashCode();
        }

        /// <summary>
        /// Determines whether two space weather instances are equal.
        /// </summary>
        public static bool operator ==(SpaceWeather left, SpaceWeather right) =>
            left.Equals(right);

        /// <summary>
        /// Determines whether two space weather instances are not equal.
        /// </summary>
        public static bool operator !=(SpaceWeather left, SpaceWeather right) =>
            !left.Equals(right);

        /// <summary>
        /// Returns a string representation of this space weather.
        /// </summary>
        public override string ToString() =>
            $"SpaceWeather(F10.7avg={F107Average:F1}, F10.7daily={F107Daily:F1}, " +
            $"Ap={string.Join(",", ApIndices?.Select(a => a.ToString("F0")) ?? new[] { "null" })})";
    }
}
