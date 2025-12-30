// #######################################################################
// MSIS® (NRL-SOF-014-1) SOFTWARE
// NRLMSIS® empirical atmospheric model software. Use is governed by the
// Open Source Academic Research License Agreement contained in the file
// nrlmsis2.1_license.txt, which is part of this software package. BY
// USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
// CONDITIONS OF THE LICENSE.
// #######################################################################

using System.Collections.Generic;

namespace NRLMSIS.Core
{
    /// <summary>
    /// Interface for atmospheric model calculations.
    /// Implementations of this interface can calculate atmospheric parameters
    /// at specified locations and times given space weather conditions.
    /// </summary>
    /// <remarks>
    /// This interface enables:
    /// - Dependency injection in applications
    /// - Mocking for unit tests
    /// - Alternative implementations (cached, simplified, etc.)
    /// - Decorator pattern (logging, validation, etc.)
    /// 
    /// Thread Safety:
    /// - Implementations should be thread-safe for Calculate() calls
    /// - Each thread should use its own Context instance
    /// - Configuration should be immutable
    /// 
    /// Performance:
    /// - For single calculations, use Calculate(location, weather)
    /// - For batch calculations, use CalculateBatch() (reuses cached values)
    /// - For manual control, use Calculate(location, weather, context)
    /// </remarks>
    public interface IAtmosphericModel
    {
        /// <summary>
        /// Calculate atmospheric parameters at a given location and time.
        /// A temporary context is created for this calculation.
        /// </summary>
        /// <param name="location">The spatio-temporal location for the calculation</param>
        /// <param name="weather">Solar and geomagnetic activity conditions</param>
        /// <returns>Complete atmospheric state including temperature and densities</returns>
        /// <exception cref="ArgumentException">Thrown when input parameters are invalid</exception>
        AtmosphericState Calculate(AtmosphericLocation location, SpaceWeather weather);

        /// <summary>
        /// Calculate atmospheric parameters at a given location and time using a provided context.
        /// This overload allows reuse of cached intermediate results for better performance.
        /// </summary>
        /// <param name="location">The spatio-temporal location for the calculation</param>
        /// <param name="weather">Solar and geomagnetic activity conditions</param>
        /// <param name="context">Computation context for caching (reuse across calls for best performance)</param>
        /// <returns>Complete atmospheric state including temperature and densities</returns>
        /// <exception cref="ArgumentException">Thrown when input parameters are invalid</exception>
        /// <remarks>
        /// Performance tip: Reuse the same context across multiple calculations at different altitudes
        /// but the same time/location to avoid recalculating basis functions (3x faster).
        /// </remarks>
        AtmosphericState Calculate(AtmosphericLocation location, SpaceWeather weather, NRLMSIS.Context context);

        /// <summary>
        /// Calculate atmospheric parameters at multiple locations.
        /// This method optimizes calculations by reusing cached results when possible.
        /// </summary>
        /// <param name="locations">Sequence of locations to calculate</param>
        /// <param name="weather">Space weather conditions (same for all locations)</param>
        /// <returns>Sequence of atmospheric states corresponding to each location</returns>
        /// <remarks>
        /// Performance: This method is significantly faster than individual Calculate() calls
        /// when locations share the same time/position (e.g., altitude profile).
        /// Basis functions are calculated once and reused for all locations.
        /// </remarks>
        IEnumerable<AtmosphericState> CalculateBatch(
            IEnumerable<AtmosphericLocation> locations,
            SpaceWeather weather);
    }
}
