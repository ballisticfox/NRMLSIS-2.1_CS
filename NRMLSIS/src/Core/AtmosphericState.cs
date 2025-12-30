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
    /// Complete atmospheric state at a specific location and time.
    /// This is the output from MSIS calculations.
    /// </summary>
    /// <remarks>
    /// This class contains the complete result of an MSIS calculation, including:
    /// - Input parameters (location, weather) for context
    /// - Temperature at the specified altitude
    /// - Exospheric temperature (limiting value at high altitude)
    /// - Densities for all atmospheric species
    /// 
    /// Memory footprint: ~180 bytes (includes SpeciesDensities object)
    /// 
    /// For performance-critical code where allocation matters:
    /// - Reuse instances if possible
    /// - Or use the legacy API with pre-allocated arrays
    /// - Or wait for Phase 3 zero-allocation API
    /// </remarks>
    public sealed class AtmosphericState
    {
        /// <summary>
        /// Gets the location where this state was calculated.
        /// </summary>
        public AtmosphericLocation Location { get; init; }

        /// <summary>
        /// Gets the space weather conditions used for this calculation.
        /// </summary>
        public SpaceWeather Weather { get; init; }

        /// <summary>
        /// Gets the temperature at the specified altitude (Kelvin).
        /// </summary>
        public double Temperature { get; init; }

        /// <summary>
        /// Gets the exospheric temperature (Kelvin).
        /// This is the limiting temperature at very high altitudes (>500 km).
        /// It represents the temperature the atmosphere asymptotically approaches.
        /// </summary>
        public double ExosphericTemperature { get; init; }

        /// <summary>
        /// Gets the atmospheric densities for each species.
        /// </summary>
        public SpeciesDensities Densities { get; init; }

        /// <summary>
        /// Creates a new atmospheric state.
        /// </summary>
        public AtmosphericState()
        {
            Densities = new SpeciesDensities();
        }

        /// <summary>
        /// Returns a string representation of this atmospheric state.
        /// </summary>
        public override string ToString() =>
            $"AtmosphericState(T={Temperature:F1}K, Tex={ExosphericTemperature:F1}K, " +
            $"ρ={Densities.TotalMassDensity:E3} kg/m³)";
    }
}
