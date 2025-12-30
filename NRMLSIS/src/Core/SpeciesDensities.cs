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
    /// Number densities for atmospheric species in SI units (m^-3 for number densities, kg/m^3 for mass density).
    /// Missing values are represented as Constants.DMissing (9.999e-38).
    /// </summary>
    /// <remarks>
    /// This class provides named properties for atmospheric species densities,
    /// replacing the legacy array-based approach with magic indices.
    /// 
    /// Benefits:
    /// - IntelliSense support
    /// - Compile-time safety (can't use wrong index)
    /// - Self-documenting code
    /// - Easy to extend with computed properties
    /// 
    /// Memory footprint: ~100 bytes
    /// 
    /// For performance-critical code:
    /// - These are simple property getters/setters (inline-able)
    /// - No performance penalty vs array access
    /// - Can convert to/from arrays for legacy compatibility
    /// 
    /// Species included:
    /// - Major: N2, O2, O, Ar (always present)
    /// - Minor: He, H, N (altitude dependent)
    /// - Trace: Anomalous O, NO (high altitude)
    /// </remarks>
    public sealed class SpeciesDensities
    {
        /// <summary>Gets or sets the total mass density (kg/m³)</summary>
        public double TotalMassDensity { get; set; }

        /// <summary>Gets or sets the N2 number density (m⁻³)</summary>
        public double N2 { get; set; }

        /// <summary>Gets or sets the O2 number density (m⁻³)</summary>
        public double O2 { get; set; }

        /// <summary>Gets or sets the O number density (m⁻³)</summary>
        public double O { get; set; }

        /// <summary>Gets or sets the He number density (m⁻³)</summary>
        public double He { get; set; }

        /// <summary>Gets or sets the H number density (m⁻³)</summary>
        public double H { get; set; }

        /// <summary>Gets or sets the Ar number density (m⁻³)</summary>
        public double Ar { get; set; }

        /// <summary>Gets or sets the N number density (m⁻³)</summary>
        public double N { get; set; }

        /// <summary>
        /// Gets or sets the anomalous oxygen number density (m⁻³).
        /// This represents a separate oxygen population with different scale height,
        /// primarily important at high altitudes (>200 km).
        /// </summary>
        public double AnomalousO { get; set; }

        /// <summary>Gets or sets the NO number density (m⁻³)</summary>
        public double NO { get; set; }

        /// <summary>
        /// Converts this species densities to the legacy array format.
        /// Array indices: [0]=Mass, [1]=N2, [2]=O2, [3]=O, [4]=He, [5]=H, [6]=Ar, [7]=N, [8]=AnomalousO, [9]=NO
        /// </summary>
        /// <returns>Array with 10 elements in MSIS internal order</returns>
        public double[] ToArray()
        {
            return new double[]
            {
                TotalMassDensity, N2, O2, O, He, H, Ar, N, AnomalousO, NO
            };
        }

        /// <summary>
        /// Creates species densities from a legacy array format.
        /// Array must have 10 elements in MSIS internal order.
        /// </summary>
        /// <param name="array">Array with [0]=Mass, [1]=N2, [2]=O2, [3]=O, [4]=He, [5]=H, [6]=Ar, [7]=N, [8]=AnomalousO, [9]=NO</param>
        /// <returns>SpeciesDensities instance with values from array</returns>
        /// <exception cref="ArgumentNullException">Thrown if array is null</exception>
        /// <exception cref="ArgumentException">Thrown if array doesn't have 10 elements</exception>
        public static SpeciesDensities FromArray(double[] array)
        {
            if (array == null)
                throw new ArgumentNullException(nameof(array));
            if (array.Length != 10)
                throw new ArgumentException("Array must have 10 elements", nameof(array));

            return new SpeciesDensities
            {
                TotalMassDensity = array[0],
                N2 = array[1],
                O2 = array[2],
                O = array[3],
                He = array[4],
                H = array[5],
                Ar = array[6],
                N = array[7],
                AnomalousO = array[8],
                NO = array[9]
            };
        }

        /// <summary>
        /// Checks if a density value is missing (i.e., equals Constants.DMissing).
        /// Missing values occur when a species is disabled in configuration or
        /// doesn't exist at the specified altitude.
        /// </summary>
        /// <param name="density">Density value to check</param>
        /// <returns>True if the value represents missing data, false otherwise</returns>
        public static bool IsMissing(double density) =>
            Math.Abs(density - NRLMSIS.Constants.DMissing) < 1e-40;

        /// <summary>
        /// Returns a string representation of the densities.
        /// Shows key species for brevity.
        /// </summary>
        public override string ToString() =>
            $"Densities(ρ={TotalMassDensity:E3}, N2={N2:E3}, O2={O2:E3}, O={O:E3})";
    }
}
