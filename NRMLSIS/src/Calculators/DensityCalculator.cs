// #######################################################################
// MSIS® (NRL-SOF-014-1) SOFTWARE
// NRLMSIS® empirical atmospheric model software. Use is governed by the
// Open Source Academic Research License Agreement contained in the file
// nrlmsis2.1_license.txt, which is part of this software package. BY
// USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
// CONDITIONS OF THE LICENSE.
// #######################################################################

using NRLMSIS.Core;
using NRLMSIS.Infrastructure;

namespace NRLMSIS.Calculators
{
    /// <summary>
    /// Density calculation wrapper.
    /// Delegates to existing MSISCalculator implementation.
    /// </summary>
    internal static class DensityCalculator
    {
        /// <summary>
        /// Calculates densities for all enabled species.
        /// </summary>
        public static SpeciesDensities CalculateAll(
            double geopotentialHeight,
            double temperature,
            Context context)
        {
            var densities = new double[10];
            MSISCalculator.CalculateDensities(geopotentialHeight, temperature, context, densities);
            return SpeciesDensities.FromArray(densities);
        }
    }
}
