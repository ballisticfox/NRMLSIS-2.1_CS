// #######################################################################
// MSIS® (NRL-SOF-014-1) SOFTWARE
// NRLMSIS® empirical atmospheric model software. Use is governed by the
// Open Source Academic Research License Agreement contained in the file
// nrlmsis2.1_license.txt, which is part of this software package. BY
// USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
// CONDITIONS OF THE LICENSE.
// #######################################################################
using NRLMSIS.Infrastructure;

namespace NRLMSIS.Calculators
{
    /// <summary>
    /// Temperature calculation wrapper.
    /// Delegates to existing MSISCalculator implementation.
    /// </summary>
    internal static class TemperatureCalculator
    {
        /// <summary>
        /// Calculates temperature at a given geopotential height.
        /// </summary>
        public static double Calculate(double geopotentialHeight, Context context)
        {
            return MSISCalculator.CalculateTemperature(geopotentialHeight, context);
        }
    }
}
