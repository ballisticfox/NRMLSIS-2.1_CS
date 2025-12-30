// #######################################################################
// MSIS® (NRL-SOF-014-1) SOFTWARE
// NRLMSIS® empirical atmospheric model software. Use is governed by the
// Open Source Academic Research License Agreement contained in the file
// nrlmsis2.1_license.txt, which is part of this software package. BY
// USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
// CONDITIONS OF THE LICENSE.
// #######################################################################

using System;
using System.Collections.Generic;

namespace NRLMSIS.Configuration
{
    /// <summary>
    /// Temperature model configuration parameters.
    /// Contains coefficients loaded from parameter file.
    /// </summary>
    public sealed class TemperatureConfiguration
    {
        /// <summary>Gets the temperature profile coefficients.</summary>
        public IReadOnlyList<double[]> Coefficients { get; }

        /// <summary>Creates a new temperature configuration.</summary>
        public TemperatureConfiguration(List<double[]> coefficients)
        {
            if (coefficients == null)
                throw new ArgumentNullException(nameof(coefficients));

            var copies = new List<double[]>(coefficients.Count);
            foreach (var coeff in coefficients)
            {
                copies.Add((double[])coeff.Clone());
            }
            Coefficients = copies.AsReadOnly();
        }
    }
}
