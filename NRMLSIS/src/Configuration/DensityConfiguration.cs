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
    /// Density model configuration parameters for a single species.
    /// Contains coefficients loaded from parameter file.
    /// </summary>
    public sealed class DensityConfiguration
    {
        /// <summary>Gets the species number (1-10).</summary>
        public int SpeciesNumber { get; }

        /// <summary>Gets the density profile coefficients.</summary>
        public IReadOnlyList<double[]> Coefficients { get; }

        /// <summary>Creates a new density configuration.</summary>
        public DensityConfiguration(int speciesNumber, List<double[]> coefficients)
        {
            if (speciesNumber < 1 || speciesNumber > 10)
                throw new ArgumentOutOfRangeException(nameof(speciesNumber), "Species number must be 1-10");
            ArgumentNullException.ThrowIfNull(coefficients);

            SpeciesNumber = speciesNumber;

            var copies = new List<double[]>(coefficients.Count);
            foreach (var coeff in coefficients)
            {
                copies.Add((double[])coeff.Clone());
            }
            Coefficients = copies.AsReadOnly();
        }
    }
}
