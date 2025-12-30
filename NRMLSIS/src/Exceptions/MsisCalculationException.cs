// #######################################################################
// MSIS® (NRL-SOF-014-1) SOFTWARE
// NRLMSIS® empirical atmospheric model software. Use is governed by the
// Open Source Academic Research License Agreement contained in the file
// nrlmsis2.1_license.txt, which is part of this software package. BY
// USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
// CONDITIONS OF THE LICENSE.
// #######################################################################

using System;
using NRLMSIS.Core;

namespace NRLMSIS.Models
{
    /// <summary>
    /// Exception thrown when MSIS calculation fails.
    /// Contains context about the failed calculation.
    /// </summary>
    public class MsisCalculationException : Exception
    {
        /// <summary>Gets the location where the calculation failed.</summary>
        public AtmosphericLocation? Location { get; }

        /// <summary>Gets the space weather conditions when the calculation failed.</summary>
        public SpaceWeather? Weather { get; }

        public MsisCalculationException(string message) : base(message)
        {
        }

        public MsisCalculationException(string message, Exception innerException)
            : base(message, innerException)
        {
        }

        public MsisCalculationException(
            string message,
            AtmosphericLocation location,
            SpaceWeather weather,
            Exception innerException = null)
            : base(message, innerException)
        {
            Location = location;
            Weather = weather;
        }

        public override string ToString()
        {
            var str = base.ToString();
            if (Location.HasValue)
                str += $"\nLocation: {Location.Value}";
            if (Weather.HasValue)
                str += $"\nWeather: {Weather.Value}";
            return str;
        }
    }
}
