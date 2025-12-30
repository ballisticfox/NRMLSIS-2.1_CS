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
using NRLMSIS.Core;
using NRLMSIS.Configuration;
using NRLMSIS.Calculators;
using NRLMSIS.Infrastructure;

namespace NRLMSIS.Models
{
    /// <summary>
    /// NRLMSIS 2.1 implementation of the atmospheric model interface.
    /// Thread-safe for concurrent calculations when each thread uses its own Context.
    /// </summary>
    public sealed class NRLMSIS_Model : IAtmosphericModel
    {
        private readonly MSISConfiguration _config;

        /// <summary>Creates a new NRLMSIS 2.1 model instance.</summary>
        public NRLMSIS_Model(MSISConfiguration config)
        {
            _config = config ?? throw new ArgumentNullException(nameof(config));
            MsisConfigurationFactory.ApplyToInitialization(_config);
        }

        /// <summary>Creates a new NRLMSIS 2.1 model with default configuration.</summary>
        public static NRLMSIS_Model CreateDefault()
        {
            var config = MSISConfiguration.CreateDefault();
            return new NRLMSIS_Model(config);
        }

        /// <summary>Calculate atmospheric parameters at a single location.</summary>
        public AtmosphericState Calculate(AtmosphericLocation location, SpaceWeather weather)
        {
            using var context = new Context();
            return Calculate(location, weather, context);
        }

        /// <summary>Calculate atmospheric parameters using a provided context.</summary>
        public AtmosphericState Calculate(
            AtmosphericLocation location, 
            SpaceWeather weather, 
            Context context)
        {
            if (context == null)
                throw new ArgumentNullException(nameof(context));

            ValidateLocation(location);
            ValidateWeather(weather);

            MSISCalculator.Calculate(
                location.DayOfYear,
                location.UniversalTimeSeconds,
                location.AltitudeKm,
                location.LatitudeDegrees,
                location.LongitudeDegrees,
                weather.F107Average,
                weather.F107Daily,
                weather.ApIndices,
                context,
                out double temperature,
                out double[] densityArray,
                out double exosphericTemperature);

            var densities = SpeciesDensities.FromArray(densityArray);

            return new AtmosphericState
            {
                Location = location,
                Weather = weather,
                Temperature = temperature,
                ExosphericTemperature = exosphericTemperature,
                Densities = densities
            };
        }

        /// <summary>Calculate atmospheric parameters at multiple locations.</summary>
        public IEnumerable<AtmosphericState> CalculateBatch(
            IEnumerable<AtmosphericLocation> locations,
            SpaceWeather weather)
        {
            if (locations == null)
                throw new ArgumentNullException(nameof(locations));

            ValidateWeather(weather);

            using var context = new Context();

            foreach (var location in locations)
            {
                yield return Calculate(location, weather, context);
            }
        }

        private static void ValidateLocation(AtmosphericLocation location)
        {
            InputValidator.ValidateDayOfYear(location.DayOfYear);
            InputValidator.ValidateUniversalTime(location.UniversalTimeSeconds);
            InputValidator.ValidateAltitude(location.AltitudeKm);
            InputValidator.ValidateLatitude(location.LatitudeDegrees);
            InputValidator.ValidateLongitude(location.LongitudeDegrees);
        }

        private static void ValidateWeather(SpaceWeather weather)
        {
            InputValidator.ValidateF107Average(weather.F107Average);
            InputValidator.ValidateF107Daily(weather.F107Daily);
            InputValidator.ValidateApArray(weather.ApIndices);
        }
    }
}
