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
using NRLMSIS.Configuration;
using NRLMSIS.Calculators;
using NRLMSIS.Infrastructure;

namespace NRLMSIS.Models
{
    /// <summary>
    /// Orchestrates the MSIS calculation workflow.
    /// Coordinates various calculation steps and manages context/caching.
    /// </summary>
    internal sealed class MSISOrchestrator
    {
        private readonly MSISConfiguration _config;

        public MSISOrchestrator(MSISConfiguration config)
        {
            _config = config ?? throw new ArgumentNullException(nameof(config));
        }

        /// <summary>Orchestrates a complete MSIS calculation.</summary>
        public AtmosphericState Calculate(
            AtmosphericLocation location,
            SpaceWeather weather,
            Context context)
        {
            try
            {
                double geopotentialHeight = _config.UseGeodeticAltitude
                    ? Utilities.Alt2Gph(location.LatitudeDegrees, location.AltitudeKm)
                    : location.AltitudeKm;

                var basisFunctions = GetOrCalculateBasisFunctions(location, weather, context);
                var temperatureProfile = GetOrCalculateTemperatureProfile(basisFunctions, context);
                GetOrCalculateDensityProfiles(basisFunctions, temperatureProfile, context);

                if (geopotentialHeight < AltitudeRegimes.ZetaB)
                {
                    UpdateSplineWeightsIfNeeded(geopotentialHeight, context);
                }

                double temperature = MSISCalculator.CalculateTemperature(geopotentialHeight, context);
                var densities = CalculateDensities(geopotentialHeight, temperature, context);

                return new AtmosphericState
                {
                    Location = location,
                    Weather = weather,
                    Temperature = temperature,
                    ExosphericTemperature = temperatureProfile.Tex,
                    Densities = densities
                };
            }
            catch (Exception ex) when (!(ex is ArgumentException || ex is ArgumentOutOfRangeException))
            {
                throw new MsisCalculationException(
                    $"Calculation failed at location {location}",
                    location,
                    weather,
                    ex);
            }
        }

        private double[] GetOrCalculateBasisFunctions(
            AtmosphericLocation location,
            SpaceWeather weather,
            Context context)
        {
            bool needsRecalculation = context.IsLocationTimeChanged(
                location.DayOfYear,
                location.UniversalTimeSeconds,
                location.LatitudeDegrees,
                location.LongitudeDegrees,
                weather.F107Daily,
                weather.F107Average,
                weather.ApIndices);

            if (!needsRecalculation)
                return context.BasisFunctions;

            BasisFunctions.Globe(
                location.DayOfYear,
                location.UniversalTimeSeconds,
                location.LatitudeDegrees,
                location.LongitudeDegrees,
                weather.F107Average,
                weather.F107Daily,
                weather.ApIndices,
                out double[] gf);

            Array.Copy(gf, context.BasisFunctions, gf.Length);
            context.UpdateLocationTimeCache(
                location.DayOfYear,
                location.UniversalTimeSeconds,
                location.LatitudeDegrees,
                location.LongitudeDegrees,
                weather.F107Daily,
                weather.F107Average,
                weather.ApIndices);

            return context.BasisFunctions;
        }

        private TemperatureProfile GetOrCalculateTemperatureProfile(
            double[] basisFunctions,
            Context context)
        {
            TemperatureFunction.TemperatureParameters(basisFunctions, out var tpro);
            CopyTemperatureProfile(tpro, context.TemperatureProfile);
            return context.TemperatureProfile;
        }

        private void GetOrCalculateDensityProfiles(
            double[] basisFunctions,
            TemperatureProfile temperatureProfile,
            Context context)
        {
            for (int ispec = 2; ispec <= 10; ispec++)
            {
                if (Initialization.SpecFlag[ispec - 1])
                {
                    DensityProfile.DensityParameters(ispec, basisFunctions, temperatureProfile, out var dpro);
                    CopyDensityProfile(dpro, context.DensityProfiles[ispec]);
                }
            }
        }

        private void UpdateSplineWeightsIfNeeded(double geopotentialHeight, Context context)
        {
            if (!context.IsAltitudeChanged(geopotentialHeight))
                return;

            int maxOrder = geopotentialHeight < AltitudeRegimes.ZetaF ? 5 : 6;

            Utilities.BSpline(
                geopotentialHeight,
                Constants.NodesTN,
                Constants.Nd + 2,
                maxOrder,
                Initialization.EtaTN,
                out double[,] splineWeights,
                out int splineIndex);

            Array.Copy(splineWeights, context.SplineWeights, splineWeights.Length);
            context.SplineIndex = splineIndex;
            context.UpdateAltitudeCache(geopotentialHeight);
        }

        private SpeciesDensities CalculateDensities(
            double geopotentialHeight,
            double temperature,
            Context context)
        {
            var densities = new double[10];
            MSISCalculator.CalculateDensities(geopotentialHeight, temperature, context, densities);
            return SpeciesDensities.FromArray(densities);
        }

        private static void CopyTemperatureProfile(TemperatureProfile source, TemperatureProfile dest)
        {
            Array.Copy(source.Cf, dest.Cf, source.Cf.Length);
            Array.Copy(source.Beta, dest.Beta, source.Beta.Length);
            Array.Copy(source.Gamma, dest.Gamma, source.Gamma.Length);

            dest.TZetaF = source.TZetaF;
            dest.TZetaA = source.TZetaA;
            dest.DlnTdzA = source.DlnTdzA;
            dest.LnDTotF = source.LnDTotF;
            dest.Tex = source.Tex;
            dest.Tgb0 = source.Tgb0;
            dest.Tb0 = source.Tb0;
            dest.Sigma = source.Sigma;
            dest.SigmaSq = source.SigmaSq;
            dest.B = source.B;
            dest.CVs = source.CVs;
            dest.CVb = source.CVb;
            dest.CWs = source.CWs;
            dest.CWb = source.CWb;
            dest.VZetaF = source.VZetaF;
            dest.VZetaA = source.VZetaA;
            dest.WZetaA = source.WZetaA;
            dest.VZeta0 = source.VZeta0;
        }

        private static void CopyDensityProfile(DensityParameters source, DensityParameters dest)
        {
            if (source.Cf != null && dest.Cf != null)
                Array.Copy(source.Cf, dest.Cf, Math.Min(source.Cf.Length, dest.Cf.Length));
            if (source.Mi != null && dest.Mi != null)
                Array.Copy(source.Mi, dest.Mi, source.Mi.Length);
            if (source.ZetaMi != null && dest.ZetaMi != null)
                Array.Copy(source.ZetaMi, dest.ZetaMi, source.ZetaMi.Length);
            if (source.AMi != null && dest.AMi != null)
                Array.Copy(source.AMi, dest.AMi, source.AMi.Length);
            if (source.WMi != null && dest.WMi != null)
                Array.Copy(source.WMi, dest.WMi, source.WMi.Length);
            if (source.XMi != null && dest.XMi != null)
                Array.Copy(source.XMi, dest.XMi, source.XMi.Length);

            dest.LnPhiF = source.LnPhiF;
            dest.LnDRef = source.LnDRef;
            dest.ZetaM = source.ZetaM;
            dest.HML = source.HML;
            dest.HMU = source.HMU;
            dest.C = source.C;
            dest.ZetaC = source.ZetaC;
            dest.HC = source.HC;
            dest.R = source.R;
            dest.ZetaR = source.ZetaR;
            dest.HR = source.HR;
            dest.ZRef = source.ZRef;
            dest.IzRef = source.IzRef;
            dest.TRef = source.TRef;
            dest.ZMin = source.ZMin;
            dest.ZHyd = source.ZHyd;
            dest.ISpec = source.ISpec;
        }
    }
}
