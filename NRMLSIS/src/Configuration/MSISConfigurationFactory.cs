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
using System.Linq;
using NRLMSIS.Infrastructure;

namespace NRLMSIS.Configuration
{
    /// <summary>
    /// Helper class to convert between legacy Initialization and immutable MsisConfiguration.
    /// Bridge for backward compatibility.
    /// </summary>
    public static class MsisConfigurationFactory
    {
        /// <summary>
        /// Creates a MsisConfiguration from the current Initialization state.
        /// </summary>
        public static MSISConfiguration FromInitialization()
        {
            if (!Initialization.InitFlag)
                throw new InvalidOperationException("Initialization.Initialize must be called first");

            var builder = new MSISConfiguration.Builder
            {
                SpeciesEnabled = (bool[])Initialization.SpecFlag.Clone(),
                MassWeights = (double[])Initialization.MassWgt.Clone(),
                LegacySwitches = Initialization.SwLeg.Select(f => (double)f).ToArray(),
                UseGeodeticAltitude = Initialization.ZAltFlag,
                UseNrlmsise00N2 = Initialization.N2RFlag,
                TemperatureEta = (double[,])Initialization.EtaTN.Clone(),
                O1Eta = (double[,])Initialization.EtaO1.Clone(),
                NOEta = (double[,])Initialization.EtaNO.Clone()
            };

            // Create temperature configuration from Initialization.TN (BasisSubset)
            var tempCoeffs = new List<double[]>();
            if (Initialization.TN?.Beta != null)
            {
                for (int i = 0; i < Initialization.TN.Beta.GetLength(0); i++)
                {
                    var row = new double[Initialization.TN.Beta.GetLength(1)];
                    for (int j = 0; j < row.Length; j++)
                    {
                        row[j] = Initialization.TN.Beta[i, j];
                    }
                    tempCoeffs.Add(row);
                }
            }
            builder.Temperature = new TemperatureConfiguration(tempCoeffs);

            // Create species configurations from BasisSubsets
            var speciesConfigs = new List<DensityConfiguration>();
            var speciesSubsets = new[] {
                (2, Initialization.N2),
                (3, Initialization.O2),
                (4, Initialization.O1),
                (5, Initialization.HE),
                (6, Initialization.H1),
                (7, Initialization.AR),
                (8, Initialization.N1),
                (9, Initialization.OA),
                (10, Initialization.NO)
            };

            foreach (var (specNum, subset) in speciesSubsets)
            {
                if (subset?.Beta != null)
                {
                    var coeffs = new List<double[]>();
                    for (int i = 0; i < subset.Beta.GetLength(0); i++)
                    {
                        var row = new double[subset.Beta.GetLength(1)];
                        for (int j = 0; j < row.Length; j++)
                        {
                            row[j] = subset.Beta[i, j];
                        }
                        coeffs.Add(row);
                    }
                    speciesConfigs.Add(new DensityConfiguration(specNum, coeffs));
                }
            }
            builder.Species = speciesConfigs;

            return builder.Build();
        }

        /// <summary>
        /// Applies a MsisConfiguration to the Initialization static state.
        /// For backward compatibility with code that still uses Initialization directly.
        /// Note: Does not modify InitFlag (it's read-only and managed by Initialization.Initialize).
        /// </summary>
        public static void ApplyToInitialization(MSISConfiguration config)
        {
            if (config == null)
                throw new ArgumentNullException(nameof(config));

            Initialization.SpecFlag = config.SpeciesEnabled.ToArray();
            Initialization.MassWgt = config.MassWeights.ToArray();
            Initialization.SwLeg = config.LegacySwitches.Select(d => (float)d).ToArray();
            Initialization.ZAltFlag = config.UseGeodeticAltitude;
            Initialization.N2RFlag = config.UseNrlmsise00N2;

            if (config.TemperatureEta != null)
                Initialization.EtaTN = (double[,])config.TemperatureEta.Clone();
            if (config.O1Eta != null)
                Initialization.EtaO1 = (double[,])config.O1Eta.Clone();
            if (config.NOEta != null)
                Initialization.EtaNO = (double[,])config.NOEta.Clone();
        }
    }
}