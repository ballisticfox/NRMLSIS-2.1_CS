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
            if (!NRLMSIS.Initialization.InitFlag)
                throw new InvalidOperationException("Initialization.Initialize must be called first");

            var builder = new MSISConfiguration.Builder
            {
                SpeciesEnabled = (bool[])NRLMSIS.Initialization.SpecFlag.Clone(),
                MassWeights = (double[])NRLMSIS.Initialization.MassWgt.Clone(),
                LegacySwitches = NRLMSIS.Initialization.SwLeg.Select(f => (double)f).ToArray(),
                UseGeodeticAltitude = NRLMSIS.Initialization.ZAltFlag,
                UseNrlmsise00N2 = NRLMSIS.Initialization.N2RFlag,
                TemperatureEta = (double[,])NRLMSIS.Initialization.EtaTN.Clone(),
                O1Eta = (double[,])NRLMSIS.Initialization.EtaO1.Clone(),
                NOEta = (double[,])NRLMSIS.Initialization.EtaNO.Clone()
            };

            // Create temperature configuration from Initialization.TN (BasisSubset)
            var tempCoeffs = new List<double[]>();
            if (NRLMSIS.Initialization.TN?.Beta != null)
            {
                for (int i = 0; i < NRLMSIS.Initialization.TN.Beta.GetLength(0); i++)
                {
                    var row = new double[NRLMSIS.Initialization.TN.Beta.GetLength(1)];
                    for (int j = 0; j < row.Length; j++)
                    {
                        row[j] = NRLMSIS.Initialization.TN.Beta[i, j];
                    }
                    tempCoeffs.Add(row);
                }
            }
            builder.Temperature = new TemperatureConfiguration(tempCoeffs);

            // Create species configurations from BasisSubsets
            var speciesConfigs = new List<DensityConfiguration>();
            var speciesSubsets = new[] {
                (2, NRLMSIS.Initialization.N2),
                (3, NRLMSIS.Initialization.O2),
                (4, NRLMSIS.Initialization.O1),
                (5, NRLMSIS.Initialization.HE),
                (6, NRLMSIS.Initialization.H1),
                (7, NRLMSIS.Initialization.AR),
                (8, NRLMSIS.Initialization.N1),
                (9, NRLMSIS.Initialization.OA),
                (10, NRLMSIS.Initialization.NO)
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

            NRLMSIS.Initialization.SpecFlag = config.SpeciesEnabled.ToArray();
            NRLMSIS.Initialization.MassWgt = config.MassWeights.ToArray();
            NRLMSIS.Initialization.SwLeg = config.LegacySwitches.Select(d => (float)d).ToArray();
            NRLMSIS.Initialization.ZAltFlag = config.UseGeodeticAltitude;
            NRLMSIS.Initialization.N2RFlag = config.UseNrlmsise00N2;

            if (config.TemperatureEta != null)
                NRLMSIS.Initialization.EtaTN = (double[,])config.TemperatureEta.Clone();
            if (config.O1Eta != null)
                NRLMSIS.Initialization.EtaO1 = (double[,])config.O1Eta.Clone();
            if (config.NOEta != null)
                NRLMSIS.Initialization.EtaNO = (double[,])config.NOEta.Clone();
        }
    }
}