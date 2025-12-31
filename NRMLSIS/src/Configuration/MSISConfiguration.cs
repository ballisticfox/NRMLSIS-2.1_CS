// #######################################################################
// MSIS® (NRL-SOF-014-1) SOFTWARE
// NRLMSIS® empirical atmospheric model software. Use is governed by the
// Open Source Academic Research License Agreement contained in the file
// nrlmsis2.1_license.txt, which is part of this software package. BY
// USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
// CONDITIONS OF THE LICENSE.
// #######################################################################


using NRLMSIS.Infrastructure;

namespace NRLMSIS.Configuration
{
    /// <summary>
    /// Immutable configuration for NRLMSIS 2.1 model calculations.
    /// Create once during initialization and reuse for all calculations.
    /// Thread-safe due to immutability.
    /// </summary>
    public sealed class MSISConfiguration
    {
        /// <summary>Gets which atmospheric species are enabled for calculation (10 elements).</summary>
        public IReadOnlyList<bool> SpeciesEnabled { get; }

        /// <summary>Gets the mass weights for each species.</summary>
        public IReadOnlyList<double> MassWeights { get; }

        /// <summary>Gets the legacy NRLMSISE-00 compatibility switches (25 elements).</summary>
        public IReadOnlyList<double> LegacySwitches { get; }

        /// <summary>Gets whether to treat altitude input as geodetic altitude (true) or geopotential height (false).</summary>
        public bool UseGeodeticAltitude { get; }

        /// <summary>Gets whether to use NRLMSISE-00 N2 model (true) or MSIS 2.1 N2 model (false).</summary>
        public bool UseNrlmsise00N2 { get; }

        /// <summary>Gets the temperature profile parameters.</summary>
        public TemperatureConfiguration Temperature { get; }

        /// <summary>Gets the density profile parameters for each species.</summary>
        public IReadOnlyList<DensityConfiguration> Species { get; }

        /// <summary>Gets the B-spline eta for temperature nodes.</summary>
        public double[,]? TemperatureEta { get; }

        /// <summary>Gets the B-spline eta for O1 nodes.</summary>
        public double[,]? O1Eta { get; }

        /// <summary>Gets the B-spline eta for NO nodes.</summary>
        public double[,]? NOEta { get; }

        private MSISConfiguration(Builder builder)
        {
            if (builder.SpeciesEnabled == null || builder.SpeciesEnabled.Length != 10)
                throw new InvalidOperationException("Must specify exactly 10 species flags");
            if (builder.MassWeights == null || builder.MassWeights.Length != 10)
                throw new InvalidOperationException("Must specify exactly 10 mass weights");
            if (builder.LegacySwitches == null || builder.LegacySwitches.Length != 25)
                throw new InvalidOperationException("Must specify exactly 25 legacy switches");

            SpeciesEnabled = Array.AsReadOnly((bool[])builder.SpeciesEnabled.Clone());
            MassWeights = Array.AsReadOnly((double[])builder.MassWeights.Clone());
            LegacySwitches = Array.AsReadOnly((double[])builder.LegacySwitches.Clone());
            
            UseGeodeticAltitude = builder.UseGeodeticAltitude;
            UseNrlmsise00N2 = builder.UseNrlmsise00N2;

            Temperature = builder.Temperature ?? throw new InvalidOperationException("Temperature configuration is required");

            if (builder.Species == null || builder.Species.Count == 0)
                throw new InvalidOperationException("Species configurations are required");
            Species = builder.Species.AsReadOnly();

            if (builder.TemperatureEta != null)
                TemperatureEta = (double[,])builder.TemperatureEta.Clone();
            if (builder.O1Eta != null)
                O1Eta = (double[,])builder.O1Eta.Clone();
            if (builder.NOEta != null)
                NOEta = (double[,])builder.NOEta.Clone();
        }

        /// <summary>Creates a configuration with default settings.</summary>
        public static MSISConfiguration CreateDefault()
        {
            // If Initialization has been run, use its values
            if (Initialization.InitFlag)
            {
                return MsisConfigurationFactory.FromInitialization();
            }

            // Otherwise create minimal configuration
            // Note: This will have empty temperature/species configs and may not work for calculations
            // User should call Initialization.Initialize() first
            var builder = new Builder
            {
                SpeciesEnabled = new bool[10] { true, true, true, true, true, true, true, true, true, true },
                MassWeights = new double[10],
                LegacySwitches = Enumerable.Repeat(1.0, 25).ToArray(),
                UseGeodeticAltitude = true,
                UseNrlmsise00N2 = false,
                Temperature = new TemperatureConfiguration(new List<double[]>()),
                Species = new List<DensityConfiguration>()
            };

            for (int i = 0; i < 10; i++)
            {
                if (builder.SpeciesEnabled[i] && i > 0)
                    builder.MassWeights[i] = Constants.SpecMass[i + 1];
            }

            // Add empty species configs to satisfy validation
            for (int i = 2; i <= 10; i++)
            {
                builder.Species.Add(new DensityConfiguration(i, new List<double[]>()));
            }

            return builder.Build();
        }

        /// <summary>Builder for creating MsisConfiguration instances.</summary>
        public sealed class Builder
        {
            public bool[] SpeciesEnabled { get; set; }
            public double[] MassWeights { get; set; }
            public double[] LegacySwitches { get; set; }
            public bool UseGeodeticAltitude { get; set; }
            public bool UseNrlmsise00N2 { get; set; }
            public TemperatureConfiguration? Temperature { get; set; }
            public List<DensityConfiguration> Species { get; set; }
            public double[,]? TemperatureEta { get; set; }
            public double[,]? O1Eta { get; set; }
            public double[,]? NOEta { get; set; }

            public Builder()
            {
                SpeciesEnabled = new bool[10];
                MassWeights = new double[10];
                LegacySwitches = new double[25];
                Species = new List<DensityConfiguration>();
                UseGeodeticAltitude = true;
                UseNrlmsise00N2 = false;
            }

            public MSISConfiguration Build()
            {
                Validate();
                return new MSISConfiguration(this);
            }

            private void Validate()
            {
                if (SpeciesEnabled == null || SpeciesEnabled.Length != 10)
                    throw new InvalidOperationException("Must specify exactly 10 species flags");
                if (MassWeights == null || MassWeights.Length != 10)
                    throw new InvalidOperationException("Must specify exactly 10 mass weights");
                if (LegacySwitches == null || LegacySwitches.Length != 25)
                    throw new InvalidOperationException("Must specify exactly 25 legacy switches");
                if (Temperature == null)
                    throw new InvalidOperationException("Temperature configuration is required");
                if (Species == null || Species.Count == 0)
                    throw new InvalidOperationException("At least one species configuration is required");
            }
        }
    }
}