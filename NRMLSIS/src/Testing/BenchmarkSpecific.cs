using System;
using NRLMSIS.Core;
using NRLMSIS.Models;
using NRLMSIS.Infrastructure;
using NRLMSIS.Calculators;
using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Running;
using BenchmarkDotNet.Configs;
using BenchmarkDotNet.Jobs;

namespace NRLMSIS.Testing
{
    /// <summary>
    /// DEEP DIVE Benchmarks - Drilling into the two hotspots:
    /// 1. TemperatureFunction.TemperatureParameters (~8.4 μs)
    /// 2. DensityProfile.DensityParameters (all species) (~24.8 μs)
    /// 
    /// FIXED VERSION - Removed crashing GeoMag/UtDep benchmarks
    /// Run with: dotnet run -c Release
    /// </summary>
    [Config(typeof(BenchmarkConfig))]
    [MemoryDiagnoser]
    public class MSISDeepDiveBenchmarks
    {
        private class BenchmarkConfig : ManualConfig
        {
            public BenchmarkConfig()
            {
                AddJob(Job.Default
                    .WithWarmupCount(3)
                    .WithIterationCount(10)
                );
            }
        }

        private NRLMSIS_Model _model;
        private AtmosphericLocation _location;
        private SpaceWeather _weather;
        private Context _context;
        private Context _cachedContext;
        
        // Pre-calculated values
        private double[] _basisFunctions;
        private TemperatureProfile _temperatureProfile;
        private double _geopotentialHeight;
        private double _temperature;
        
        // For deeper testing - extracted subsets
        private double[] _geomagPrimary;
        private double[,] _geomagSecondary;
        private double[] _utDependent;

        [GlobalSetup]
        public void GlobalSetup()
        {
            Console.WriteLine("Initializing for DEEP DIVE benchmarking...");
            
            // ADJUST THIS PATH TO YOUR SYSTEM
            Initialization.Initialize(parmPath: "", parmFile: "msis21.parm");
            _model = NRLMSIS_Model.CreateDefault();

            _location = new AtmosphericLocation(
                dayOfYear: 172,
                universalTimeSeconds: 43200,
                latitudeDegrees: 45.0,
                longitudeDegrees: -75.0,
                altitudeKm: 400.0
            );

            _weather = new SpaceWeather(
                f107Average: 150.0,
                f107Daily: 150.0,
                apIndices: new double[] { 4, 0, 0, 0, 0, 0, 0 }
            );

            _geopotentialHeight = Utilities.Alt2Gph(_location.LatitudeDegrees, _location.AltitudeKm);
            _context = new Context();
            _cachedContext = new Context();
            
            _model.Calculate(_location, _weather, _cachedContext);
            
            // Pre-calculate basis functions
            BasisFunctions.Globe(
                _location.DayOfYear,
                _location.UniversalTimeSeconds,
                _location.LatitudeDegrees,
                _location.LongitudeDegrees,
                _weather.F107Average,
                _weather.F107Daily,
                _weather.ApIndices,
                out _basisFunctions
            );
            
            // Pre-calculate temperature profile
            TemperatureFunction.TemperatureParameters(_basisFunctions, out _temperatureProfile);
            _temperature = MSISCalculator.CalculateTemperature(_geopotentialHeight, _cachedContext);
            
            // Extract basis function subsets for deeper testing
            _geomagPrimary = BasisFunctionExtractor.ExtractGeomagneticPrimary(_basisFunctions);
            _geomagSecondary = BasisFunctionExtractor.ExtractGeomagneticSecondary(_basisFunctions);
            _utDependent = BasisFunctionExtractor.ExtractUniversalTimeDependent(_basisFunctions);
            
            Console.WriteLine("Setup complete - ready for deep dive!");
            Console.WriteLine();
        }

        [GlobalCleanup]
        public void GlobalCleanup()
        {
            _context?.Dispose();
            _cachedContext?.Dispose();
        }

        // ====================================================================
        // BASELINE BENCHMARKS (for reference)
        // ====================================================================

        [Benchmark(Baseline = true, Description = "BASELINE: Full calculation")]
        public AtmosphericState Baseline_FullCalculation()
        {
            return _model.Calculate(_location, _weather);
        }

        // ====================================================================
        // HOTSPOT 1: TemperatureFunction.TemperatureParameters (~8.4 μs)
        // ====================================================================

        [Benchmark(Description = "HOTSPOT 1: TemperatureFunction (total)")]
        public TemperatureProfile Hotspot1_TemperatureFunction_Total()
        {
            TemperatureFunction.TemperatureParameters(_basisFunctions, out var tempProfile);
            return tempProfile;
        }

        [Benchmark(Description = "  └─ BasisFunction extraction (3 calls)")]
        public void Hotspot1_BasisFunctionExtraction()
        {
            var geomagPrimary = BasisFunctionExtractor.ExtractGeomagneticPrimary(_basisFunctions);
            var geomagSecondary = BasisFunctionExtractor.ExtractGeomagneticSecondary(_basisFunctions);
            var utDependent = BasisFunctionExtractor.ExtractUniversalTimeDependent(_basisFunctions);
        }

        [Benchmark(Description = "  └─ BasisFunctions.SFluxMod()")]
        public double Hotspot1_SFluxMod()
        {
            // SFluxMod is called multiple times in temperature calculation
            var subset = Initialization.TN;
            int columnIndex = Constants.ItEx;
            double normalizationFactor = 1.0 / subset.Beta[0, columnIndex - subset.Bl];
            return BasisFunctions.SFluxMod(columnIndex, _basisFunctions, subset, normalizationFactor);
        }

        // NOTE: GeoMag and UtDep benchmarks removed due to array indexing issues
        // These functions are called internally but require complex parameter arrays
        // The overhead from these functions is captured in the total TemperatureFunction benchmark

        // ====================================================================
        // HOTSPOT 2: DensityProfile.DensityParameters (all) (~24.8 μs)
        // ====================================================================

        [Benchmark(Description = "HOTSPOT 2: DensityProfile all species (total)")]
        public void Hotspot2_DensityProfile_AllSpecies()
        {
            // Species: 2=N2, 3=O2, 4=O, 5=He, 6=H, 7=Ar, 8=N, 9=O*, 10=NO
            for (int ispec = 2; ispec <= 10; ispec++)
            {
                if (Initialization.SpecFlag[ispec - 1])
                {
                    DensityProfile.DensityParameters(
                        ispec,
                        _basisFunctions,
                        _temperatureProfile,
                        out var densityParams
                    );
                }
            }
        }

        [Benchmark(Description = "  └─ DensityProfile: N2 (species 2)")]
        public DensityParameters Hotspot2_DensityProfile_N2()
        {
            DensityProfile.DensityParameters(2, _basisFunctions, _temperatureProfile, out var dp);
            return dp;
        }

        [Benchmark(Description = "  └─ DensityProfile: O2 (species 3)")]
        public DensityParameters Hotspot2_DensityProfile_O2()
        {
            DensityProfile.DensityParameters(3, _basisFunctions, _temperatureProfile, out var dp);
            return dp;
        }

        [Benchmark(Description = "  └─ DensityProfile: O (species 4)")]
        public DensityParameters Hotspot2_DensityProfile_O()
        {
            DensityProfile.DensityParameters(4, _basisFunctions, _temperatureProfile, out var dp);
            return dp;
        }

        [Benchmark(Description = "  └─ DensityProfile: He (species 5)")]
        public DensityParameters Hotspot2_DensityProfile_He()
        {
            DensityProfile.DensityParameters(5, _basisFunctions, _temperatureProfile, out var dp);
            return dp;
        }

        [Benchmark(Description = "  └─ DensityProfile: H (species 6)")]
        public DensityParameters Hotspot2_DensityProfile_H()
        {
            DensityProfile.DensityParameters(6, _basisFunctions, _temperatureProfile, out var dp);
            return dp;
        }

        [Benchmark(Description = "  └─ DensityProfile: Ar (species 7)")]
        public DensityParameters Hotspot2_DensityProfile_Ar()
        {
            DensityProfile.DensityParameters(7, _basisFunctions, _temperatureProfile, out var dp);
            return dp;
        }

        [Benchmark(Description = "  └─ DensityProfile: N (species 8)")]
        public DensityParameters Hotspot2_DensityProfile_N()
        {
            DensityProfile.DensityParameters(8, _basisFunctions, _temperatureProfile, out var dp);
            return dp;
        }

        [Benchmark(Description = "  └─ DensityProfile: O* (species 9)")]
        public DensityParameters Hotspot2_DensityProfile_OAnomolous()
        {
            DensityProfile.DensityParameters(9, _basisFunctions, _temperatureProfile, out var dp);
            return dp;
        }

        [Benchmark(Description = "  └─ DensityProfile: NO (species 10)")]
        public DensityParameters Hotspot2_DensityProfile_NO()
        {
            DensityProfile.DensityParameters(10, _basisFunctions, _temperatureProfile, out var dp);
            return dp;
        }

        // ====================================================================
        // COMPONENT ANALYSIS: What's inside DensityProfile?
        // ====================================================================

        [Benchmark(Description = "  └─ BasisFunction extraction (per species)")]
        public void Hotspot2_BasisFunctionExtraction_PerSpecies()
        {
            // Each species does this extraction
            var geomagPrimary = BasisFunctionExtractor.ExtractGeomagneticPrimary(_basisFunctions);
            var geomagSecondary = BasisFunctionExtractor.ExtractGeomagneticSecondary(_basisFunctions);
            var utDependent = BasisFunctionExtractor.ExtractUniversalTimeDependent(_basisFunctions);
        }

        // ====================================================================
        // MATH PRIMITIVES - Testing the core operations
        // ====================================================================

        [Benchmark(Description = "Math: Array allocation (double[50])")]
        public double[] Math_ArrayAllocation()
        {
            return new double[50];  // Typical size in MSIS
        }

        [Benchmark(Description = "Math: Array.Copy (50 elements)")]
        public void Math_ArrayCopy()
        {
            double[] source = _basisFunctions;
            double[] dest = new double[50];
            Array.Copy(source, dest, 50);
        }

        // ====================================================================
        // COMPARISON BENCHMARKS
        // ====================================================================

        [Benchmark(Description = "COMPARE: BasisFunctions.Globe()")]
        public double[] Compare_BasisFunctions()
        {
            BasisFunctions.Globe(
                _location.DayOfYear,
                _location.UniversalTimeSeconds,
                _location.LatitudeDegrees,
                _location.LongitudeDegrees,
                _weather.F107Average,
                _weather.F107Daily,
                _weather.ApIndices,
                out double[] gf
            );
            return gf;
        }

        [Benchmark(Description = "COMPARE: DensityCalculator.CalculateAll()")]
        public SpeciesDensities Compare_DensityCalculation()
        {
            return DensityCalculator.CalculateAll(
                _geopotentialHeight,
                _temperature,
                _cachedContext
            );
        }
    }

    class BenchmarkSpecific
    {
        public static void RunBenchmark()
        {
            Console.WriteLine("╔══════════════════════════════════════════════════════════════════════╗");
            Console.WriteLine("║              NRLMSIS 2.1 DEEP DIVE Benchmarking                      ║");
            Console.WriteLine("║                  Hotspot Analysis Suite                              ║");
            Console.WriteLine("╚══════════════════════════════════════════════════════════════════════╝");
            Console.WriteLine();

            #if DEBUG
            Console.ForegroundColor = ConsoleColor.Red;
            Console.WriteLine("⚠️  WARNING: Running in DEBUG mode!");
            Console.WriteLine("   Run with: dotnet run -c Release");
            Console.ResetColor();
            Console.WriteLine();
            return;
            #endif

            Console.WriteLine("🔍 Deep Dive Focus Areas:");
            Console.WriteLine();
            Console.WriteLine("HOTSPOT 1: TemperatureFunction.TemperatureParameters (~8.4 μs)");
            Console.WriteLine("  • Testing internal function calls");
            Console.WriteLine("  • SFluxMod breakdown");
            Console.WriteLine("  • Basis function extraction overhead");
            Console.WriteLine("  • (Note: GeoMag/UtDep removed due to parameter complexity)");
            Console.WriteLine();
            Console.WriteLine("HOTSPOT 2: DensityProfile.DensityParameters (~24.8 μs for all species)");
            Console.WriteLine("  • Per-species breakdown (N2, O2, O, He, H, Ar, N, O*, NO)");
            Console.WriteLine("  • Identify most expensive species");
            Console.WriteLine("  • Extraction overhead analysis");
            Console.WriteLine();
            Console.WriteLine("This will take 5-7 minutes...");
            Console.WriteLine();

            var summary = BenchmarkRunner.Run<MSISDeepDiveBenchmarks>();

            Console.WriteLine();
            Console.WriteLine("╔══════════════════════════════════════════════════════════════════════╗");
            Console.WriteLine("║                    Deep Dive Complete!                               ║");
            Console.WriteLine("╚══════════════════════════════════════════════════════════════════════╝");
        }
    }
}