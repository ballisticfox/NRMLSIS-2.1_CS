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
    /// NRLMSIS 2.1 Performance Benchmarking with BenchmarkDotNet
    ///
    /// VERIFIED - All method signatures correct!
    ///
    /// Run with: dotnet run -c Release
    ///
    /// This benchmarks your exact test case from the simple script:
    /// - June 21 (DOY 172), 12:00 UT
    /// - 45°N, 75°W, 400 km altitude
    /// - F10.7 = 150, Ap = 4
    /// </summary>
    [Config(typeof(BenchmarkConfig))]
    [MemoryDiagnoser]
    public class MSISBenchmarks
    {
        private class BenchmarkConfig : ManualConfig
        {
            public BenchmarkConfig()
            {
                AddJob(Job.Default
                    .WithWarmupCount(3)      // 3 warmup runs
                    .WithIterationCount(10)  // 10 measured runs
                );
            }
        }

        private NRLMSIS_Model _model;
        private AtmosphericLocation _location;
        private SpaceWeather _weather;
        private Context _context;
        private Context _cachedContext;
        
        // Pre-calculated values for isolated benchmarks
        private double[] _basisFunctions;
        private TemperatureProfile _temperatureProfile;
        private double _geopotentialHeight;
        private double _temperature;

        [GlobalSetup]
        public void GlobalSetup()
        {
            Console.WriteLine("Initializing NRLMSIS model...");
            
            // Initialize model (loads msis21.parm)
            Initialization.Initialize(parmPath: "", parmFile: "msis21.parm");
            _model = NRLMSIS_Model.CreateDefault();

            // Create test location (same as simple script)
            _location = new AtmosphericLocation(
                dayOfYear: 172,
                universalTimeSeconds: 43200,
                latitudeDegrees: 45.0,
                longitudeDegrees: -75.0,
                altitudeKm: 400.0
            );

            // Create space weather (same as simple script)
            _weather = new SpaceWeather(
                f107Average: 150.0,
                f107Daily: 150.0,
                apIndices: new double[] { 4, 0, 0, 0, 0, 0, 0 }
            );

            // Calculate geopotential height
            _geopotentialHeight = Utilities.Alt2Gph(_location.LatitudeDegrees, _location.AltitudeKm);

            // Create contexts
            _context = new Context();
            _cachedContext = new Context();
            
            // Warm up cached context
            var warmupState = _model.Calculate(_location, _weather, _cachedContext);
            
            // Pre-calculate basis functions for isolated tests
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
            
            // Pre-calculate temperature
            _temperature = MSISCalculator.CalculateTemperature(_geopotentialHeight, _cachedContext);

            Console.WriteLine("Setup complete!");
            Console.WriteLine($"Test: DOY {_location.DayOfYear}, {_location.AltitudeKm} km, F10.7={_weather.F107Daily}");
            Console.WriteLine($"Expected temp: ~{warmupState.Temperature:F0} K");
            Console.WriteLine($"Expected O density: ~{warmupState.Densities.O:E2} m⁻³");
            Console.WriteLine();
        }

        [GlobalCleanup]
        public void GlobalCleanup()
        {
            _context?.Dispose();
            _cachedContext?.Dispose();
        }

        // ====================================================================
        // FULL CALCULATION BENCHMARKS
        // ====================================================================

        [Benchmark(Baseline = true, Description = "Full calc (no context)")]
        public AtmosphericState FullCalculation_NoContextReuse()
        {
            // This is what your simple script does
            return _model.Calculate(_location, _weather);
        }

        [Benchmark(Description = "Full calc (context, no cache)")]
        public AtmosphericState FullCalculation_ContextNoCache()
        {
            _context.Reset();
            return _model.Calculate(_location, _weather, _context);
        }

        [Benchmark(Description = "Full calc (cached)")]
        public AtmosphericState FullCalculation_Cached()
        {
            // This is what you get with proper context reuse
            return _model.Calculate(_location, _weather, _cachedContext);
        }

        // ====================================================================
        // COMPONENT BENCHMARKS - Where does time go?
        // ====================================================================

        [Benchmark(Description = "BasisFunctions.Globe()")]
        public double[] BenchmarkBasisFunctions()
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

        [Benchmark(Description = "TemperatureFunction params")]
        public TemperatureProfile BenchmarkTemperatureProfile()
        {
            // CORRECT: TemperatureFunction, not TemperatureProfile
            TemperatureFunction.TemperatureParameters(_basisFunctions, out var tempProfile);
            return tempProfile;
        }

        [Benchmark(Description = "DensityProfile (O2 only)")]
        public DensityParameters BenchmarkDensityProfile_O2()
        {
            DensityProfile.DensityParameters(
                2,  // O2 species
                _basisFunctions,
                _temperatureProfile,
                out var densityParams
            );
            return densityParams;
        }

        [Benchmark(Description = "DensityProfile (all species)")]
        public void BenchmarkDensityProfile_AllSpecies()
        {
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

        // ====================================================================
        // EVALUATION BENCHMARKS - Final calculations
        // ====================================================================

        [Benchmark(Description = "Temperature evaluation")]
        public double BenchmarkTemperatureEvaluation()
        {
            return MSISCalculator.CalculateTemperature(_geopotentialHeight, _cachedContext);
        }

        [Benchmark(Description = "Density evaluation (all)")]
        public SpeciesDensities BenchmarkDensityEvaluation()
        {
            // CORRECT: DensityCalculator.CalculateAll, not Calculate
            return DensityCalculator.CalculateAll(
                _geopotentialHeight,
                _temperature,
                _cachedContext
            );
        }

        // ====================================================================
        // UTILITY BENCHMARKS
        // ====================================================================

        [Benchmark(Description = "Alt2Gph conversion")]
        public double BenchmarkAlt2Gph()
        {
            return Utilities.Alt2Gph(_location.LatitudeDegrees, _location.AltitudeKm);
        }

        [Benchmark(Description = "BSpline calculation")]
        public void BenchmarkBSpline()
        {
            int maxOrder = _geopotentialHeight < AltitudeRegimes.ZetaF ? 5 : 6;

            // CORRECT: All parameters including Initialization.EtaTN
            Utilities.BSpline(
                _geopotentialHeight,
                Constants.NodesTN,
                Constants.Nd + 2,
                maxOrder,
                Initialization.EtaTN,  // eta coefficients
                out double[,] S,        // 2D spline array
                out int i               // interval index
            );
        }
    }

    class Benchmark
    {
        public static void RunBenchmark()
        {
            Console.WriteLine("╔══════════════════════════════════════════════════════════════════════╗");
            Console.WriteLine("║          NRLMSIS 2.1 Performance Benchmarking Suite                  ║");
            Console.WriteLine("╚══════════════════════════════════════════════════════════════════════╝");
            Console.WriteLine();

            #if DEBUG
            Console.ForegroundColor = ConsoleColor.Red;
            Console.WriteLine("   WARNING: Running in DEBUG mode!");
            Console.WriteLine("   Benchmarks will be inaccurate.");
            Console.WriteLine("   Run with: dotnet run -c Release");
            Console.ResetColor();
            Console.WriteLine();
            return;
            #endif

            Console.WriteLine("Benchmark Suite:");
            Console.WriteLine("  • Full calculations (3 variations)");
            Console.WriteLine("  • Component methods (4 tests)");
            Console.WriteLine("  • Evaluation methods (2 tests)");
            Console.WriteLine("  • Utilities (2 tests)");
            Console.WriteLine();
            Console.WriteLine("This will take 3-5 minutes...");
            Console.WriteLine();

            var summary = BenchmarkRunner.Run<MSISBenchmarks>();

            Console.WriteLine();
            Console.WriteLine("╔══════════════════════════════════════════════════════════════════════╗");
            Console.WriteLine("║                      Benchmarking Complete!                          ║");
            Console.WriteLine("╚══════════════════════════════════════════════════════════════════════╝");
            Console.WriteLine();
        }
    }
}