using System.Diagnostics;
using NRLMSIS.Core;
using NRLMSIS.Models;
using NRLMSIS.Infrastructure;

namespace NRLMSIS.Testing
{
    /// <summary>
    /// Performance evaluation script for NRLMSIS 2.1
    /// Tests two scenarios:
    /// 1. Fixed weather, varying locations (best-case caching)
    /// 2. Varying weather and locations (worst-case)
    /// </summary>
    class PerformanceEvaluation
    {
        private const int SAMPLES_PER_SET = 1000;
        private static readonly Random random = new Random(42); // Fixed seed for reproducibility

        public static void RunTest(string[] args)
        {
            Console.WriteLine("========================================================================");
            Console.WriteLine("           NRLMSIS 2.1 Performance Evaluation - 1000 Samples           ");
            Console.WriteLine("========================================================================");
            Console.WriteLine();

            // ========================================================================
            // INITIALIZATION
            // ========================================================================
            Console.WriteLine("Initializing model...");
            var initTimer = Stopwatch.StartNew();
            
            Initialization.Initialize(parmPath: "../", parmFile: "msis21.parm");
            var model = NRLMSIS_Model.CreateDefault();
            
            initTimer.Stop();
            Console.WriteLine($"Initialization complete: {initTimer.Elapsed.TotalMilliseconds:F2} ms");
            Console.WriteLine();

            // ========================================================================
            // TEST SET 1: Fixed Weather, Varying Locations
            // ========================================================================
            Console.WriteLine("========================================================================");
            Console.WriteLine("  TEST SET 1: Fixed Weather, Varying Locations (Best-case Caching)");
            Console.WriteLine("========================================================================");
            Console.WriteLine();
            Console.WriteLine("Scenario: Simulates scanning atmospheric grid at a fixed time");
            Console.WriteLine("          (e.g., mapping global atmosphere at one instant)");
            Console.WriteLine();

            var fixedWeather = CreateFixedWeather();
            var varyingLocations = GenerateVaryingLocations(SAMPLES_PER_SET);

            var set1Results = RunTestSet(
                model, 
                varyingLocations, 
                Enumerable.Repeat(fixedWeather, SAMPLES_PER_SET).ToList(),
                "Set 1"
            );

            // ========================================================================
            // TEST SET 2: Varying Weather and Locations
            // ========================================================================
            Console.WriteLine();
            Console.WriteLine("========================================================================");
            Console.WriteLine("  TEST SET 2: Varying Weather & Locations (Worst-case, No Caching)");
            Console.WriteLine("========================================================================");
            Console.WriteLine();
            Console.WriteLine("Scenario: Simulates time-series simulation with changing conditions");
            Console.WriteLine("          (e.g., rocket trajectory through time)");
            Console.WriteLine();

            var varyingWeather = GenerateVaryingWeather(SAMPLES_PER_SET);
            var varyingLocations2 = GenerateVaryingLocations(SAMPLES_PER_SET);

            var set2Results = RunTestSet(
                model,
                varyingLocations2,
                varyingWeather,
                "Set 2"
            );

            // ========================================================================
            // COMPARATIVE ANALYSIS
            // ========================================================================
            Console.WriteLine();
            Console.WriteLine("========================================================================");
            Console.WriteLine("                         Comparative Analysis                          ");
            Console.WriteLine("========================================================================");
            Console.WriteLine();

            PrintComparison(set1Results, set2Results);

            // ========================================================================
            // MEMORY ANALYSIS
            // ========================================================================
            Console.WriteLine();
            Console.WriteLine("========================================================================");
            Console.WriteLine("                          Memory Analysis                              ");
            Console.WriteLine("========================================================================");
            Console.WriteLine();

            PrintMemoryAnalysis();
        }

        // ========================================================================
        // TEST EXECUTION
        // ========================================================================

        static TestResults RunTestSet(
            NRLMSIS_Model model,
            List<AtmosphericLocation> locations,
            List<SpaceWeather> weatherConditions,
            string setName)
        {
            var results = new TestResults { SetName = setName };
            var sampleTimes = new List<double>();
            var sampleTimer = Stopwatch.StartNew();

            Console.WriteLine($"Running {setName}: {SAMPLES_PER_SET} calculations...");

            // Warm-up run (JIT compilation, cache warming)
            Console.WriteLine("  Warming up (10 samples)...");
            for (int i = 0; i < 10; i++)
            {
                model.Calculate(locations[i], weatherConditions[i]);
            }

            // Clear GC before actual test
            GC.Collect();
            GC.WaitForPendingFinalizers();
            GC.Collect();

            var gcGen0Before = GC.CollectionCount(0);
            var gcGen1Before = GC.CollectionCount(1);
            var gcGen2Before = GC.CollectionCount(2);
            var memoryBefore = GC.GetTotalMemory(false);

            Console.WriteLine("  Running timed test...");
            var totalTimer = Stopwatch.StartNew();

            using var context = new Context();

            for (int i = 0; i < SAMPLES_PER_SET; i++)
            {
                sampleTimer.Restart();
                
                // var state = model.Calculate(locations[i], weatherConditions[i]);
                var state = model.Calculate(locations[i], weatherConditions[i], context);
                
                sampleTimer.Stop();
                sampleTimes.Add(sampleTimer.Elapsed.TotalMilliseconds);

                // Store first result for validation
                if (i == 0)
                {
                    results.FirstTemperature = state.Temperature;
                    results.FirstDensity = state.Densities.TotalMassDensity;
                }

                // Progress indicator
                if ((i + 1) % 200 == 0)
                    Console.WriteLine($"    Progress: {i + 1}/{SAMPLES_PER_SET}");
            }

            totalTimer.Stop();

            var gcGen0After = GC.CollectionCount(0);
            var gcGen1After = GC.CollectionCount(1);
            var gcGen2After = GC.CollectionCount(2);
            var memoryAfter = GC.GetTotalMemory(false);

            // Calculate statistics
            results.TotalTime = totalTimer.Elapsed;
            results.SampleTimes = sampleTimes;
            results.SamplesProcessed = SAMPLES_PER_SET;
            results.GCGen0 = gcGen0After - gcGen0Before;
            results.GCGen1 = gcGen1After - gcGen1Before;
            results.GCGen2 = gcGen2After - gcGen2Before;
            results.MemoryDelta = memoryAfter - memoryBefore;

            Console.WriteLine($"  Completed in {results.TotalTime.TotalMilliseconds:F2} ms");
            Console.WriteLine();

            PrintTestResults(results);

            return results;
        }

        // ========================================================================
        // RESULT PRINTING
        // ========================================================================

        static void PrintTestResults(TestResults results)
        {
            if (results.SampleTimes == null)
                return;

            var sortedTimes = results.SampleTimes.OrderBy(t => t).ToList();
            double minTime = sortedTimes.First();
            double maxTime = sortedTimes.Last();
            double avgTime = results.SampleTimes.Average();
            double medianTime = sortedTimes[sortedTimes.Count / 2];
            double p95Time = sortedTimes[(int)(sortedTimes.Count * 0.95)];
            double p99Time = sortedTimes[(int)(sortedTimes.Count * 0.99)];
            double stdDev = CalculateStdDev(results.SampleTimes, avgTime);

            Console.WriteLine($"Results for {results.SetName}:");
            Console.WriteLine("------------------------------------------------------------------------");
            Console.WriteLine($"  Total Time:                {results.TotalTime.TotalMilliseconds:F2} ms");
            Console.WriteLine($"  Total Time (seconds):      {results.TotalTime.TotalSeconds:F3} s");
            Console.WriteLine();
            Console.WriteLine("  Per-Sample Statistics:");
            Console.WriteLine($"    Average:                 {avgTime:F4} ms");
            Console.WriteLine($"    Median:                  {medianTime:F4} ms");
            Console.WriteLine($"    Min:                     {minTime:F4} ms");
            Console.WriteLine($"    Max:                     {maxTime:F4} ms");
            Console.WriteLine($"    95th Percentile:         {p95Time:F4} ms");
            Console.WriteLine($"    99th Percentile:         {p99Time:F4} ms");
            Console.WriteLine($"    Std Deviation:           {stdDev:F4} ms");
            Console.WriteLine($"    Coeff of Variation:      {(stdDev / avgTime * 100):F2}%");
            Console.WriteLine();
            Console.WriteLine("  Throughput:");
            Console.WriteLine($"    Samples/Second:          {results.SamplesProcessed / results.TotalTime.TotalSeconds:F1}");
            Console.WriteLine($"    Calculations/Minute:     {results.SamplesProcessed / results.TotalTime.TotalSeconds * 60:F0}");
            Console.WriteLine();
            Console.WriteLine("  GC Activity:");
            Console.WriteLine($"    Gen 0 Collections:       {results.GCGen0}");
            Console.WriteLine($"    Gen 1 Collections:       {results.GCGen1}");
            Console.WriteLine($"    Gen 2 Collections:       {results.GCGen2}");
            Console.WriteLine($"    Memory Delta:            {results.MemoryDelta / 1024.0 / 1024.0:F2} MB");
            Console.WriteLine($"    Memory per Sample:       {results.MemoryDelta / (double)results.SamplesProcessed / 1024.0:F2} KB");
            Console.WriteLine();
            Console.WriteLine("  Sample Output (validation):");
            Console.WriteLine($"    First Temperature:       {results.FirstTemperature:F2} K");
            Console.WriteLine($"    First Density:           {results.FirstDensity:E4} kg/m³");
            Console.WriteLine();
        }

        static void PrintComparison(TestResults set1, TestResults set2)
        {
            if (set1.SampleTimes == null || set2.SampleTimes == null)
                return;

            double set1Avg = set1.SampleTimes.Average();
            double set2Avg = set2.SampleTimes.Average();
            double speedRatio = set2Avg / set1Avg;
            double timeRatio = set2.TotalTime.TotalMilliseconds / set1.TotalTime.TotalMilliseconds;

            Console.WriteLine("Set 1 vs Set 2 Comparison:");
            Console.WriteLine("------------------------------------------------------------------------");
            Console.WriteLine($"                                 Set 1          Set 2          Ratio");
            Console.WriteLine("------------------------------------------------------------------------");
            Console.WriteLine($"  Total Time (ms):            {set1.TotalTime.TotalMilliseconds,9:F2}    {set2.TotalTime.TotalMilliseconds,9:F2}    {timeRatio,7:F3}x");
            Console.WriteLine($"  Average Time (ms):          {set1Avg,9:F4}    {set2Avg,9:F4}    {speedRatio,7:F3}x");
            Console.WriteLine($"  Throughput (samples/s):     {set1.SamplesProcessed / set1.TotalTime.TotalSeconds,9:F1}    {set2.SamplesProcessed / set2.TotalTime.TotalSeconds,9:F1}");
            Console.WriteLine($"  Gen 0 GC:                   {set1.GCGen0,9}    {set2.GCGen0,9}");
            Console.WriteLine($"  Gen 1 GC:                   {set1.GCGen1,9}    {set2.GCGen1,9}");
            Console.WriteLine($"  Gen 2 GC:                   {set2.GCGen2,9}    {set2.GCGen2,9}");
            Console.WriteLine($"  Memory Delta (MB):          {set1.MemoryDelta / 1024.0 / 1024.0,9:F2}    {set2.MemoryDelta / 1024.0 / 1024.0,9:F2}");
            Console.WriteLine();
        }

        static void PrintMemoryAnalysis()
        {
            long totalMemory = GC.GetTotalMemory(false);
            Console.WriteLine("Overall Memory Statistics:");
            Console.WriteLine("------------------------------------------------------------------------");
            Console.WriteLine($"  Total Managed Memory:      {totalMemory / 1024.0 / 1024.0:F2} MB");
            Console.WriteLine($"  Total GC Gen 0:            {GC.CollectionCount(0)}");
            Console.WriteLine($"  Total GC Gen 1:            {GC.CollectionCount(1)}");
            Console.WriteLine($"  Total GC Gen 2:            {GC.CollectionCount(2)}");
            Console.WriteLine();

            // Force GC to see actual memory usage
            GC.Collect();
            GC.WaitForPendingFinalizers();
            GC.Collect();

            long memoryAfterGC = GC.GetTotalMemory(false);
            Console.WriteLine($"  After Full GC:             {memoryAfterGC / 1024.0 / 1024.0:F2} MB");
            Console.WriteLine($"  Collectable Garbage:       {(totalMemory - memoryAfterGC) / 1024.0 / 1024.0:F2} MB");
        }

        // ========================================================================
        // DATA GENERATION
        // ========================================================================

        static SpaceWeather CreateFixedWeather()
        {
            // Moderate solar activity conditions
            return new SpaceWeather(
                f107Average: 150.0,
                f107Daily: 145.0,
                apIndices: new double[] { 15, 12, 10, 8, 6, 4, 3 }
            );
        }

        static List<SpaceWeather> GenerateVaryingWeather(int count)
        {
            var weatherList = new List<SpaceWeather>();
            
            for (int i = 0; i < count; i++)
            {
                // Generate varying solar and geomagnetic conditions
                double f107a = 80 + random.NextDouble() * 170;  // 80-250
                double f107 = f107a * (0.8 + random.NextDouble() * 0.4); // ±20% of average
                
                var apIndices = new double[7];
                for (int j = 0; j < 7; j++)
                {
                    apIndices[j] = random.NextDouble() * 100; // 0-100
                }

                weatherList.Add(new SpaceWeather(
                    f107Average: f107a,
                    f107Daily: f107,
                    apIndices: apIndices
                ));
            }

            return weatherList;
        }

        static List<AtmosphericLocation> GenerateVaryingLocations(int count)
        {
            var locations = new List<AtmosphericLocation>();

            for (int i = 0; i < count; i++)
            {
                // Generate random location across globe and altitude range
                int dayOfYear = 1 + random.Next(365);
                double utSec = random.NextDouble() * 86400; // 0-24 hours in seconds
                double lat = -90 + random.NextDouble() * 180; // -90 to +90
                double lon = -180 + random.NextDouble() * 360; // -180 to +180
                double alt = random.NextDouble() * 500; // 0-500 km

                locations.Add(new AtmosphericLocation(
                    dayOfYear: dayOfYear,
                    universalTimeSeconds: utSec,
                    latitudeDegrees: lat,
                    longitudeDegrees: lon,
                    altitudeKm: alt
                ));
            }

            return locations;
        }

        // ========================================================================
        // UTILITIES
        // ========================================================================

        static double CalculateStdDev(List<double> values, double mean)
        {
            double sumSquaredDiff = values.Sum(v => Math.Pow(v - mean, 2));
            return Math.Sqrt(sumSquaredDiff / values.Count);
        }

        // ========================================================================
        // RESULT STORAGE
        // ========================================================================

        class TestResults
        {
            public string? SetName { get; set; }
            public TimeSpan TotalTime { get; set; }
            public List<double>? SampleTimes { get; set; }
            public int SamplesProcessed { get; set; }
            public int GCGen0 { get; set; }
            public int GCGen1 { get; set; }
            public int GCGen2 { get; set; }
            public long MemoryDelta { get; set; }
            public double FirstTemperature { get; set; }
            public double FirstDensity { get; set; }
        }
    }
}