using System.Diagnostics;
using System.Globalization;
using NRLMSIS.Core;
using NRLMSIS.Models;
using NRLMSIS.Infrastructure;

namespace NRLMSIS.Testing
{
    class Testing
    {
        public static void RunTest(string[] args)
        {
            const int NRec = 200;
            int iyd;
            float sec, alt, glat, glong, stl, f107a, f107, apd;
            float[] ap = new float[7];
            string? dummy;

            // ========================================================================
            // Performance Tracking
            // ========================================================================
            var fullRunTimer = Stopwatch.StartNew();
            var initTimer = Stopwatch.StartNew();

            Initialization.Initialize(parmPath: "../", parmFile: "msis21.parm");
            var model = NRLMSIS_Model.CreateDefault();

            initTimer.Stop();
            var initializationTime = initTimer.Elapsed;

            // ========================================================================
            // Calculation Loop with Timing
            // ========================================================================
            var calculationTimer = Stopwatch.StartNew();
            int processedRecords = 0;

            using (StreamReader reader = new("../msis2.1_test_in.txt"))
            using (StreamWriter writer = new("../msis2.1_test_out.txt"))
            {
                dummy = reader.ReadLine();
                writer.WriteLine($"{"iyd",7}{"sec",7}{"alt",7}{"glat",7}{"glong",7}{"stl",7}" +
                                $"{"f107a",7}{"f107",7}{"Ap",7}" +
                                $"{"He",13}{"O",13}{"N2",13}{"O2",13}{"Ar",13}" +
                                $"{"rho",13}{"H",13}{"N",13}{"O*",13}{"NO",13}{"T",8}");

                for (int i = 0; i < NRec; i++)
                {
                    string? line = reader.ReadLine();
                    if (string.IsNullOrWhiteSpace(line))
                        break;

                    string[] parts = line.Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                    iyd = int.Parse(parts[0]);
                    sec = float.Parse(parts[1], CultureInfo.InvariantCulture);
                    alt = float.Parse(parts[2], CultureInfo.InvariantCulture);
                    glat = float.Parse(parts[3], CultureInfo.InvariantCulture);
                    glong = float.Parse(parts[4], CultureInfo.InvariantCulture);
                    stl = float.Parse(parts[5], CultureInfo.InvariantCulture);
                    f107a = float.Parse(parts[6], CultureInfo.InvariantCulture);
                    f107 = float.Parse(parts[7], CultureInfo.InvariantCulture);
                    apd = float.Parse(parts[8], CultureInfo.InvariantCulture);

                    ap[0] = apd;
                    int dayOfYear = iyd % 1000;

                    var location = new AtmosphericLocation(
                        dayOfYear: dayOfYear,
                        universalTimeSeconds: sec,
                        latitudeDegrees: glat,
                        longitudeDegrees: glong,
                        altitudeKm: alt
                    );

                    var weather = new SpaceWeather(
                        f107Average: f107a,
                        f107Daily: f107,
                        apIndices: new double[] { apd, 0, 0, 0, 0, 0, 0 }
                    );

                    var state = model.Calculate(location, weather);
                    var d = state.Densities;

                    // Convert to CGS, but preserve missing values
                    // Missing value is 9.999e-38 and should NOT be scaled
                    double he = ConvertToCGS(d.He, false);
                    double o = ConvertToCGS(d.O, false);
                    double n2 = ConvertToCGS(d.N2, false);
                    double o2 = ConvertToCGS(d.O2, false);
                    double ar = ConvertToCGS(d.Ar, false);
                    double rho = ConvertToCGS(d.TotalMassDensity, true);
                    double h = ConvertToCGS(d.H, false);
                    double n = ConvertToCGS(d.N, false);
                    double ostar = ConvertToCGS(d.AnomalousO, false);
                    double no = ConvertToCGS(d.NO, false);

                    writer.WriteLine($"{iyd,7}{(int)sec,7}{alt,7:F1}{glat,7:F1}{glong,7:F1}{stl,7:F2}" +
                                    $"{f107a,7:F1}{f107,7:F1}{ap[0],7:F1}" +
                                    $"{he,13:E4}{o,13:E4}{n2,13:E4}{o2,13:E4}{ar,13:E4}" +
                                    $"{rho,13:E4}{h,13:E4}{n,13:E4}{ostar,13:E4}{no,13:E4}" +
                                    $"{state.Temperature,8:F2}");

                    processedRecords++;
                }
            }

            calculationTimer.Stop();
            fullRunTimer.Stop();

            // ========================================================================
            // Performance Report
            // ========================================================================
            var calculationTime = calculationTimer.Elapsed;
            var fullRunTime = fullRunTimer.Elapsed;
            var averageTimePerSample = calculationTime.TotalMilliseconds / processedRecords;

            Console.WriteLine("========================================================================");
            Console.WriteLine("                    NRLMSIS 2.1 Performance Report                     ");
            Console.WriteLine("========================================================================");
            Console.WriteLine();
            Console.WriteLine($"Records Processed:           {processedRecords}");
            Console.WriteLine($"Output File:                 msis2.1_test_out.txt");
            Console.WriteLine();
            Console.WriteLine("------------------------------------------------------------------------");
            Console.WriteLine("                         Timing Breakdown                              ");
            Console.WriteLine("------------------------------------------------------------------------");
            Console.WriteLine();
            Console.WriteLine($"Initialization Time:         {FormatTimeSpan(initializationTime)}");
            Console.WriteLine($"  - Model loading            ({initializationTime.TotalMilliseconds:F2} ms)");
            Console.WriteLine();
            Console.WriteLine($"Calculation Time (200 pts):  {FormatTimeSpan(calculationTime)}");
            Console.WriteLine($"  - Pure computation         ({calculationTime.TotalMilliseconds:F2} ms)");
            Console.WriteLine();
            Console.WriteLine($"Full Run Time:               {FormatTimeSpan(fullRunTime)}");
            Console.WriteLine($"  - Init + Calculations      ({fullRunTime.TotalMilliseconds:F2} ms)");
            Console.WriteLine();
            Console.WriteLine("------------------------------------------------------------------------");
            Console.WriteLine("                      Per-Sample Performance                           ");
            Console.WriteLine("------------------------------------------------------------------------");
            Console.WriteLine();
            Console.WriteLine($"Average Time per Sample:     {averageTimePerSample:F3} ms");
            Console.WriteLine($"Throughput:                  {1000.0 / averageTimePerSample:F1} samples/second");
            Console.WriteLine();
            Console.WriteLine("------------------------------------------------------------------------");
            Console.WriteLine("                         Performance Metrics                           ");
            Console.WriteLine("------------------------------------------------------------------------");
            Console.WriteLine();

            // Calculate samples per second for different scenarios
            double throughputWithInit = 1000.0 * processedRecords / fullRunTime.TotalMilliseconds;
            double throughputWithoutInit = 1000.0 * processedRecords / calculationTime.TotalMilliseconds;

            Console.WriteLine($"Batch Throughput (w/ init):  {throughputWithInit:F1} samples/sec");
            Console.WriteLine($"Batch Throughput (no init):  {throughputWithoutInit:F1} samples/sec");
            Console.WriteLine();

            // Estimate for larger batches
            double estimate1000 = calculationTime.TotalMilliseconds * 1000.0 / processedRecords;
            double estimate10000 = calculationTime.TotalMilliseconds * 10000.0 / processedRecords;

            Console.WriteLine("Estimated Time for Larger Batches (excluding init):");
            Console.WriteLine($"  - 1,000 samples:           {estimate1000:F2} ms ({estimate1000/1000.0:F2} seconds)");
            Console.WriteLine($"  - 10,000 samples:          {estimate10000:F2} ms ({estimate10000/1000.0:F2} seconds)");
            Console.WriteLine();
            Console.WriteLine("========================================================================");
            Console.WriteLine();

            // Memory information (optional, requires GC stats)
            Console.WriteLine("Memory Statistics:");
            Console.WriteLine($"  Gen 0 Collections:         {GC.CollectionCount(0)}");
            Console.WriteLine($"  Gen 1 Collections:         {GC.CollectionCount(1)}");
            Console.WriteLine($"  Gen 2 Collections:         {GC.CollectionCount(2)}");
            Console.WriteLine($"  Total Memory:              {GC.GetTotalMemory(false) / 1024.0 / 1024.0:F2} MB");
            Console.WriteLine();
            Console.WriteLine("========================================================================");
        }

        // Helper to format TimeSpan for display
        static string FormatTimeSpan(TimeSpan ts)
        {
            if (ts.TotalSeconds >= 1.0)
                return $"{ts.TotalSeconds:F3} seconds";
            else if (ts.TotalMilliseconds >= 1.0)
                return $"{ts.TotalMilliseconds:F2} ms";
            else
                return $"{ts.TotalMilliseconds * 1000.0:F2} μs";
        }

        // Helper to convert SI to CGS, preserving missing values
        static double ConvertToCGS(double value, bool isMassDensity)
        {
            if (SpeciesDensities.IsMissing(value))
                return Constants.DMissing;  // Return as-is

            return isMassDensity ? value / 1e3 : value / 1e6;
        }
    }
}