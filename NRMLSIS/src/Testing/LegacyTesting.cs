using System;
using System.IO;
using System.Globalization;
using NRLMSIS.Core;
using NRLMSIS.Models;
using NRLMSIS.Infrastructure;

namespace NRLMSIS.Testing
{
    class Testing
    {
        public static void Main(string[] args)
        {
            const int NRec = 200;
            int iyd;
            float sec, alt, glat, glong, stl, f107a, f107, apd;
            float[] ap = new float[7];
            string dummy;

            Initialization.Initialize(parmPath: "../", parmFile: "msis21.parm");
            var model = NRLMSIS_Model.CreateDefault();


            using (StreamReader reader = new StreamReader("../msis2.1_test_in.txt"))
            using (StreamWriter writer = new StreamWriter("../msis2.1_test_out.txt"))
            {
                dummy = reader.ReadLine();
                writer.WriteLine($"{"iyd",7}{"sec",7}{"alt",7}{"glat",7}{"glong",7}{"stl",7}" +
                                $"{"f107a",7}{"f107",7}{"Ap",7}" +
                                $"{"He",13}{"O",13}{"N2",13}{"O2",13}{"Ar",13}" +
                                $"{"rho",13}{"H",13}{"N",13}{"O*",13}{"NO",13}{"T",8}");

                for (int i = 0; i < NRec; i++)
                {
                    string line = reader.ReadLine();
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
                }
            }

            Console.WriteLine($"Successfully processed {NRec} records.");
            Console.WriteLine("Output written to msis2.1_test_out.txt");
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