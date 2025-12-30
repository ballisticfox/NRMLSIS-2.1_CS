// #######################################################################
// MSIS® (NRL-SOF-014-1) SOFTWARE
// NRLMSIS® empirical atmospheric model software. Use is governed by the
// Open Source Academic Research License Agreement contained in the file
// nrlmsis2.1_license.txt, which is part of this software package. BY
// USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
// CONDITIONS OF THE LICENSE.
// #######################################################################

// ===========================================================================
// NRLMSIS 2.1:
// Neutral atmosphere empirical model from the surface to lower exosphere
// ===========================================================================

// ==================================================================================================
// Testing: Test program for NRLMSIS 2.1
// ==================================================================================================

using System;
using System.IO;
using System.Globalization;

namespace NRLMSIS
{
    class Testing
    {
        static void Main(string[] args)
        {
            const int NRec = 200;
            int iyd, mass = 48; // mass parameter (not used in 2.1, kept for compatibility)
            float sec, alt, glat, glong, stl, f107a, f107, apd;
            float[] ap = new float[7];
            float[] d;
            float[] t;
            string dummy;

            // Initialize model, set the path here.
            Initialization.Initializationialize(parmPath: "../", parmFile: "msis21.parm");

            // Open input and output files, loop through records, and call model
            try
            {
                using (StreamReader reader = new StreamReader("../msis2.1_test_in.txt"))
                using (StreamWriter writer = new StreamWriter("../msis2.1_test_out.txt"))
                {
                    // Read and skip header line
                    dummy = reader.ReadLine();

                    // Write output header
                    writer.WriteLine($"{"iyd",7}{"sec",7}{"alt",7}{"glat",7}{"glong",7}{"stl",7}" +
                                   $"{"f107a",7}{"f107",7}{"Ap",7}" +
                                   $"{"He",13}{"O",13}{"N2",13}{"O2",13}{"Ar",13}" +
                                   $"{"rho",13}{"H",13}{"N",13}{"O*",13}{"NO",13}{"T",8}");

                    // Process 200 records
                    for (int i = 0; i < NRec; i++)
                    {
                        string line = reader.ReadLine();
                        if (string.IsNullOrWhiteSpace(line))
                            break;

                        // Parse input line
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

                        // Set ap array (only first element used in this test)
                        ap[0] = apd;

                        // Call legacy interface (gtd8d)
                        LegacyInterface.Calculate(iyd, sec, alt, glat, glong, stl, f107a, f107, ap, mass, out d, out t);

                        // Write output line
                        // Format: 2 integers (7 chars), 3 floats (7.1), 1 float (7.2), 3 floats (7.1),
                        //         10 exponentials (13.4), 1 float (8.2)
                        // Legacy output: d[0]=He, d[1]=O, d[2]=N2, d[3]=O2, d[4]=Ar, d[5]=rho, d[6]=H, d[7]=N, d[8]=O*, d[9]=NO
                        writer.WriteLine($"{iyd,7}{(int)sec,7}{alt,7:F1}{glat,7:F1}{glong,7:F1}{stl,7:F2}" +
                                       $"{f107a,7:F1}{f107,7:F1}{ap[0],7:F1}" +
                                       $"{d[0],13:E4}{d[1],13:E4}{d[2],13:E4}{d[3],13:E4}{d[4],13:E4}" +
                                       $"{d[5],13:E4}{d[6],13:E4}{d[7],13:E4}{d[8],13:E4}{d[9],13:E4}" +
                                       $"{t[1],8:F2}");
                    }
                }

                Console.WriteLine($"Successfully processed {NRec} records.");
                Console.WriteLine("Output written to msis2.1_test_out.txt");
            }
            catch (FileNotFoundException ex)
            {
                Console.WriteLine($"Error: {ex.Message}");
                Console.WriteLine("Make sure msis2.1_test_in.txt and msis21.parm are in the correct location.");
            }
            // catch (Exception ex)
            // {
            //     Console.WriteLine($"Error: {ex.Message}");
            // }
        }
    }
}