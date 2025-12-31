using NRLMSIS.Testing;
BenchmarkSpecific.RunBenchmark();
// using NRLMSIS.Core;
// using NRLMSIS.Models;
// using NRLMSIS.Infrastructure;
// using BenchmarkDotNet.Attributes;
// using BenchmarkDotNet.Running;

// // Initialize once at startup (loads parameters from file)
// Initialization.Initialize(parmPath: "../", parmFile: "msis21.parm");

// // Create model instance (reusable, thread-safe with separate contexts)
// var model = NRLMSIS_Model.CreateDefault();

// // Define location
// var location = new AtmosphericLocation(
//     dayOfYear: 172,           // June 21
//     universalTimeSeconds: 43200,  // 12:00 UT
//     latitudeDegrees: 45.0,
//     longitudeDegrees: -75.0,
//     altitudeKm: 400.0
// );

// // Define space weather
// var weather = new SpaceWeather(
//     f107Average: 150.0,
//     f107Daily: 150.0,
//     apIndices: new double[] { 4, 0, 0, 0, 0, 0, 0 }
// );

// // Calculate
// var state = model.Calculate(location, weather);

// // Access results
// Console.WriteLine($"Temperature: {state.Temperature:F2} K");
// Console.WriteLine($"O density: {state.Densities.O:E3} m^-3");
// Console.WriteLine($"Mass density: {state.Densities.TotalMassDensity:E3}");