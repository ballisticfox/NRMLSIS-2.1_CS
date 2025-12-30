"""
Compare MSIS output files and report discrepancies
"""

import sys

def parse_line(line):
    """Parse a data line into numerical values"""
    parts = line.split()
    if len(parts) < 19:
        return None
    
    try:
        # First 9 values (iyd, sec, alt, glat, glong, stl, f107a, f107, Ap)
        metadata = [float(x) for x in parts[:9]]
        # Next 10 values are the densities (He, O, N2, O2, Ar, rho, H, N, O*, NO)
        densities = [float(x) for x in parts[9:19]]
        # Last value is temperature
        temp = float(parts[19])
        return metadata, densities, temp
    except (ValueError, IndexError):
        return None

def compare_files(ref_file, test_file, output_file=None):
    """Compare two MSIS output files"""
    
    print(f"Comparing files:")
    print(f"  Reference: {ref_file}")
    print(f"  Test:      {test_file}")
    if output_file:
        print(f"  Output:    {output_file}")
    print()
    
    output_lines = []

    with open(ref_file, 'r') as f_ref, open(test_file, 'r') as f_test:
        # Skip headers
        ref_header = f_ref.readline()
        test_header = f_test.readline()
        
        # Create output header
        header = f"{'Line':>5} {'iyd':>7} {'sec':>7} {'alt':>7} {'glat':>7} {'glong':>7} " + \
                f"{'He_diff':>13} {'O_diff':>13} {'N2_diff':>13} {'O2_diff':>13} {'Ar_diff':>13} " + \
                f"{'rho_diff':>13} {'H_diff':>13} {'N_diff':>13} {'O*_diff':>13} {'NO_diff':>13} {'T_diff':>10}\n"
        output_lines.append(header)

        line_num = 1
        max_errors = {'He': 0, 'O': 0, 'N2': 0, 'O2': 0, 'Ar': 0, 
                      'rho': 0, 'H': 0, 'N': 0, 'O*': 0, 'NO': 0, 'T': 0}
        max_rel_errors = {'He': 0, 'O': 0, 'N2': 0, 'O2': 0, 'Ar': 0, 
                          'rho': 0, 'H': 0, 'N': 0, 'O*': 0, 'NO': 0, 'T': 0}
        
        species_names = ['He', 'O', 'N2', 'O2', 'Ar', 'rho', 'H', 'N', 'O*', 'NO']
        
        total_lines = 0
        lines_with_discrepancies = 0
        
        for ref_line, test_line in zip(f_ref, f_test):
            ref_data = parse_line(ref_line)
            test_data = parse_line(test_line)
            
            if ref_data is None or test_data is None:
                continue
            
            total_lines += 1
            ref_meta, ref_dens, ref_temp = ref_data
            test_meta, test_dens, test_temp = test_data
            
            has_discrepancy = False
            discrepancies = []
            differences = []
            
            # Compare densities
            for i, species in enumerate(species_names):
                ref_val = ref_dens[i]
                test_val = test_dens[i]
                
                abs_diff = test_val - ref_val  # Signed difference
                
                # Calculate relative error (avoid division by zero)
                if abs(ref_val) > 1e-40:
                    rel_error = abs(abs_diff) / abs(ref_val)
                else:
                    rel_error = 0 if abs(abs_diff) < 1e-40 else float('inf')
                
                # Update max errors
                if abs(abs_diff) > max_errors[species]:
                    max_errors[species] = abs(abs_diff)
                if rel_error > max_rel_errors[species]:
                    max_rel_errors[species] = rel_error
                
                differences.append(abs_diff)

                # Flag significant discrepancies (>0.1% relative error or >1e-10 absolute for small values)
                if rel_error > 0.001 or (abs(ref_val) < 1e-30 and abs(abs_diff) > 1e-10):
                    has_discrepancy = True
                    discrepancies.append(f"  {species:3s}: ref={ref_val:.6e}, test={test_val:.6e}, "
                                       f"diff={abs_diff:.6e}, rel_err={rel_error:.6e}")
            
            # Compare temperature
            temp_diff = test_temp - ref_temp  # Signed difference
            
            temp_rel_error = abs(temp_diff) / ref_temp if ref_temp != 0 else 0

            if abs(temp_diff) > max_errors['T']:
                max_errors['T'] = abs(temp_diff)
            if temp_rel_error > max_rel_errors['T']:
                max_rel_errors['T'] = temp_rel_error
            
            if temp_rel_error > 0.001:
                has_discrepancy = True
                discrepancies.append(f"  T  : ref={ref_temp:.2f}, test={test_temp:.2f}, "
                                   f"diff={temp_diff:.6e}, rel_err={temp_rel_error:.6e}")
            
            # Create output line with differences
            output_line = f"{line_num:>5} {int(ref_meta[0]):>7} {int(ref_meta[1]):>7} " + \
                         f"{ref_meta[2]:>7.1f} {ref_meta[3]:>7.1f} {ref_meta[4]:>7.1f} "
            for diff in differences:
                output_line += f"{diff:>13.6e} "
            output_line += f"{temp_diff:>10.3f}\n"
            output_lines.append(output_line)

            # Print line discrepancies
            if has_discrepancy:
                lines_with_discrepancies += 1
                print(f"Line {line_num}: alt={ref_meta[2]:.1f} km, lat={ref_meta[3]:.1f}°, lon={ref_meta[4]:.1f}°")
                for disc in discrepancies:
                    print(disc)
                print()
            
            line_num += 1
    
    # Write output file if requested
    if output_file:
        with open(output_file, 'w') as f_out:
            f_out.writelines(output_lines)
        print(f"Difference file written to: {output_file}\n")

    # Print summary
    print("="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total lines compared: {total_lines}")
    print(f"Lines with discrepancies (>0.1% rel error): {lines_with_discrepancies}")
    print()
    
    print("Maximum Absolute Errors:")
    for species in species_names + ['T']:
        print(f"  {species:3s}: {max_errors[species]:.6e}")
    print()
    
    print("Maximum Relative Errors:")
    for species in species_names + ['T']:
        if max_rel_errors[species] == float('inf'):
            print(f"  {species:3s}: inf (division by near-zero)")
        else:
            print(f"  {species:3s}: {max_rel_errors[species]:.6e} ({max_rel_errors[species]*100:.4f}%)")


if __name__ == '__main__':
    # if len(sys.argv) != 3:
    #     print("Usage: python compare_msis_output.py <reference_file> <test_file>")
    #     print("Example: python compare_msis_output.py msis2.1_test_ref_dp.txt msis2.1_test_out.txt")
    #     sys.exit(1)
    
    compare_files("msis2.1_test_ref_dp.txt", "NRMLSIS/msis2.1_test_out.txt", "msis2.1_test_compare.txt")