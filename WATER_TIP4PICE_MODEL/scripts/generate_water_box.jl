import Pkg
Pkg.add("Packmol")
using Packmol

# Create a water molecule template
function write_water_template()
    open("water.pdb", "w") do f
        write(f, """ATOM      1  OW  SOL     1       0.000   0.000   0.000
ATOM      2  HW1 SOL     1       0.957   0.000   0.000
ATOM      3  HW2 SOL     1      -0.240   0.927   0.000
END
""")
    end
    return "water.pdb"
end

# Main function to generate water box
function generate_water_box()
    println("Generating water box using Packmol API...")
    
    # Write water template
    water_template = write_water_template()
    
    # Create a new simulation
    sim = Simulation()
    
    # Add water molecules with better tolerance for improved packing
    add_structure!(sim, water_template, 5500, 
                  inside=Box([0.0, 0.0, 0.0], [54.0, 54.0, 54.0]))
    
    # Set parameters for improved packing
    set_tolerance!(sim, 2.2)  # Slightly increased tolerance
    set_seed!(sim, 12345)     # Random seed for reproducibility
    set_maxit!(sim, 30)       # Increase max iterations for better convergence
    
    # Run the simulation
    println("Running Packmol to generate water box...")
    output_file = "water_box.pdb"
    success = run_packmol!(sim, output_file)
    
    if success
        println("Water box successfully generated: $output_file")
    else
        error("Failed to generate water box!")
    end
end

# Run the generation
generate_water_box() 