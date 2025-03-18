import os

def main():
    """Main function to run the visualizations."""
    try:
        print("Starting RMSD stability visualization script...")
        
        # Generate synthetic data
        print("Generating synthetic data...")
        data = generate_synthetic_data()
        
        # Create visualization
        print("Creating visualization...")
        output_path = os.path.join(output_dir, "rmsd_stability_despite_drift.png")
        plot_rmsd_energy_stability(data, output_path)
        
        print(f"RMSD stability visualization created successfully!")
        print(f"Output file saved to {output_path}")
    except Exception as e:
        print(f"Error occurred: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 