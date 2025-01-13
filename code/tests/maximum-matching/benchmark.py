import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re

def run_makefile():
    """Run the Makefile and check for successful execution"""
    try:
        subprocess.run(['make', 'run_benchmark'], check=True)
        print("Makefile executed successfully")
    except subprocess.CalledProcessError as e:
        print(f"Error executing Makefile: {e}")
        raise
    except FileNotFoundError:
        print("Makefile not found in current directory")
        raise

def clean_filename(filename):
    """Remove input/ prefix and .in suffix from filename"""
    # Remove 'input/' prefix if it exists
    if filename.startswith('input/'):
        filename = filename[6:]  # len('input/') = 6
    
    # Remove '.in' suffix if it exists
    if filename.endswith('.in'):
        filename = filename[:-3]  # len('.in') = 3
    
    return filename

def should_include_file(filename):
    """Check if the input file should be included based on its number"""
    # Extract number from filename using regex
    match = re.search(r'input/(\d+)', filename)
    if match:
        number = int(match.group(1))
        return number >= 50
    return True  # Include files that don't match the pattern

def read_and_process_csv(filename='maximum_matching.csv'):
    """Read and process the CSV file, extracting the summary section"""
    try:
        # Read all content to find where summary starts
        with open(filename, 'r') as f:
            lines = f.readlines()

        df = pd.read_csv(filename)

        # Filter out files with numbers < 50
        df = df[df['Input'].apply(should_include_file)]
        
        # Clean filename
        df['Input File'] = df['Input'].apply(clean_filename)
        
        # Convert average time to milliseconds if needed
        df['Average Time (ms)'] = pd.to_numeric(df['Average Time (ms)'])
        
        return df
    except FileNotFoundError:
        print(f"CSV file {filename} not found")
        raise
    except Exception as e:
        print(f"Error processing CSV file: {e}")
        raise

def create_plot(df):
    """Create a line plot with different colors for each algorithm"""
    
    # Create the figure and axis
    plt.figure(figsize=(12, 6))
    
    # Create line plot for each algorithm
    for algorithm in df['Algorithm'].unique():
        algorithm_data = df[df['Algorithm'] == algorithm]
        plt.plot(algorithm_data['Input'].apply(clean_filename), 
                algorithm_data['Average Time (ms)'],
                marker='o',
                label=algorithm,
                linewidth=2)
    
    # Customize the plot
    plt.title('Algorithm Performance Comparison', fontsize=14, pad=20)
    plt.xlabel('Input File', fontsize=8)
    plt.ylabel('Average Time (ms)', fontsize=12)
    plt.legend(title='Algorithm', title_fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')

    plt.tick_params(axis='x', labelsize=6)
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('benchmark_results.png', dpi=300, bbox_inches='tight')
    print("Plot saved as benchmark_results.png")

def main():
    try:
        # Run the Makefile
        run_makefile()
        
        # Process the CSV and create the plot
        df = read_and_process_csv()
        create_plot(df)
        
    except Exception as e:
        print(f"An error occurred: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    main()
