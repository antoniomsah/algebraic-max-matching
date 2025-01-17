import subprocess
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re
from pathlib import Path

def run_benchmark():
    """Run the benchmark command and return the output"""
    try:
        result = subprocess.run(['make', 'run_benchmark'], capture_output=True, text=True)
        if result.returncode != 0:
            raise subprocess.CalledProcessError(result.returncode, 'make run_benchmark', result.stdout, result.stderr)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running benchmark: {e}")
        print(f"Error output: {e.stderr}")
        return None

def parse_filename(filename):
    """Parse the filename to extract vertices"""
    pattern = r'input/(\d+)-random\d+'
    match = re.match(pattern, filename)
    if match:
        vertices = int(match.group(1))
        return vertices
    return None

def create_boxplots(df):
    """Create boxplots for each test configuration"""
    # Get unique combinations of vertices and matching numbers
    configs = df.groupby(['vertices']).groups.keys()
    
    for vertices in configs:
        # Filter data for current configuration
        config_data = df[
            (df['vertices'] == vertices)
        ]
        
        # Create boxplot
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=config_data, x='Algorithm', y='Average Time (ms)')
        
        plt.title(f'Algorithm Performance Comparison\nVertices: {vertices}')
        plt.xlabel('Algorithm')
        plt.ylabel('Time (milliseconds)')
        plt.xticks(rotation=45)
        
        # Save plot
        output_dir = Path('benchmark_plots')
        output_dir.mkdir(exist_ok=True)
        plt.savefig(
            output_dir / f'boxplot_v{vertices}.png',
            bbox_inches='tight'
        )
        plt.close()

def main():
    # Run benchmark
    print("Running benchmark...")
    run_benchmark()
    
    # Read CSV data
    try:
        df = pd.read_csv('perfect_matching.csv')
    except FileNotFoundError:
        print("Error: perfect_matching.csv not found")
        return
    
    # Extract information from filenames
    df[['vertices']] = pd.DataFrame(
        [parse_filename(f) for f in df['Input']],  # Updated column name
        index=df.index
    )
    
    # Create boxplots
    print("Creating boxplots...")
    create_boxplots(df)
    print("Boxplots have been saved to the 'benchmark_plots' directory")

if __name__ == "__main__":
    main()
