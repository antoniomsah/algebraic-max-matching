import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import re

# Comprehensive sorting
def parse_input(x):
    match = re.match(r'input/(\d+)-random(\d+)', x)
    return (int(match.group(1)), int(match.group(2))) if match else (0, 0)

def generate_plots(csv_path):
    # Create benchmark_plots directory if it doesn't exist
    os.makedirs('benchmark_plots', exist_ok=True)
    
    # Read the CSV file
    df = pd.read_csv(csv_path)

    df['SortValue'] = df['Input'].apply(parse_input)
    df = df.sort_values('SortValue')
    
    # Extract unique vertex counts and test groups
    df['Vertices'] = df['Input'].str.extract('(\d+)-')[0].astype(int)
    
    # Individual plots for each vertex count
    vertex_groups = df.groupby('Vertices')
    for vertices, group in vertex_groups:
        plt.figure(figsize=(10, 6))
        
        unique_algorithms = group['Algorithm'].unique()
        unique_algorithms.sort()
        colors = plt.cm.Set3(np.linspace(0, 1, len(unique_algorithms)))
        bar_width = 0.8 / len(unique_algorithms)
        
        for i, algorithm in enumerate(unique_algorithms):
            algo_data = group[group['Algorithm'] == algorithm]
            bars = plt.bar(
                [x + i * bar_width for x in range(len(algo_data))], 
                algo_data['Average Time (ms)'],
                width=bar_width, 
                label=algorithm,
                color=colors[i]
            )
            
            bar_tops = [b.get_y() + b.get_height() for b in bars]
            bar_centers = [b.get_x() + b.get_width() / 2 for b in bars]
            
            if len(bar_tops) > 1:
                plt.plot(bar_centers, bar_tops, color=colors[i], linestyle='--', linewidth=2)
        
        plt.title(f'Average Time for {vertices} Vertices')
        plt.xlabel('Test Cases')
        plt.ylabel('Average Time (ms)')
        plt.legend(title='Algorithms')
        plt.tight_layout()
        
        # Save individual plot
        output_filename = os.path.join('benchmark_plots', f'perfect_matching_{vertices}_vertices.png')
        plt.savefig(output_filename)
        plt.close()
    
    # Comprehensive single plot
    plt.figure(figsize=(15, 8))
    unique_algorithms = df['Algorithm'].unique()
    unique_algorithms.sort()
    colors = plt.cm.Set3(np.linspace(0, 1, len(unique_algorithms)))

    unique_inputs = df['Input'].unique()
    x_labels = [f"{parse_input(x)[0]} (test{parse_input(x)[1]+1})" for x in unique_inputs]
    
    for i, algorithm in enumerate(unique_algorithms):
        algo_data = df[df['Algorithm'] == algorithm]
        plt.scatter(algo_data['Input'], algo_data['Average Time (ms)'], 
                    label=algorithm, color=colors[i], marker='o')
        
        # Interpolate with a line
        sorted_data = algo_data.sort_values('SortValue')
        plt.plot(sorted_data['Input'], sorted_data['Average Time (ms)'], 
                 color=colors[i], linestyle='--')
    
    plt.title('Perfect Matching Time Across Different Inputs')
    plt.xlabel('Input')
    plt.ylabel('Average Time (ms)')
    plt.legend(title='Algorithms')
    plt.xticks(range(len(unique_inputs)), x_labels, rotation=45, ha='right', fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join('benchmark_plots', 'perfect_matching_plot.png'))
    plt.close()

# Usage
generate_plots('perfect_matching.csv')
