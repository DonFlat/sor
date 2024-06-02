import pandas as pd
import matplotlib.pyplot as plt

output_name = "sor_das6_node8_iter10_large_matrix"

proc_num = 8
# Load the data from the CSV file
df = pd.read_csv(f"./experiment_data/sor_source_data/{output_name}.csv")

# Ensure 'Problem Size' is the x-axis
df.set_index("Matrix size", inplace=True)

# Plotting
plt.figure(figsize=(10, 6))

# Iterate through the columns (excluding 'Problem Size' since it's now the index) and plot
for column in df.columns:
    plt.plot(df.index, df[column], label=column, marker='o')

plt.title('MPI Performance Comparison')
plt.xlabel('matrix size (n x n)')
plt.ylabel('Time (ms)')
plt.legend()
plt.grid(True)

# Show the plot
plt.savefig(f"./experiment_data/sor_graph/{output_name}.png")
